#!/usr/bin/env python3
"""Step 4 of the telescope skill: preflight and submit a range of jobs.

These jobs run for days on a shared allocation. Submit only the range the user
asked for, and only when they have asked. There is no "submit everything else"
mode on purpose.

Preflight aborts the whole batch by default: a range where some manifests are
not ready usually means the previous step is still running, and submitting the
rest anyway wastes an allocation. Pass --allow-partial to submit what is ready.
"""

import argparse
import grp
import math
import os
import re
import shlex
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (add_config_arg, load_config, parse_ids, manifest_path,  # noqa: E402
                     read_manifest_rows, sample_inputs, job_name, queued_job_names, human,
                     STRANDED_MODES)

DIRECTIVES = ("partition", "time", "cpus-per-task", "mem", "account")


def directives_of(path):
    got = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#SBATCH"):
                if line.strip() and not line.startswith("#") and not line.startswith("\n"):
                    break
                continue
            m = re.match(r"#SBATCH\s+--([a-z-]+)=(.*)", line.strip())
            if m and m.group(1) in DIRECTIVES:
                got[m.group(1)] = m.group(2).strip()
    return got


def check_environment(cfg, step):
    """Split into fatal and advisory.

    Only some of the --gtf/--genome/--transcript paths are actually opened by
    the BOWTIE and TELESCOPE workflows; run.py accepts the rest and never reads
    them. Blocking a batch on an unread path would be a false alarm, so those
    are warnings. Everything the job genuinely cannot run without is fatal.
    """
    fatal, warn = [], []
    for rel in ("run.py", "lib"):
        p = os.path.join(cfg["code_dir"], rel)
        if not os.path.exists(p):
            fatal.append(f"code_dir is missing {rel}: {p}")

    herv = os.path.join(cfg["ref_dir"], cfg["herv_gtf"])
    if not os.path.isfile(herv):
        (fatal if step == 2 else warn).append(f"herv_gtf missing: {herv}")

    if step == 1:
        # bowtie2 index is a prefix, not a file: check one of its parts.
        idx = os.path.join(cfg["index_dir"], cfg["bowtie2_idx"])
        if not os.path.isfile(idx + ".1.bt2") and not os.path.isfile(idx + ".1.bt2l"):
            fatal.append(f"bowtie2 index not found at prefix: {idx}")

    for key in ("gtf", "genome", "transcript"):
        p = os.path.join(cfg["ref_dir"], cfg[key])
        if not os.path.isfile(p):
            warn.append(f"--{key} path does not exist (unread by this workflow, "
                        f"but fix the config): {p}")
    return fatal, warn


def bowtie2_pairs(cfg, sample):
    """Read pairs bowtie2 aligned for a sample (first line of its step-1 log), or None."""
    fp = os.path.join(cfg["project"], "output", "BOWTIE", sample, f"{sample}_bowties2.log")
    try:
        with open(fp) as f:
            return int(f.readline().split()[0])
    except (OSError, ValueError, IndexError):
        return None


def counts_done(cfg, sample):
    c = os.path.join(cfg["project"], "output", "TELESCOPE", sample, f"{sample}-TE_counts.tsv")
    return os.path.isfile(c) and os.path.getsize(c) > 0


def plan_step2(cfg, rows):
    """--mem / --cpus-per-task / concurrency for one step-2 job, from its samples.

    Each remaining sample's peak is estimated as base + per_mpairs * M pairs.
    The samples run s2_concurrent at a time, so the job must hold the sum of
    the largest that can coincide. Memory on nextgen comes with cpus
    (max_mem_per_cpu_mb each), so cpus follow memory; if even max_cpus_per_job
    is too little, fewer samples run at once.
    -> (plan dict, problem or None)
    """
    todo = [r["sample"] for r in rows if not counts_done(cfg, r["sample"])]
    pairs = {s: bowtie2_pairs(cfg, s) for s in todo}
    missing = [s for s, p in pairs.items() if p is None]
    if missing:
        return None, (f"no bowtie2 read count for {', '.join(missing[:3])} "
                      f"(output/BOWTIE/<s>/<s>_bowties2.log) -- cannot size step 2")
    est = sorted(((float(cfg["s2_mem_base_gb"]) + float(cfg["s2_mem_gb_per_mpairs"]) * p / 1e6, s)
                  for s, p in pairs.items()), reverse=True)
    per_cpu = int(cfg["max_mem_per_cpu_mb"])
    cap_mb = int(cfg["max_cpus_per_job"]) * per_cpu
    floor_mb = mem_mb(cfg["s2_mem"])
    k = max(1, min(int(cfg["s2_concurrent"]), len(est)))
    while True:
        need_mb = int(math.ceil(sum(g for g, _ in est[:k]) * float(cfg["s2_mem_margin"]) * 1024))
        if need_mb <= cap_mb or k == 1:
            break
        k -= 1
    if need_mb > cap_mb:
        return None, (f"{est[0][1]} alone needs ~{need_mb / 1024:.0f}G, more than "
                      f"{cfg['max_cpus_per_job']} cpus x {per_cpu} MB")
    want_mb = max(need_mb, floor_mb)
    cpus = max(int(cfg["s2_cpus"]), int(math.ceil(want_mb / per_cpu)))
    want_mb = cpus * per_cpu             # use all the memory the cpus pay for
    return {"samples": len(est), "concurrency": k, "cpus": cpus, "mem_mb": want_mb,
            "need_gb": round(need_mb / 1024, 1),
            "largest": "%s %.0fM pairs" % (est[0][1][:8], pairs[est[0][1]] / 1e6) if est else "",
            "est_gb": ";".join("%s:%.1f" % (s, g) for g, s in est)}, None


def mem_mb(spec):
    m = re.match(r"^(\d+)([GgMm]?)$", str(spec))
    if not m:
        return 0
    return int(m.group(1)) * (1024 if m.group(2) in ("G", "g", "") else 1)


def sbatch_command(cfg, script, extra=()):
    """sbatch, run as cfg submit_group (default: the account) via sg.

    A job's $PFSDIR is created under /fs/scratch/<job's Unix group>/ -- by
    default the user's *primary* group, whatever --account says. Submitting
    through `sg <group>` makes the job's group, and so its scratch, match.
    """
    cmd = ["sbatch"] + list(extra) + [script]
    group = cfg.get("submit_group") or cfg["account"]
    try:
        gid = grp.getgrnam(group).gr_gid
    except KeyError:
        sys.exit(f"submit_group {group!r} is not a Unix group on this system "
                 f"(set submit_group in the config)")
    if gid == os.getgid():
        return cmd                        # already the job's group
    if gid not in os.getgroups():
        sys.exit(f"you are not a member of group {group!r}; set submit_group in the config")
    return ["sg", group, "-c", " ".join(shlex.quote(c) for c in cmd)]


def preflight(cfg, step, ids, queued):
    """-> (ready_ids, {id: [problems]}, info)"""
    ready, problems, info = [], {}, {}
    for i in ids:
        probs = []
        script = os.path.join(cfg["project"], f"jobs_step{step}", f"job_{i}.sh")
        mpath = manifest_path(cfg, i)

        if not os.path.isfile(script):
            probs.append(f"job script missing: {script}")
        else:
            r = subprocess.run(["bash", "-n", script], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            if r.returncode != 0:
                probs.append(f"bash -n failed: {r.stderr.strip()[:200]}")

        if not os.path.isfile(mpath):
            probs.append(f"manifest missing: {mpath}")
            problems[i] = probs
            continue

        rows = read_manifest_rows(mpath)
        if not rows:
            probs.append("manifest is empty")

        nbytes = 0
        if step == 1:
            unmoded = 0
            for r in rows:
                s = r["sample"]
                if r["bam"] and (r["fq1"] or r["fq2"]):
                    probs.append(f"{s}: row has both bam and fq1/fq2")
                paths = sample_inputs(r)
                if not paths or (not r["bam"] and len(paths) < 2):
                    probs.append(f"{s}: needs fq1+fq2 or bam")
                for p in paths:
                    if not os.path.isfile(p) or os.path.getsize(p) == 0:
                        probs.append(f"missing {'BAM' if r['bam'] else 'FASTQ'} for {s}: {p}")
                    else:
                        nbytes += os.path.getsize(p)
                if r["stranded_mode"] and r["stranded_mode"] not in STRANDED_MODES:
                    probs.append(f"{s}: stranded_mode {r['stranded_mode']!r} is not one of {', '.join(STRANDED_MODES)}")
                elif not r["stranded_mode"]:
                    unmoded += 1
            if unmoded:
                info[i] = f"{unmoded} sample(s) without stranded_mode (measured in the job)"
        else:
            tmp = os.path.join(cfg["project"], f"tmp_{i}")
            done = 0
            for r in rows:
                s = r["sample"]
                bam = os.path.join(tmp, f"{s}.bam")
                counts = os.path.join(cfg["project"], "output", "TELESCOPE",
                                      s, f"{s}-TE_counts.tsv")
                if os.path.isfile(counts) and os.path.getsize(counts) > 0:
                    done += 1
                    continue
                if not os.path.isfile(bam) or os.path.getsize(bam) == 0:
                    probs.append(f"missing BAM for {s}: {bam}")
                else:
                    nbytes += os.path.getsize(bam)
                    try:
                        with open(os.path.join(tmp, f"{s}.stranded_mode")) as f:
                            mode = f.read().strip()
                    except OSError:
                        mode = ""
                    if mode not in STRANDED_MODES:
                        probs.append(f"{s}: no valid {s}.stranded_mode next to its BAM "
                                     f"(step 1 predates the stranded-mode change? re-run step 1)")
            if done == len(rows):
                probs.append(f"all {done} sample(s) already have counts - nothing to do")
            elif done:
                info[i] = f"{done}/{len(rows)} already quantified"

        jn = job_name(cfg, step, i)
        if jn in queued:
            probs.append(f"already in the queue as {jn}")

        if probs:
            problems[i] = probs
        else:
            ready.append(i)
            info[i] = (info.get(i, "") + f" {human(nbytes)} input").strip()
    return ready, problems, info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_config_arg(ap)
    ap.add_argument("--step", choices=["1", "2"], required=True)
    ap.add_argument("--ids", required=True,
                    help="manifest ids to submit: '233-332', '5', '1,4,7-9'")
    ap.add_argument("--dry-run", action="store_true", help="preflight only")
    ap.add_argument("--allow-partial", action="store_true",
                    help="submit the manifests that pass instead of aborting")
    ap.add_argument("--max", type=int, help="cap the number submitted")
    ap.add_argument("--no-auto-size", action="store_true",
                    help="step 2: use the script's --mem/--cpus instead of sizing from read counts")
    args = ap.parse_args()

    cfg = load_config(args.config)
    step = int(args.step)
    ids = parse_ids(args.ids)
    if not ids:
        sys.exit("--ids selected nothing")

    fatal, warn = check_environment(cfg, step)
    for w in warn:
        print(f"WARNING: {w}", file=sys.stderr)
    if fatal:
        print("ERROR: the run environment is not usable; no job would survive:",
              file=sys.stderr)
        for e in fatal:
            print(f"  {e}", file=sys.stderr)
        return 1

    queued = queued_job_names()
    ready, problems, info = preflight(cfg, step, ids, queued)

    print(f"step {step}, manifests {min(ids)}..{max(ids)} ({len(ids)} requested)")
    print(f"  ready   : {len(ready)}")
    print(f"  blocked : {len(problems)}")

    if problems:
        print("\nblocked:")
        for i in sorted(problems)[:30]:
            for p in problems[i][:3]:
                print(f"  manifest_{i}: {p}")
        if len(problems) > 30:
            print(f"  ... and {len(problems) - 30} more")

    # A batch whose scripts disagree on resources is almost always a half-applied
    # template edit; catching it here is cheaper than after 100 jobs land.
    seen = {}
    for i in ready:
        d = directives_of(os.path.join(cfg["project"], f"jobs_step{step}", f"job_{i}.sh"))
        seen.setdefault(tuple(sorted(d.items())), []).append(i)
    if len(seen) > 1:
        print("\nERROR: the selected scripts do not agree on resources:", file=sys.stderr)
        for k, v in seen.items():
            print(f"  {dict(k)}  <- {len(v)} script(s), e.g. job_{v[0]}.sh", file=sys.stderr)
        print("  regenerate them with make_jobs.py --force before submitting.",
              file=sys.stderr)
        return 1
    if seen:
        print(f"\nresources: {dict(next(iter(seen)))}")

    if problems and not args.allow_partial:
        print("\nnothing submitted. Fix the above, or pass --allow-partial to "
              "submit the ready manifests only.", file=sys.stderr)
        return 1
    if not ready:
        print("\nnothing to submit.", file=sys.stderr)
        return 1

    if args.max and len(ready) > args.max:
        print(f"\n--max {args.max}: submitting the first {args.max} of {len(ready)}")
        ready = ready[:args.max]

    plans = {}
    if step == 2 and cfg["s2_auto_size"] and not args.no_auto_size:
        print("\nstep-2 resources sized per manifest (est. peak = %.1f GB + %.2f GB per M pairs, x%.2f):"
              % (float(cfg["s2_mem_base_gb"]), float(cfg["s2_mem_gb_per_mpairs"]), float(cfg["s2_mem_margin"])))
        print("  %-9s %7s %5s %10s %6s %8s  %s" % ("manifest", "samples", "conc", "need", "cpus", "mem", "largest"))
        unsized = []
        for i in ready:
            plan, why = plan_step2(cfg, read_manifest_rows(manifest_path(cfg, i)))
            if why:
                unsized.append(i)
                print(f"  manifest_{i}: {why}", file=sys.stderr)
                continue
            plans[i] = plan
            print("  %-9s %7d %5d %9.1fG %6d %7.0fG  %s" % (
                f"_{i}", plan["samples"], plan["concurrency"], plan["need_gb"], plan["cpus"],
                plan["mem_mb"] / 1024, plan["largest"]))
        if unsized:
            print("\nnothing submitted: some manifests cannot be sized (or pass --no-auto-size).",
                  file=sys.stderr)
            return 1

    group = cfg.get("submit_group") or cfg["account"]
    print(f"\nsubmitting as group {group} (its scratch holds the jobs' $PFSDIR)")

    if args.dry_run:
        print(f"\n--dry-run: would submit {len(ready)} job(s): "
              f"{ready[:10]}{'...' if len(ready) > 10 else ''}")
        return 0

    ok, bad, jids = 0, 0, []
    for i in ready:
        script = os.path.join(cfg["project"], f"jobs_step{step}", f"job_{i}.sh")
        extra = []
        if i in plans:
            p = plans[i]
            extra = [f"--cpus-per-task={p['cpus']}", f"--mem={p['mem_mb']}M",
                     f"--export=ALL,S2_CONCURRENT_OVERRIDE={p['concurrency']}"]
        r = subprocess.run(sbatch_command(cfg, script, extra),
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        m = re.search(r"Submitted batch job (\d+)", r.stdout or "")
        if r.returncode == 0 and m:
            ok += 1
            jids.append(int(m.group(1)))
            if i in plans:
                record_plan(cfg, i, m.group(1), plans[i])
        else:
            bad += 1
            print(f"FAILED manifest_{i}: {(r.stderr or r.stdout).strip()}", file=sys.stderr)

    print(f"\nsubmitted={ok} failed={bad}"
          + (f" range={min(jids)}..{max(jids)}" if jids else ""))
    if jids:
        # sacct takes a comma list, not a range, and job ids are not guaranteed
        # contiguous -- so keep the exact list for later.
        rec = os.path.join(cfg["project"], f".last_submit_step{step}.txt")
        with open(rec, "w") as f:
            f.write(",".join(str(j) for j in jids) + "\n")
        print(f"job ids recorded in {rec}")
        print("\ntrack with:\n"
              f"  sacct -j $(cat {rec}) -X "
              "-o JobID,JobName%18,State,Elapsed,ExitCode -S $(date +%F) -E now\n"
              f"  squeue -u $USER -h -o '%T' | sort | uniq -c")
    return 0 if bad == 0 else 1


def record_plan(cfg, idx, jobid, plan):
    """Append what was asked for, so estimates can be checked against MaxRSS later."""
    fp = os.path.join(cfg["project"], "logs", "step2_resources.tsv")
    os.makedirs(os.path.dirname(fp), exist_ok=True)
    cols = ["manifest", "job", "samples", "concurrency", "need_gb", "cpus", "mem_mb", "est_gb"]
    new = not os.path.exists(fp)
    with open(fp, "a") as f:
        if new:
            f.write("\t".join(cols) + "\n")
        row = dict(plan, manifest=idx, job=jobid)
        f.write("\t".join(str(row[c]) for c in cols) + "\n")


if __name__ == "__main__":
    sys.exit(main())
