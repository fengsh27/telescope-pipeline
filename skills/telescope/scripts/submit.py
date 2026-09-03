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
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (add_config_arg, load_config, parse_ids, manifest_path,  # noqa: E402
                     read_manifest, job_name, queued_job_names, human)

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
            r = subprocess.run(["bash", "-n", script], capture_output=True, text=True)
            if r.returncode != 0:
                probs.append(f"bash -n failed: {r.stderr.strip()[:200]}")

        if not os.path.isfile(mpath):
            probs.append(f"manifest missing: {mpath}")
            problems[i] = probs
            continue

        rows = read_manifest(mpath)
        if not rows:
            probs.append("manifest is empty")

        nbytes = 0
        if step == 1:
            for s, f1, f2 in rows:
                for p in (f1, f2):
                    if not p or not os.path.isfile(p) or os.path.getsize(p) == 0:
                        probs.append(f"missing FASTQ for {s}: {p or '(no path)'}")
                    else:
                        nbytes += os.path.getsize(p)
        else:
            tmp = os.path.join(cfg["project"], f"tmp_{i}")
            done = 0
            for s, _, _ in rows:
                bam = os.path.join(tmp, f"{s}.bam")
                counts = os.path.join(cfg["project"], "output", "TELESCOPE",
                                      s, f"{s}-TE_counts.tsv")
                if os.path.isfile(counts) and os.path.getsize(counts) > 0:
                    done += 1
                    continue
                if not os.path.isfile(bam) or os.path.getsize(bam) == 0:
                    probs.append(f"missing BAM for {s}: {bam}")
                elif not os.path.isfile(bam + ".bai"):
                    probs.append(f"missing index: {bam}.bai")
                else:
                    nbytes += os.path.getsize(bam)
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

    if args.dry_run:
        print(f"\n--dry-run: would submit {len(ready)} job(s): "
              f"{ready[:10]}{'...' if len(ready) > 10 else ''}")
        return 0

    ok, bad, jids = 0, 0, []
    for i in ready:
        script = os.path.join(cfg["project"], f"jobs_step{step}", f"job_{i}.sh")
        r = subprocess.run(["sbatch", script], capture_output=True, text=True)
        m = re.search(r"Submitted batch job (\d+)", r.stdout or "")
        if r.returncode == 0 and m:
            ok += 1
            jids.append(int(m.group(1)))
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


if __name__ == "__main__":
    sys.exit(main())
