#!/usr/bin/env python3
"""Decide every sample's Telescope --stranded_mode before any job is submitted.

    BAM samples  : RSeQC infer_experiment.py on the BAM's own alignments (~20 s)
    FASTQ samples: kallisto pseudo-alignment of 1M reads + infer_experiment.py
                   (~1-2 min; the kallisto index is built once per project)

Both give the same two fractions and the same call (code/lib/strandedness.py):
RF / FR / None, or UNCLEAR for mixed or weakly stranded libraries, which
Telescope then runs as None. Results land in <project>/qc/strandedness.tsv,
cached by (input path, size, mtime), and make_manifests.py copies them into the
manifests' stranded_mode column.

To force a value for a sample, add it to qc/strandedness_overrides.tsv
(sample<TAB>stranded_mode<TAB>reason); overrides win over measured calls.

Runs with the immuno_telescope env's tools (samtools, infer_experiment.py,
check_strandedness, kallisto). Keep this file Python 3.6 compatible: array
workers run under that env's interpreter.
"""

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (add_config_arg, load_config, discover_samples, input_stamp,  # noqa: E402
                     load_strand_table, write_strand_table, strand_table_path, overrides_path,
                     load_overrides, human, SKILL_DIR, STRAND_COLS)

TOOLS = {"bam": ("samtools", "infer_experiment.py"),
         "fastq": ("check_strandedness", "kallisto", "infer_experiment.py")}


def strand_lib(cfg):
    """code/lib/strandedness.py -- the same module step 1 uses. Prefer the
    project's code_dir so the pre-check and the jobs agree."""
    for d in (os.path.join(cfg["code_dir"], "lib"),
              os.path.join(os.path.dirname(os.path.dirname(SKILL_DIR)), "code", "lib")):
        if os.path.isfile(os.path.join(d, "strandedness.py")):
            sys.path.insert(0, d)
            import strandedness
            return strandedness
    sys.exit("cannot find lib/strandedness.py in code_dir or this repository")


def ref(cfg, key):
    return os.path.join(cfg["ref_dir"], cfg[key])


def check_one(job):
    """job = (sample, kind, [paths], cfg) -> table row"""
    sample, kind, paths, cfg = job
    lib = strand_lib(cfg)
    row = {c: "" for c in STRAND_COLS}
    row["sample"] = sample
    try:
        row["input"], row["size"], row["mtime"] = input_stamp(paths)
        if kind == "bam":
            r = lib.infer_bam(paths[0], os.path.join(cfg["project"], "qc", "genes.protein_coding.bed12"),
                              int(cfg["strand_n_reads"]))
        else:
            r = lib.infer_fastq(paths[0], paths[1], ref(cfg, "gtf"), ref(cfg, "transcript"),
                                cfg["kallisto_index"], int(cfg["strand_n_reads"]),
                                workdir=os.path.join(cfg["project"], "qc"))
        for k in ("method", "layout", "frac_rf", "frac_fr", "frac_undetermined", "n_reads",
                  "stranded_mode", "telescope_mode", "note", "centers", "n_libraries"):
            v = r.get(k, "")
            row[k] = "" if v is None else v
        row["status"] = "OK"
    except Exception as e:  # recorded, not fatal: one bad file must not stop the sweep
        row["status"] = "ERROR"
        row["detail"] = " ".join(str(e).split())[:400]
    return row


def prepare(cfg, kinds):
    """Build what every worker would otherwise race to build."""
    lib = strand_lib(cfg)
    os.makedirs(os.path.join(cfg["project"], "qc"), exist_ok=True)
    if "bam" in kinds:
        lib.bed_from_gtf(ref(cfg, "gtf"), os.path.join(cfg["project"], "qc", "genes.protein_coding.bed12"))
    if "fastq" in kinds:
        if not os.path.isfile(cfg["kallisto_index"]):
            print(f"building kallisto index once: {cfg['kallisto_index']} (~10 min)")
        lib.ensure_kallisto_index(ref(cfg, "transcript"), cfg["kallisto_index"])


ARRAY_TMPL = """#!/bin/bash
#SBATCH --job-name={prefix}_strand
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --array=0-{last}%{throttle}
#SBATCH --output={qc}/logs/strand_%A_%a.out
#SBATCH --error={qc}/logs/strand_%A_%a.err

set -euo pipefail
module load {conda_module}
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate {conda_env}

python {script} --config {config} --worker \\
  --task-id "${{SLURM_ARRAY_TASK_ID}}" --array-chunk {chunk} -j {cpus}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_config_arg(ap)
    ap.add_argument("-j", "--jobs", type=int, default=4, help="parallel checks (default: %(default)s)")
    ap.add_argument("--slurm", action="store_true", help="write a Slurm array script instead of checking here")
    ap.add_argument("--submit", action="store_true", help="with --slurm, actually sbatch the array")
    ap.add_argument("--array-chunk", type=int, default=50, help="samples per array task (default: %(default)s)")
    ap.add_argument("--array-throttle", type=int, default=10)
    ap.add_argument("--array-time", default="04:00:00")
    ap.add_argument("--array-cpus", type=int, default=8)
    ap.add_argument("--array-mem", default="32G")
    ap.add_argument("--sample-only", default="", help="comma-separated sample ids to check")
    ap.add_argument("--recheck", action="store_true", help="ignore cached calls")
    ap.add_argument("--merge", action="store_true", help="fold array results into qc/strandedness.tsv and report")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--task-id", type=int, help=argparse.SUPPRESS)
    args = ap.parse_args()

    cfg = load_config(args.config)
    qc = os.path.join(cfg["project"], "qc")
    table_fp = strand_table_path(cfg)
    pending_fp = os.path.join(qc, "strand_pending.tsv")
    parts = os.path.join(qc, "strand_parts")
    kind = cfg["input_type"]

    # ---- worker: one slice of the pending list -> its own part file ----
    if args.worker:
        with open(pending_fp) as f:
            todo = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
        mine = todo[args.task_id * args.array_chunk:(args.task_id + 1) * args.array_chunk]
        prepare(cfg, {kind})
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            rows = list(ex.map(check_one, [(s, kind, p.split(";"), cfg) for s, p in mine]))
        os.makedirs(parts, exist_ok=True)
        write_strand_table(os.path.join(parts, f"part_{args.task_id}.tsv"), rows)
        print(f"task {args.task_id}: {len(rows)} checked, {sum(r['status'] != 'OK' for r in rows)} failed")
        return 0

    cache = {} if args.recheck else load_strand_table(table_fp)

    if args.merge:
        if not os.path.isdir(parts):
            sys.exit(f"no array results in {parts}")
        merged = dict(cache)
        n = 0
        for fn in sorted(os.listdir(parts)):
            if fn.endswith(".tsv"):
                for s, r in load_strand_table(os.path.join(parts, fn)).items():
                    merged[s] = r
                    n += 1
        write_strand_table(table_fp, merged.values())
        shutil.rmtree(parts)
        print(f"merged {n} result(s) -> {table_fp}")
        return report(cfg, merged)

    samples, _, notes = discover_samples(cfg)
    for n_ in notes:
        print(f"note: {n_}")
    if args.sample_only:
        want = {s.strip() for s in args.sample_only.split(",") if s.strip()}
        samples = {s: p for s, p in samples.items() if s in want}
    if not samples:
        sys.exit("no samples found (check samples_dir / bam_table in the config)")

    pending = []
    for s, paths in sorted(samples.items()):
        c = cache.get(s)
        try:
            stamp = input_stamp(paths)
        except OSError:
            pending.append((s, paths))
            continue
        if c and (c["input"], c["size"], c["mtime"]) == stamp and c["status"] == "OK":
            continue
        pending.append((s, paths))

    nbytes = sum(os.path.getsize(p) for _, ps in pending for p in ps if os.path.exists(p))
    print(f"input type    : {kind}")
    print(f"samples       : {len(samples)}")
    print(f"already called: {len(samples) - len(pending)} (cached)")
    print(f"to check      : {len(pending)}  ({human(nbytes)} of input)")
    if not pending:
        return report(cfg, {s: cache[s] for s in samples if s in cache})
    if args.dry_run:
        print("\n--dry-run: nothing checked")
        return 0

    missing_tools = [t for t in TOOLS[kind] if not shutil.which(t)]
    if missing_tools and not args.slurm:
        sys.exit(f"not on PATH: {', '.join(missing_tools)}\n"
                 f"  conda activate {cfg['conda_env']}  (or use --slurm --submit)")

    os.makedirs(qc, exist_ok=True)
    if args.slurm:
        with open(pending_fp, "w") as f:
            for s, ps in pending:
                f.write(f"{s}\t{';'.join(ps)}\n")
        ntasks = (len(pending) + args.array_chunk - 1) // args.array_chunk
        os.makedirs(os.path.join(qc, "logs"), exist_ok=True)
        script = os.path.join(qc, "strandedness_array.sh")
        with open(script, "w") as f:
            f.write(ARRAY_TMPL.format(
                prefix=cfg["job_prefix"], account=cfg["account"], partition=cfg["partition"],
                time=args.array_time, cpus=args.array_cpus, mem=args.array_mem,
                last=ntasks - 1, throttle=args.array_throttle, qc=qc,
                conda_module=cfg["conda_module"], conda_env=cfg["conda_env"],
                script=os.path.abspath(__file__), config=os.path.abspath(args.config),
                chunk=args.array_chunk))
        os.chmod(script, 0o750)
        print(f"\narray script  : {script}  ({ntasks} task(s) x {args.array_chunk} samples)")
        if kind == "fastq" and not os.path.isfile(cfg["kallisto_index"]):
            print(f"note          : the first task builds {cfg['kallisto_index']}; the others wait for it")
        if not args.submit:
            print("not submitted (pass --submit to sbatch it)")
            return 0
        p = subprocess.Popen(["sbatch", script], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True)
        out, err = p.communicate()
        sys.stdout.write(out)
        sys.stderr.write(err)
        if p.returncode != 0:
            return 1
        print(f"\nwhen the array finishes:\n  python {os.path.abspath(__file__)} --config {args.config} --merge")
        return 0

    prepare(cfg, {kind})
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        rows = list(ex.map(check_one, [(s, kind, ps, cfg) for s, ps in pending]))
    merged = dict(cache)
    for r in rows:
        merged[r["sample"]] = r
    write_strand_table(table_fp, merged.values())
    print(f"\nwrote {table_fp}")
    return report(cfg, {s: merged[s] for s in samples if s in merged})


def report(cfg, rows):
    from collections import Counter
    ok = [r for r in rows.values() if r["status"] == "OK"]
    bad = [r for r in rows.values() if r["status"] != "OK"]
    print(f"\ncalled : {len(ok)}")
    for mode, n in sorted(Counter(r["telescope_mode"] for r in ok).items()):
        print(f"  telescope --stranded_mode {mode:<5} {n}")
    unclear = [r for r in ok if r["stranded_mode"] == "UNCLEAR"]
    if unclear:
        print(f"\nUNCLEAR ({len(unclear)}; will run as None unless overridden in "
              f"{os.path.relpath(overrides_path(cfg), cfg['project'])}):")
        for r in unclear[:20]:
            print(f"  {r['sample']}  RF {r['frac_rf']}  FR {r['frac_fr']}  {r['note']}")
    by_center = Counter((r["centers"], r["telescope_mode"]) for r in ok if r["centers"])
    if by_center:
        print("\nby sequencing center (BAM @RG CN):")
        for (c, m), n in sorted(by_center.items()):
            print(f"  {c:<24} {m:<5} {n}")
    ov = load_overrides(cfg)
    if ov:
        print(f"\noverrides in effect: {len(ov)}")
    if bad:
        print(f"\nFAILED ({len(bad)}):")
        for r in bad[:20]:
            print(f"  {r['sample']}: {r['detail']}")
        return 1
    print("\nnext: make_manifests.py --require-strand")
    return 0


if __name__ == "__main__":
    sys.exit(main())
