#!/usr/bin/env python3
"""Step 1 of the telescope skill: verify every FASTQ is a complete gzip stream.

`gzip -t` decompresses the whole file, so a full cohort is tens of TB of I/O.
Run it as a Slurm array (--slurm --submit) for anything larger than a pilot.

Results land in <project>/qc/fastq_check.tsv and are cached by (path, size,
mtime): re-running only checks files that are new or have changed.
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import add_config_arg, load_config, human, SKILL_DIR  # noqa: E402

HEADER = ["path", "size", "mtime", "status", "detail"]


def discover(samples_dir, suffixes):
    """All FASTQ under samples_dir, skipping azcopy's stranded temp files."""
    found = []
    for fn in sorted(os.listdir(samples_dir)):
        # .azDownload-<guid>-<name> files are partial or duplicate transfers.
        if fn.startswith("."):
            continue
        if not any(fn.endswith(s) for s in suffixes):
            continue
        p = os.path.join(samples_dir, fn)
        if os.path.isfile(p):
            found.append(p)
    return found


def sample_of(path, suffixes):
    base = os.path.basename(path)
    for s in suffixes:
        if base.endswith(s):
            return base[:-len(s)]
    return base


def load_cache(fpath):
    """-> {path: (size, mtime, status, detail)}"""
    cache = {}
    if not os.path.isfile(fpath):
        return cache
    with open(fpath) as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[0] == "path":
                continue
            if len(parts) < 4:
                continue
            path, size, mtime, status = parts[0], parts[1], parts[2], parts[3]
            detail = parts[4] if len(parts) > 4 else ""
            cache[path] = (size, mtime, status, detail)
    return cache


def write_rows(fpath, rows):
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    tmp = fpath + ".tmp"
    with open(tmp, "w") as f:
        f.write("\t".join(HEADER) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp, fpath)


def check_one(path):
    """-> (path, size, mtime, status, detail)"""
    try:
        st = os.stat(path)
    except OSError as e:
        return (path, 0, 0, "BAD", f"stat failed: {e}")
    if st.st_size == 0:
        return (path, 0, int(st.st_mtime), "BAD", "empty file")
    try:
        p = subprocess.run(["gzip", "-t", path], capture_output=True, text=True)
    except OSError as e:
        return (path, st.st_size, int(st.st_mtime), "ERROR", f"could not run gzip: {e}")
    if p.returncode == 0:
        return (path, st.st_size, int(st.st_mtime), "OK", "")
    detail = " ".join((p.stderr or "").split())[:400] or f"gzip rc={p.returncode}"
    return (path, st.st_size, int(st.st_mtime), "BAD", detail)


ARRAY_TMPL = """#!/bin/bash
#SBATCH --job-name={prefix}_fqchk
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --array=0-{last}%{throttle}
#SBATCH --output={qc}/logs/fqchk_%A_%a.out
#SBATCH --error={qc}/logs/fqchk_%A_%a.err

set -euo pipefail
module load {conda_module}
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate {conda_env}

python {script} --config {config} --worker \\
  --filelist {qc}/pending.txt \\
  --task-id "${{SLURM_ARRAY_TASK_ID}}" \\
  --array-chunk {chunk} \\
  --part-dir {qc}/parts \\
  -j {cpus}
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_config_arg(ap)
    ap.add_argument("-j", "--jobs", type=int, default=4,
                    help="parallel gzip -t processes (default: %(default)s)")
    ap.add_argument("--slurm", action="store_true",
                    help="write a Slurm array script instead of checking locally")
    ap.add_argument("--submit", action="store_true",
                    help="with --slurm, actually sbatch the array")
    ap.add_argument("--array-chunk", type=int, default=20,
                    help="files per array task (default: %(default)s)")
    ap.add_argument("--array-throttle", type=int, default=20,
                    help="max concurrent array tasks (default: %(default)s)")
    ap.add_argument("--array-time", default="04:00:00")
    ap.add_argument("--array-cpus", type=int, default=4)
    ap.add_argument("--array-mem", default="8G")
    ap.add_argument("--sample-only", default="",
                    help="comma-separated sample ids to check")
    ap.add_argument("--recheck", action="store_true",
                    help="ignore the cache and re-test everything")
    ap.add_argument("--merge", action="store_true",
                    help="fold qc/parts/*.tsv into qc/fastq_check.tsv and report")
    # worker-side
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--filelist", help=argparse.SUPPRESS)
    ap.add_argument("--task-id", type=int, help=argparse.SUPPRESS)
    ap.add_argument("--part-dir", help=argparse.SUPPRESS)
    args = ap.parse_args()

    cfg = load_config(args.config)
    qc = os.path.join(cfg["project"], "qc")
    out_fpath = os.path.join(qc, "fastq_check.tsv")
    suffixes = (cfg["fastq_suffix_1"], cfg["fastq_suffix_2"])

    # ---- worker: check one slice of pending.txt, write its own part file ----
    if args.worker:
        with open(args.filelist) as f:
            allf = [ln.strip() for ln in f if ln.strip()]
        lo = args.task_id * args.array_chunk
        mine = allf[lo:lo + args.array_chunk]
        if not mine:
            print(f"task {args.task_id}: nothing to do")
            return 0
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            rows = list(ex.map(check_one, mine))
        os.makedirs(args.part_dir, exist_ok=True)
        write_rows(os.path.join(args.part_dir, f"part_{args.task_id}.tsv"), rows)
        bad = sum(1 for r in rows if r[3] != "OK")
        print(f"task {args.task_id}: {len(rows)} checked, {bad} bad")
        return 0

    cache = {} if args.recheck else load_cache(out_fpath)

    # ---- merge mode: fold worker parts back into the main table ----
    if args.merge:
        pdir = os.path.join(qc, "parts")
        if not os.path.isdir(pdir):
            sys.exit(f"no parts directory at {pdir}")
        merged = dict(cache)
        n = 0
        for fn in sorted(os.listdir(pdir)):
            if not fn.endswith(".tsv"):
                continue
            for path, vals in load_cache(os.path.join(pdir, fn)).items():
                merged[path] = vals
                n += 1
        rows = [(p,) + v for p, v in sorted(merged.items())]
        write_rows(out_fpath, rows)
        print(f"merged {n} rows from {pdir} -> {out_fpath}")
        return report(rows, suffixes)

    if not os.path.isdir(cfg["samples_dir"]):
        sys.exit(f"samples_dir does not exist: {cfg['samples_dir']}")
    files = discover(cfg["samples_dir"], suffixes)
    if not files:
        sys.exit(f"no {'/'.join(suffixes)} files found in {cfg['samples_dir']}")

    if args.sample_only:
        want = {s.strip() for s in args.sample_only.split(",") if s.strip()}
        files = [p for p in files if sample_of(p, suffixes) in want]
        if not files:
            sys.exit(f"none of the named samples matched files in {cfg['samples_dir']}")

    pending = []
    for p in files:
        try:
            st = os.stat(p)
        except OSError:
            pending.append(p)
            continue
        c = cache.get(p)
        if c and c[0] == str(st.st_size) and c[1] == str(int(st.st_mtime)) and c[2] == "OK":
            continue
        pending.append(p)

    total = sum(os.path.getsize(p) for p in pending if os.path.exists(p))
    print(f"files found   : {len(files)}")
    print(f"already OK    : {len(files) - len(pending)} (cached)")
    print(f"to check      : {len(pending)}  ({human(total)} to decompress)")
    if not pending:
        return report([(p,) + cache[p] for p in files if p in cache], suffixes)

    os.makedirs(qc, exist_ok=True)

    if args.slurm:
        with open(os.path.join(qc, "pending.txt"), "w") as f:
            f.write("\n".join(pending) + "\n")
        ntasks = (len(pending) + args.array_chunk - 1) // args.array_chunk
        os.makedirs(os.path.join(qc, "logs"), exist_ok=True)
        script = os.path.join(qc, "fastq_check_array.sh")
        with open(script, "w") as f:
            f.write(ARRAY_TMPL.format(
                prefix=cfg["job_prefix"], account=cfg["account"],
                partition=cfg["partition"], time=args.array_time,
                cpus=args.array_cpus, mem=args.array_mem,
                last=ntasks - 1, throttle=args.array_throttle,
                qc=qc, conda_module=cfg["conda_module"], conda_env=cfg["conda_env"],
                script=os.path.abspath(__file__), config=os.path.abspath(args.config),
                chunk=args.array_chunk))
        os.chmod(script, 0o750)
        print(f"\narray script  : {script}  ({ntasks} tasks x {args.array_chunk} files)")
        if not args.submit:
            print("not submitted (pass --submit to sbatch it)")
            return 0
        r = subprocess.run(["sbatch", script], capture_output=True, text=True)
        sys.stdout.write(r.stdout)
        sys.stderr.write(r.stderr)
        if r.returncode != 0:
            return 1
        print(f"\nwhen the array finishes:\n"
              f"  python {os.path.abspath(__file__)} --config {args.config} --merge")
        return 0

    if total > 2 * 1024 ** 4:
        print(f"\nWARNING: {human(total)} of decompression on this node. "
              f"Consider --slurm --submit instead.", file=sys.stderr)

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        rows = list(ex.map(check_one, pending))

    merged = dict(cache)
    for r in rows:
        merged[r[0]] = tuple(str(x) for x in r[1:])
    allrows = [(p,) + v for p, v in sorted(merged.items())]
    write_rows(out_fpath, allrows)
    print(f"\nwrote {out_fpath}")
    return report(allrows, suffixes)


def report(rows, suffixes):
    bad = [r for r in rows if r[3] != "OK"]
    print(f"\nOK  : {len(rows) - len(bad)}")
    print(f"BAD : {len(bad)}")
    if bad:
        print("\nthese files are corrupt and must be re-downloaded:")
        for r in sorted(bad):
            print(f"  {os.path.basename(r[0])}\t{r[3]}\t{r[4]}")
        samples = sorted({sample_of(r[0], suffixes) for r in bad})
        print(f"\naffected samples ({len(samples)}): {','.join(samples)}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
