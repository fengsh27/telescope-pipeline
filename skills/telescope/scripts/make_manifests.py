#!/usr/bin/env python3
"""Step 2 of the telescope skill: group samples into N-sample manifests.

Samples are paired FASTQ discovered in samples_dir (input_type "fastq") or one
BAM per row of bam_table (input_type "bam"). Each manifest row carries the
sample's Telescope --stranded_mode from qc/strandedness.tsv (check_strandedness.py)
or qc/strandedness_overrides.tsv; a blank mode is measured inside step 1.

Manifest indices are permanent: job scripts, tmp_<N> directories, log file
names and logs_<N>.tsv are all keyed to them. Renumbering divorces finished
work from its manifest, so existing manifests are never touched unless you
pass --append (add new ones after the current maximum) or --force (rewrite).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (add_config_arg, load_config, existing_manifest_ids,  # noqa: E402
                     manifest_path, read_manifest, human, discover_samples, strand_calls,
                     MANIFEST_COLS)


def discover_pairs(samples_dir, s1, s2):
    """-> ({sample: (fq1, fq2)}, [unpaired sample ids])"""
    mates = {}
    for fn in sorted(os.listdir(samples_dir)):
        if fn.startswith("."):          # azcopy .azDownload temp files
            continue
        p = os.path.join(samples_dir, fn)
        if not os.path.isfile(p):
            continue
        if fn.endswith(s1):
            mates.setdefault(fn[:-len(s1)], {})["1"] = p
        elif fn.endswith(s2):
            mates.setdefault(fn[:-len(s2)], {})["2"] = p
    paired, unpaired = {}, []
    for s in sorted(mates):
        if "1" in mates[s] and "2" in mates[s]:
            paired[s] = (mates[s]["1"], mates[s]["2"])
        else:
            unpaired.append(s)
    return paired, unpaired


def passed_qc(project):
    """-> set of file paths marked OK in qc/fastq_check.tsv, or None if absent."""
    fpath = os.path.join(project, "qc", "fastq_check.tsv")
    if not os.path.isfile(fpath):
        return None
    ok = set()
    with open(fpath) as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[0] == "path":
                continue
            if len(parts) >= 4 and parts[3] == "OK":
                ok.add(parts[0])
    return ok


def already_in_manifests(cfg, ids):
    seen = {}
    for i in ids:
        for s, _, _ in read_manifest(manifest_path(cfg, i)):
            seen[s] = i
    return seen


def write_manifest(cfg, idx, chunk, rows, modes):
    path = manifest_path(cfg, idx)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(MANIFEST_COLS) + "\n")
        for s in chunk:
            r = dict(rows[s], sample=s, stranded_mode=modes.get(s, ""))
            f.write("\t".join(r.get(c, "") for c in MANIFEST_COLS) + "\n")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_config_arg(ap)
    ap.add_argument("--append", action="store_true",
                    help="keep existing manifests, add only new samples after the max index")
    ap.add_argument("--force", action="store_true",
                    help="delete and rewrite every manifest (only if nothing was submitted)")
    ap.add_argument("--require-qc", action="store_true",
                    help="include only samples whose FASTQ passed check_fastq.py")
    ap.add_argument("--require-strand", action="store_true",
                    help="include only samples with a stranded_mode from check_strandedness.py "
                         "(or an override); without it, missing modes are measured in step 1")
    ap.add_argument("--per-manifest", type=int,
                    help="override samples_per_manifest from the config")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    n_per = args.per_manifest or int(cfg["samples_per_manifest"])
    if n_per < 1:
        sys.exit("--per-manifest must be >= 1")

    unpaired = []
    if cfg["input_type"] == "bam":
        inputs, rows, notes = discover_samples(cfg)
        for n_ in notes:
            print(f"note: {n_}")
        gone = [s for s, ps in inputs.items() if not os.path.isfile(ps[0])]
        print(f"BAM samples    : {len(inputs)}  (from {cfg['bam_table']})")
        if gone:
            print(f"BAM missing    : {len(gone)}  {gone[:10]}{'...' if len(gone) > 10 else ''}")
            for s in gone:
                del inputs[s], rows[s]
        paired = {s: tuple(ps) for s, ps in inputs.items()}      # used for sizes below
        if args.require_qc:
            print("note: --require-qc applies to FASTQ; BAM integrity is the download's md5 check")
            args.require_qc = False
    else:
        if not os.path.isdir(cfg["samples_dir"]):
            sys.exit(f"samples_dir does not exist: {cfg['samples_dir']}")
        paired, unpaired = discover_pairs(
            cfg["samples_dir"], cfg["fastq_suffix_1"], cfg["fastq_suffix_2"])
        rows = {s: {"fq1": a, "fq2": b} for s, (a, b) in paired.items()}
        print(f"complete pairs : {len(paired)}")
        print(f"unpaired       : {len(unpaired)}"
              + (f"  {unpaired[:10]}{'...' if len(unpaired) > 10 else ''}" if unpaired else ""))

    if args.require_qc:
        ok = passed_qc(cfg["project"])
        if ok is None:
            sys.exit("--require-qc given but qc/fastq_check.tsv does not exist.\n"
                     "  run check_fastq.py first.")
        dropped = [s for s, (a, b) in paired.items() if a not in ok or b not in ok]
        for s in dropped:
            del paired[s]
        print(f"dropped by QC  : {len(dropped)}"
              + (f"  {dropped[:10]}{'...' if len(dropped) > 10 else ''}" if dropped else ""))

    modes, no_mode = strand_calls(cfg, {s: list(ps) for s, ps in paired.items()})
    print(f"stranded_mode  : {len(modes)} known"
          + (f" ({', '.join(f'{m} {sum(1 for v in modes.values() if v == m)}' for m in sorted(set(modes.values())))})"
             if modes else "")
          + (f", {len(no_mode)} without" if no_mode else ""))
    if args.require_strand and no_mode:
        for s in no_mode:
            del paired[s]
        reasons = {}
        for s, why in no_mode.items():
            reasons.setdefault(why.split(":")[0], []).append(s)
        for why, ss in reasons.items():
            print(f"dropped        : {len(ss)} ({why})  {ss[:5]}{'...' if len(ss) > 5 else ''}")
        print("                 run check_strandedness.py first, or add them to qc/strandedness_overrides.tsv")
    elif no_mode:
        print("note           : samples without a mode get it measured inside step 1")

    if not paired:
        sys.exit("no samples to write")

    existing = existing_manifest_ids(cfg)
    mdir = os.path.join(cfg["project"], "manifests")

    if existing and not (args.append or args.force):
        seen = already_in_manifests(cfg, existing)
        new = sorted(set(paired) - set(seen))
        print(f"\n{len(existing)} manifest(s) already exist "
              f"(manifest_{min(existing)}..manifest_{max(existing)}), "
              f"covering {len(seen)} sample(s).")
        print(f"samples not yet in any manifest: {len(new)}")
        if new:
            print(f"  {new[:20]}{'...' if len(new) > 20 else ''}")
            print("\nrerun with --append to add these as new manifests, "
                  "or --force to renumber everything from scratch.")
        else:
            print("\nnothing to do. Every discovered sample is already in a manifest.")
        return 0

    if args.force:
        seen, start = {}, 0
        if not args.dry_run:
            for fn in os.listdir(mdir) if os.path.isdir(mdir) else []:
                if fn.startswith("manifest_") and fn.endswith(".tsv"):
                    os.remove(os.path.join(mdir, fn))
        order = sorted(paired)
    else:
        seen = already_in_manifests(cfg, existing) if existing else {}
        start = (max(existing) + 1) if existing else 0
        order = sorted(set(paired) - set(seen))
        if not order:
            print("\nnothing to append: every discovered sample is already in a manifest.")
            return 0

    chunks = [order[i:i + n_per] for i in range(0, len(order), n_per)]
    total_bytes = sum(os.path.getsize(p) for s in order for p in paired[s]
                      if os.path.exists(p))
    kind = "BAM" if cfg["input_type"] == "bam" else "FASTQ"
    print(f"\nwriting {len(chunks)} manifest(s): "
          f"manifest_{start}..manifest_{start + len(chunks) - 1}")
    print(f"samples        : {len(order)}  ({human(total_bytes)} of {kind})")
    if chunks and len(chunks[-1]) < n_per:
        print(f"note           : manifest_{start + len(chunks) - 1} "
              f"has only {len(chunks[-1])} sample(s)")

    if args.dry_run:
        for j, chunk in enumerate(chunks[:5]):
            print(f"  manifest_{start + j}: {', '.join(chunk)}")
        if len(chunks) > 5:
            print(f"  ... {len(chunks) - 5} more")
        print("\n--dry-run: nothing written")
        return 0

    for j, chunk in enumerate(chunks):
        write_manifest(cfg, start + j, chunk, rows, modes)

    up = os.path.join(mdir, "UNPAIRED.txt")
    with open(up, "w") as f:
        f.write("\n".join(unpaired) + ("\n" if unpaired else ""))
    print(f"\nwrote {len(chunks)} manifest(s) to {mdir}")
    print(f"unpaired ids   -> {up}")
    print(f"\nnext: make_jobs.py --config {args.config} --step both")
    return 0


if __name__ == "__main__":
    sys.exit(main())
