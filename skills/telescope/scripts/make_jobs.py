#!/usr/bin/env python3
"""Step 3 of the telescope skill: render SLURM job scripts from the templates.

One script per manifest, in jobs_step1/ (BOWTIE) and jobs_step2/ (TELESCOPE).

Existing scripts are left alone by default. A queued job reads its script at
start time, so rewriting one under a pending job silently changes what runs --
use --force only for manifests that have not been submitted.
"""

import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (add_config_arg, load_config, parse_ids, SKILL_DIR,  # noqa: E402
                     existing_manifest_ids, manifest_path, job_name)

TEMPLATES = {1: "step1.sh.tmpl", 2: "step2.sh.tmpl"}


def subst(cfg, step, idx):
    return {
        "{{ID}}": str(idx),
        "{{JOB_NAME}}": job_name(cfg, step, idx),
        "{{ACCOUNT}}": cfg["account"],
        "{{PARTITION}}": cfg["partition"],
        "{{PROJECT}}": cfg["project"],
        "{{SAMPLES_DIR}}": cfg["samples_dir"],
        "{{CODE_DIR}}": cfg["code_dir"],
        "{{REF_DIR}}": cfg["ref_dir"],
        "{{INDEX_DIR}}": cfg["index_dir"],
        "{{CONDA_MODULE}}": cfg["conda_module"],
        "{{CONDA_ENV}}": cfg["conda_env"],
        "{{INTEL_MODULE}}": cfg["intel_module"],
        "{{INTEL_CMPLR_ROOT}}": cfg["intel_cmplr_root"],
        "{{SEED}}": str(cfg["seed"]),
        "{{GTF}}": cfg["gtf"],
        "{{GENOME}}": cfg["genome"],
        "{{TRANSCRIPT}}": cfg["transcript"],
        "{{HERV_GTF}}": cfg["herv_gtf"],
        "{{BOWTIE2_IDX}}": cfg["bowtie2_idx"],
        "{{TRIMGALORE_N_CORES}}": str(cfg["trimgalore_n_cores"]),
        "{{SAMPLE_CHECK_N_CORES}}": str(cfg["sample_check_n_cores"]),
        "{{S1_CPUS}}": str(cfg["s1_cpus"]),
        "{{S1_MEM}}": cfg["s1_mem"],
        "{{S1_TIME}}": cfg["s1_time"],
        "{{S2_CPUS}}": str(cfg["s2_cpus"]),
        "{{S2_MEM}}": cfg["s2_mem"],
        "{{S2_TIME}}": cfg["s2_time"],
        "{{S2_CONCURRENT}}": str(cfg["s2_concurrent"]),
    }


def render(text, mapping):
    for k, v in mapping.items():
        text = text.replace(k, v)
    left = sorted(set(re.findall(r"\{\{[A-Z0-9_]+\}\}", text)))
    if left:
        sys.exit(f"template has placeholders the config does not fill: {', '.join(left)}")
    return text


def check_mem_ceiling(cfg):
    """nextgen bills extra cpus when --mem exceeds cpus * MaxMemPerCPU."""
    warn = []
    per_cpu_mb = 4027
    for step, c, m in ((1, cfg["s1_cpus"], cfg["s1_mem"]),
                       (2, cfg["s2_cpus"], cfg["s2_mem"])):
        mt = re.match(r"^(\d+)([GgMm])$", str(m))
        if not mt:
            continue
        mb = int(mt.group(1)) * (1024 if mt.group(2) in "Gg" else 1)
        ceiling = int(c) * per_cpu_mb
        if mb > ceiling:
            warn.append(f"  step {step}: --mem={m} with {c} cpus exceeds the "
                        f"{ceiling / 1024:.1f}G ceiling; Slurm will allocate and "
                        f"bill ~{-(-mb // per_cpu_mb)} cpus instead of {c}")
    return warn


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_config_arg(ap)
    ap.add_argument("--step", choices=["1", "2", "both"], default="both")
    ap.add_argument("--ids", help="manifest ids: '5', '233-332', '1,4,7-9' (default: all)")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing scripts (never for submitted manifests)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)

    warnings = check_mem_ceiling(cfg)
    if warnings:
        print("WARNING: requested memory exceeds the per-cpu ceiling",
              file=sys.stderr)
        for w in warnings:
            print(w, file=sys.stderr)

    available = existing_manifest_ids(cfg)
    if not available:
        sys.exit(f"no manifests in {cfg['project']}/manifests\n"
                 f"  run make_manifests.py first.")
    ids = parse_ids(args.ids) if args.ids else available
    unknown = [i for i in ids if i not in set(available)]
    if unknown:
        sys.exit(f"no manifest for id(s): {unknown[:20]}")

    steps = [1, 2] if args.step == "both" else [int(args.step)]
    written, skipped = {}, {}

    for step in steps:
        tpath = os.path.join(SKILL_DIR, "templates", TEMPLATES[step])
        if not os.path.isfile(tpath):
            sys.exit(f"template missing: {tpath}")
        with open(tpath) as f:
            tmpl = f.read()
        outdir = os.path.join(cfg["project"], f"jobs_step{step}")
        if not args.dry_run:
            os.makedirs(outdir, exist_ok=True)
            os.makedirs(os.path.join(cfg["project"], "logs", f"step{step}"), exist_ok=True)
        w, s = [], []
        for i in ids:
            fp = os.path.join(outdir, f"job_{i}.sh")
            if os.path.exists(fp) and not args.force:
                s.append(i)
                continue
            if not os.path.isfile(manifest_path(cfg, i)):
                sys.exit(f"manifest_{i}.tsv missing; refusing to write job_{i}.sh")
            body = render(tmpl, subst(cfg, step, i))
            if not args.dry_run:
                with open(fp, "w") as f:
                    f.write(body)
                os.chmod(fp, 0o750)
            w.append(i)
        written[step], skipped[step] = w, s

    for step in steps:
        w, s = written[step], skipped[step]
        tag = "would write" if args.dry_run else "wrote"
        print(f"step{step}: {tag} {len(w)} script(s) -> "
              f"{cfg['project']}/jobs_step{step}")
        if s:
            print(f"       skipped {len(s)} that already exist "
                  f"(pass --force to overwrite)")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return 0

    # A script that does not parse is worse than no script: it fails minutes
    # into the allocation with a syntax error and no useful log.
    bad = []
    for step in steps:
        for i in written[step]:
            fp = os.path.join(cfg["project"], f"jobs_step{step}", f"job_{i}.sh")
            r = subprocess.run(["bash", "-n", fp], capture_output=True, text=True)
            if r.returncode != 0:
                bad.append((fp, r.stderr.strip()))
    if bad:
        print(f"\n{len(bad)} script(s) FAILED bash -n:", file=sys.stderr)
        for fp, err in bad[:10]:
            print(f"  {fp}: {err}", file=sys.stderr)
        return 1
    total = sum(len(written[s]) for s in steps)
    if total:
        print(f"\nall {total} new script(s) pass bash -n")
    print(f"\nnext: submit.py --config {args.config} --step <1|2> --ids <range>")
    return 0


if __name__ == "__main__":
    sys.exit(main())
