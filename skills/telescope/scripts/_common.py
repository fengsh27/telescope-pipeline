"""Shared helpers for the telescope skill scripts."""

import argparse
import json
import os
import re
import subprocess
import sys

DEFAULTS = {
    "partition": "nextgen",
    "account": "PDE0005",
    "job_prefix": "telescope",
    "conda_env": "immuno_telescope",
    "seed": 123456,
    "samples_per_manifest": 4,
    "s1_cpus": 24, "s1_mem": "64G", "s1_time": "96:00:00",
    "s2_cpus": 12, "s2_mem": "47G", "s2_time": "144:00:00",
    "s2_concurrent": 4,
    "trimgalore_n_cores": 6,
    "sample_check_n_cores": 6,
    "fastq_suffix_1": "_1.fastq.gz",
    "fastq_suffix_2": "_2.fastq.gz",
    "gtf": "gencode.v39.annotation.gtf",
    "genome": "GRCh38.p13.genome.fa",
    "transcript": "gencode.v39.transcripts.fa",
    "herv_gtf": "HG38_HERV_LINE_all_families_telescope_ann.gtf",
    "bowtie2_idx": "gencode.v39_bowtie2/human",
    "intel_module": "intel/2021.10.0",
    "conda_module": "miniconda3/24.1.2-py310",
    "intel_cmplr_root": (
        "/apps/spack/0.21/ascend/linux-rhel9-zen2/intel-oneapi-compilers/"
        "gcc/11.4.1/2023.2.3-nkzjnam/compiler/2023.2.3"),
}

REQUIRED = ("project", "samples_dir", "code_dir", "ref_dir", "index_dir")

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def add_config_arg(ap):
    ap.add_argument("--config", default="telescope.config.json",
                    help="path to the project's telescope config (default: %(default)s)")


def load_config(path):
    if not os.path.isfile(path):
        sys.exit(
            f"config not found: {path}\n"
            f"  cp {os.path.join(SKILL_DIR, 'config.example.json')} {path}\n"
            f"  then edit the paths in it.")
    with open(path) as f:
        try:
            cfg = json.load(f)
        except json.JSONDecodeError as e:
            sys.exit(f"{path} is not valid JSON: {e}")
    missing = [k for k in REQUIRED if not cfg.get(k)]
    if missing:
        sys.exit(f"{path} is missing required key(s): {', '.join(missing)}")
    out = dict(DEFAULTS)
    out.update(cfg)
    for k in REQUIRED:
        out[k] = os.path.abspath(os.path.expanduser(out[k]))
    return out


def parse_ids(spec):
    """'3', '5-9', '1,4,7-9' -> sorted unique list of ints."""
    ids = set()
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part.lstrip("-"):
            lo, _, hi = part.partition("-")
            lo, hi = int(lo), int(hi)
            if hi < lo:
                sys.exit(f"bad range '{part}': end is before start")
            ids.update(range(lo, hi + 1))
        else:
            ids.add(int(part))
    return sorted(ids)


def manifest_path(cfg, idx):
    return os.path.join(cfg["project"], "manifests", f"manifest_{idx}.tsv")


def existing_manifest_ids(cfg):
    d = os.path.join(cfg["project"], "manifests")
    if not os.path.isdir(d):
        return []
    rx = re.compile(r"^manifest_(\d+)\.tsv$")
    return sorted(int(m.group(1)) for m in
                  (rx.match(fn) for fn in os.listdir(d)) if m)


def read_manifest(path):
    """-> [(sample, fq1, fq2), ...]; tolerates a missing fq2 column."""
    rows = []
    with open(path) as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[0] == "sample":
                continue
            if not parts[0].strip():
                continue
            rows.append(tuple((parts + ["", ""])[:3]))
    return rows


def job_name(cfg, step, idx):
    return f"{cfg['job_prefix']}_s{step}_{idx}"


def queued_job_names(user=None):
    """Set of this user's job names currently in the queue. Empty on failure."""
    cmd = ["squeue", "-h", "-o", "%j"]
    cmd += ["-u", user or os.environ.get("USER", "")]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return set()
    if out.returncode != 0:
        return set()
    return {ln.strip() for ln in out.stdout.splitlines() if ln.strip()}


def human(n):
    for unit in ("B", "K", "M", "G", "T", "P"):
        if abs(n) < 1024 or unit == "P":
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024.0
