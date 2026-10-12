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
    # input_type "fastq": pairs discovered in samples_dir.
    # input_type "bam": one BAM per row of bam_table (a TSV, e.g. a
    # gdc-download batches/batch_NNN/files.tsv), sample name from bam_sample_col.
    "input_type": "fastq",
    "bam_table": "",
    "bam_sample_col": "id",
    "bam_path_col": "path",
    "bam2fq_exclude_flags": "0xB00",
    "strand_n_reads": 1000000,
    "kallisto_index": "",          # default <project>/qc/kallisto_gencode.idx
    "done_marker": "",             # e.g. /x/results/{sample}.done, touched by step 2
    # submit.py: run sbatch as this Unix group ('' = the account). A job's
    # $PFSDIR is created under /fs/scratch/<job's group>/, not the account.
    "submit_group": "",
    # submit.py --step 2 sizes --mem/--cpus per manifest from bowtie2 read pairs:
    # per-sample peak ~ base + per_mpairs * M pairs; the job needs the sum over
    # the samples that run at once, times margin. Calibrate from the
    # telescope_max_rss_gb column of output/TELESCOPE/logs_<N>.tsv.
    "s2_auto_size": True,
    "s2_mem_base_gb": 4.0,
    "s2_mem_gb_per_mpairs": 0.25,
    "s2_mem_margin": 1.2,
    "max_mem_per_cpu_mb": 4027,    # nextgen MaxMemPerCPU
    "max_cpus_per_job": 120,       # nextgen MaxCPUsPerNode
    "gtf": "gencode.v39.annotation.gtf",
    "genome": "GRCh38.p13.genome.chr.fa",
    "transcript": "gencode.v39.transcripts.fa",
    "herv_gtf": "HG38_HERV_LINE_all_families_telescope_ann.gtf",
    "bowtie2_idx": "gencode.v39_bowtie2/human",
    "intel_module": "intel/2021.10.0",
    "conda_module": "miniconda3/24.1.2-py310",
    "intel_cmplr_root": (
        "/apps/spack/0.21/ascend/linux-rhel9-zen2/intel-oneapi-compilers/"
        "gcc/11.4.1/2023.2.3-nkzjnam/compiler/2023.2.3"),
}

REQUIRED = ("project", "code_dir", "ref_dir", "index_dir")

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
    out = dict(DEFAULTS)
    out.update(cfg)
    need = list(REQUIRED) + (["bam_table"] if out["input_type"] == "bam" else ["samples_dir"])
    missing = [k for k in need if not out.get(k)]
    if missing:
        sys.exit(f"{path} is missing required key(s): {', '.join(missing)}")
    if out["input_type"] not in ("fastq", "bam"):
        sys.exit(f"{path}: input_type must be 'fastq' or 'bam', not {out['input_type']!r}")
    for k in need + ["samples_dir"]:
        if out.get(k):
            out[k] = os.path.abspath(os.path.expanduser(out[k]))
    if not out["kallisto_index"]:
        out["kallisto_index"] = os.path.join(out["project"], "qc", "kallisto_gencode.idx")
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


MANIFEST_COLS = ("sample", "fq1", "fq2", "bam", "stranded_mode")
STRANDED_MODES = ("None", "RF", "FR", "R", "F")


def read_manifest_rows(path):
    """-> [{sample, fq1, fq2, bam, stranded_mode}, ...]. Columns are read by
    header name; a headerless or old 3-column (sample fq1 fq2) file still works."""
    rows = []
    with open(path) as f:
        lines = [ln.rstrip("\n").split("\t") for ln in f]
    if not lines:
        return rows
    header = [c.strip().lower() for c in lines[0]]
    if header[0] == "sample":
        body = lines[1:]
    else:
        header, body = ["sample", "fq1", "fq2"], lines
    for parts in body:
        if not parts or not parts[0].strip():
            continue
        r = dict(zip(header, (p.strip() for p in parts)))
        rows.append({c: r.get(c, "") for c in MANIFEST_COLS})
    return rows


def read_manifest(path):
    """-> [(sample, fq1, fq2), ...] (kept for callers that only need FASTQ columns)."""
    return [(r["sample"], r["fq1"], r["fq2"]) for r in read_manifest_rows(path)]


def sample_inputs(row):
    """Input files of one manifest row: [bam] or [fq1, fq2]."""
    return [row["bam"]] if row.get("bam") else [p for p in (row.get("fq1"), row.get("fq2")) if p]


def discover_bams(cfg):
    """-> ({sample: bam_path}, [problems]) from cfg['bam_table']."""
    import csv
    t = cfg["bam_table"]
    if not os.path.isfile(t):
        sys.exit(f"bam_table does not exist: {t}")
    with open(t, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    sc, pc = cfg["bam_sample_col"], cfg["bam_path_col"]
    if rows and (sc not in rows[0] or pc not in rows[0]):
        sys.exit(f"{t} needs columns {sc!r} and {pc!r} (set bam_sample_col / bam_path_col); "
                 f"it has: {', '.join(rows[0])}")
    out, problems = {}, []
    for r in rows:
        s, p = r[sc].strip(), r[pc].strip()
        if not s:
            continue
        if s in out:
            problems.append(f"duplicate sample {s} in {t}")
            continue
        out[s] = os.path.abspath(p)
    return out, problems


def job_name(cfg, step, idx):
    return f"{cfg['job_prefix']}_s{step}_{idx}"


def queued_job_names(user=None):
    """Set of this user's job names currently in the queue. Empty on failure."""
    cmd = ["squeue", "-h", "-o", "%j"]
    cmd += ["-u", user or os.environ.get("USER", "")]
    try:
        out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=60)
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


# ---------------------------------------------------------------------------
# strandedness table (written by check_strandedness.py, read by make_manifests.py)
# ---------------------------------------------------------------------------

STRAND_COLS = ["sample", "input", "size", "mtime", "status", "method", "layout",
               "frac_rf", "frac_fr", "frac_undetermined", "n_reads", "stranded_mode",
               "telescope_mode", "note", "centers", "n_libraries", "detail"]


def strand_table_path(cfg):
    return os.path.join(cfg["project"], "qc", "strandedness.tsv")


def overrides_path(cfg):
    return os.path.join(cfg["project"], "qc", "strandedness_overrides.tsv")


def load_strand_table(fpath):
    """-> {sample: row dict}"""
    out = {}
    if not os.path.isfile(fpath):
        return out
    with open(fpath) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                out[parts[0]] = dict(zip(header, parts + [""] * (len(header) - len(parts))))
    return out


def write_strand_table(fpath, rows):
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    tmp = fpath + ".tmp"
    with open(tmp, "w") as f:
        f.write("\t".join(STRAND_COLS) + "\n")
        for r in sorted(rows, key=lambda r: r["sample"]):
            f.write("\t".join(str(r.get(c, "")) for c in STRAND_COLS) + "\n")
    os.replace(tmp, fpath)


def input_stamp(paths):
    """(joined paths, total size, max mtime) -- what a cached call is valid for."""
    size = mtime = 0
    for p in paths:
        st = os.stat(p)
        size += st.st_size
        mtime = max(mtime, int(st.st_mtime))
    return ";".join(paths), str(size), str(mtime)


def load_overrides(cfg):
    """qc/strandedness_overrides.tsv: sample<TAB>stranded_mode<TAB>reason (hand-edited)."""
    out = {}
    fp = overrides_path(cfg)
    if not os.path.isfile(fp):
        return out
    with open(fp) as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            if not parts[0] or parts[0].startswith("#") or (i == 0 and parts[0] == "sample"):
                continue
            mode = parts[1].strip() if len(parts) > 1 else ""
            if mode not in STRANDED_MODES:
                sys.exit(f"{fp} line {i + 1}: stranded_mode {mode!r} is not one of {', '.join(STRANDED_MODES)}")
            out[parts[0]] = mode
    return out


def strand_calls(cfg, samples):
    """-> {sample: telescope mode} for samples with a current, OK call (overrides win),
    plus {sample: reason} for the ones without."""
    table = load_strand_table(strand_table_path(cfg))
    ov = load_overrides(cfg)
    calls, missing = {}, {}
    for s, inputs in samples.items():
        if s in ov:
            calls[s] = ov[s]
            continue
        r = table.get(s)
        if not r:
            missing[s] = "not checked"
            continue
        try:
            stamp = input_stamp(inputs)
        except OSError:
            missing[s] = "input missing"
            continue
        if (r["input"], r["size"], r["mtime"]) != stamp:
            missing[s] = "input changed since it was checked"
        elif r["status"] != "OK" or r["telescope_mode"] not in STRANDED_MODES:
            missing[s] = f"check failed: {r.get('detail', '')[:80]}"
        else:
            calls[s] = r["telescope_mode"]
    return calls, missing


def discover_samples(cfg):
    """-> ({sample: [input paths]}, {sample: {'fq1','fq2'} or {'bam'}}, [notes])"""
    if cfg["input_type"] == "bam":
        bams, notes = discover_bams(cfg)
        return ({s: [p] for s, p in bams.items()}, {s: {"bam": p} for s, p in bams.items()}, notes)
    from make_manifests import discover_pairs          # FASTQ discovery lives there
    paired, unpaired = discover_pairs(cfg["samples_dir"], cfg["fastq_suffix_1"], cfg["fastq_suffix_2"])
    notes = [f"{len(unpaired)} unpaired FASTQ sample(s) ignored"] if unpaired else []
    return ({s: [a, b] for s, (a, b) in paired.items()},
            {s: {"fq1": a, "fq2": b} for s, (a, b) in paired.items()}, notes)
