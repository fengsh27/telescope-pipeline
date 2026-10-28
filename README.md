# telescope-pipeline

Quantifies HERV and LINE retrotransposon expression from paired-end RNA-seq
with [Telescope](https://github.com/mlbendall/telescope), on the OSC Slurm
cluster. Input is paired FASTQ or aligned BAM (e.g. GDC RNA-seq BAMs).

```
check_strandedness.py ──► stranded_mode per sample (qc/strandedness.tsv → manifest)
                            BAM: RSeQC on its alignments · FASTQ: kallisto + RSeQC

FASTQ pair ─┐
            ├─► step 1: BOWTIE ───────────────────► tmp_<N>/<sample>.bam ──► step 2: TELESCOPE ──► <sample>-TE_counts.tsv
BAM ────────┘   BAM → FASTQ (samtools, verified)     <sample>.stranded_mode    telescope assign --stranded_mode <mode>
                trim (Trim Galore)                                           (samples run concurrently)
                bowtie2 -k 100 | samtools view -b
```

Samples are processed in manifests of 4. Each manifest `N` gets one step-1
job and one step-2 job; `tmp_<N>/` carries the BAMs from step 1 to step 2 and
is deleted only after every sample in that manifest has a counts file.

## Repository layout

| Path | Contents |
|---|---|
| `code/run.py` | Pipeline entry point, run inside each Slurm job |
| `code/lib/worker_patch.py` | The pipeline itself (manifest-driven, FASTQ or BAM rows) |
| `code/lib/strandedness.py` | Strandedness measurement and the `--stranded_mode` call (shared by step 1 and `check_strandedness.py`) |
| `code/lib/bam_input.py` | BAM → paired FASTQ.gz with read-count verification |
| `code/lib/cmd.py` | Command-line arguments |
| `code/lib/worker.py` | Original upstream version (GDC download path); not used |
| `skills/telescope/` | Scripts, job templates and config for running a cohort — **use these** |
| `osc_job_generation.ipynb`, `jobs_step*/`, `manifests/` | Legacy job generator and example outputs; out of date |
| `requirements.txt` | Frozen `immuno_telescope` conda environment |

## Requirements

- OSC account `PDE0005`, partition `nextgen`
- Conda environment `immuno_telescope` (Python 3.6, telescope-ngs `4cf1859`,
  bowtie2 2.5.4, samtools 1.9, Trim Galore, RSeQC — see `requirements.txt`)
- Modules `miniconda3/24.1.2-py310` and `intel/2021.10.0`
- References in `/fs/ess/PDE0005/telescope_refs/data`:

| File | Role |
|---|---|
| `Indexes/gencode.v39_bowtie2/human` | bowtie2 index — **the genome reads are aligned to** |
| `REF/HG38_HERV_LINE_all_families_telescope_ann.gtf` | HERV/LINE annotation for Telescope (28,513 loci) |
| `REF/gencode.v39.annotation.gtf`, `REF/gencode.v39.transcripts.fa` | Strandedness check (gene BED for RSeQC; kallisto index for FASTQ) |
| `REF/GRCh38.p13.genome.chr.fa` | Passed as `--genome` but not read by the pipeline |

The bowtie2 index contains only the 25 primary chromosomes of GRCh38.p13
(chr1–22, X, Y, M), i.e. `GRCh38.p13.genome.chr.fa`, which is the first 25
records of GENCODE's `GRCh38.p13.genome.fa`. Scaffolds, patches and alternate
haplotypes are excluded, so 765 annotated loci that lie only on those
sequences always have zero counts.

## Running a cohort

All scripts are in `skills/telescope/scripts/`, take `--config` (default
`./telescope.config.json`); all but `check_fastq.py` accept `--dry-run`. Run
them from the project directory.

The scripts and `code/` run on Python 3.6 (the `immuno_telescope` env) and
on the system `python3` (3.9) alike; keep new code 3.6-compatible — e.g.
`subprocess.run(stdout=PIPE, stderr=PIPE, universal_newlines=True)`, not
`capture_output=`/`text=`, which need 3.7. `check_strandedness.py` also needs
the env's tools (samtools, RSeQC, kallisto) on `PATH`, so activate the env
for it.

**0. Configure.** Create a new project directory and a config in it:

```bash
cp skills/telescope/config.example.json /fs/scratch/PDE0005/projects/MYPROJECT/telescope.config.json
```

Set `project`, `code_dir`, `ref_dir`, `index_dir` and `job_prefix`, and the
input: `samples_dir` for FASTQ, or `"input_type": "bam"` with `bam_table` (a
TSV with `id` and `path` columns, e.g. a gdc-download `files.tsv`). Optional
`done_marker` (e.g. `/path/results/{sample}.done`) is touched per sample when
its counts exist. `code_dir` must point to an up-to-date copy of this repo's
`code/` — each job copies `run.py` and `lib/` from there when it starts.

Make that copy a **snapshot inside the project** rather than pointing at this
checkout: a job reads `code_dir` when it *starts*, so switching branches or
pulling while jobs are queued changes what they run.

```bash
cp -r code /fs/scratch/PDE0005/projects/MYPROJECT/code
git rev-parse --short HEAD > /fs/scratch/PDE0005/projects/MYPROJECT/code/SNAPSHOT_OF
```

A BAM cohort, e.g. one gdc-download batch (samples are named by GDC file id):

```json
{
  "project":    "/fs/scratch/PDE0005/projects/MYPROJECT",
  "code_dir":   "/fs/scratch/PDE0005/projects/MYPROJECT/code",
  "ref_dir":    "/fs/ess/PDE0005/telescope_refs/data/REF",
  "index_dir":  "/fs/ess/PDE0005/telescope_refs/data/Indexes",
  "input_type": "bam",
  "bam_table":  "/path/to/download/batches/batch_001/files.tsv",
  "job_prefix": "myproj",
  "done_marker": "/path/to/results/{sample}.done"
}
```

Leave `done_marker` empty for pilots whose inputs must not be cleaned up.

**1. Check FASTQ integrity** (`gzip -t` on every file; results are cached;
FASTQ cohorts only):

```bash
python check_fastq.py --slurm --submit     # full cohort, as a job array
python check_fastq.py -j 8                 # a few files, locally
```

**2. Decide each sample's `--stranded_mode`** (`qc/strandedness.tsv`):

```bash
conda activate immuno_telescope
python check_strandedness.py -j 8                 # a batch of BAMs: seconds each, run here
python check_strandedness.py --slurm --submit     # large FASTQ cohorts; then: --merge
```

Review the report (RF / FR / None / UNCLEAR, by sequencing center for BAMs).
Hand overrides go in `qc/strandedness_overrides.tsv`.

**3. Build manifests** (`manifests/manifest_<N>.tsv`, 4 samples each, with a
`stranded_mode` column):

```bash
python make_manifests.py --require-qc --require-strand     # --require-qc: FASTQ only
```

Manifest numbers are permanent. To add samples later use `--append`; never
renumber a project that has submitted jobs.

**4. Generate job scripts** (`jobs_step1/job_<N>.sh`, `jobs_step2/job_<N>.sh`):

```bash
python make_jobs.py --step both
```

**5. Submit**, naming the manifests explicitly:

```bash
python submit.py --step 1 --ids 0-9 --dry-run
python submit.py --step 1 --ids 0-9
# after step 1 finishes for those manifests:
python submit.py --step 2 --ids 0-9
```

`submit.py` checks inputs, queue state and resources before submitting and
aborts the batch if anything is wrong. It also:

- **submits as Unix group `submit_group`** (default: the account) via `sg`, so
  each job's `$PFSDIR` lands in that group's scratch. Without it, `$PFSDIR` is
  created under your *primary* group's scratch whatever `--account` is.
- **sizes each step-2 job from its own samples**: read pairs from the bowtie2
  logs → estimated Telescope peak memory (`s2_mem_base_gb` +
  `s2_mem_gb_per_mpairs` × M pairs, summed over samples running at once,
  × `s2_mem_margin`) → `--mem` and the cpus that memory costs on nextgen
  (4027 MB/cpu, ≤ 120 cpus). `--dry-run` prints the plan; it is recorded in
  `logs/step2_resources.tsv`. Submit step 2 once step 1 has finished.

Step 2 copies each sample to `output/TELESCOPE/<sample>/` as soon as it
finishes, and a re-submitted step 2 skips samples that already have counts,
so a killed job loses only the sample that was running. Each sample's peak
memory is logged as `telescope_max_rss_gb` for calibrating the sizing.

Default resources: step 1 — 24 CPUs, 64G, 96 h; step 2 — 12 CPUs, 47G, 144 h,
4 samples at a time.

## Outputs

Under the project directory:

| Path | Contents |
|---|---|
| `output/TELESCOPE/<sample>/<sample>-TE_counts.tsv` | Final counts per locus |
| `output/TELESCOPE/<sample>/<sample>-run_stats.tsv` | Telescope run statistics |
| `output/BOWTIE/<sample>/<sample>_bowties2.log` | Alignment summary |
| `output/{BOWTIE,TELESCOPE}/logs_<N>.tsv` | Per-manifest status and timings |
| `logs/step{1,2}/job_<N>-<jobid>.{out,err}` | Slurm logs |
| `qc/fastq_check.tsv` | FASTQ integrity results |

**Health check:** for each sample, `total_fragments` in the run stats should
equal the read-pair count on the first line of the bowtie2 log. The command is
in `skills/telescope/SKILL.md`.

## Things that must not change

- **Telescope always gets an explicit `--stranded_mode`.** Step 1 writes it
  next to each BAM (`tmp_<N>/<sample>.stranded_mode`); step 2 refuses a sample
  without one. Before this change every sample ran as `None` even though
  step 1 measured strandedness (1459 of 1503 `telescope_TCC00573` samples
  were RF), so older counts are unstranded counts — do not mix them with
  stranded results.
- **BAM → FASTQ keeps unmapped reads and drops secondary/supplementary
  records**, and is verified by read counts; samtools exits 0 even when it
  writes nothing.

- **Do not sort the step-1 BAM by coordinate.** Telescope groups a read's
  alignments by taking consecutive records with the same name. bowtie2's raw
  output is already grouped that way, so step 1 writes it unsorted and without
  an index. Coordinate-sorted BAMs make almost every alignment its own
  fragment; step 2 refuses them. Results produced before this fix (commit
  `c61d881`) are biased and must not be combined with newer results.
- **Keep the Intel runtime preload in the step-2 template.** Telescope's C
  extension needs `libimf`, `libsvml` and `libirc` via `LD_PRELOAD`
  (`LD_LIBRARY_PATH` alone does not work). The template checks this before
  starting.
- **Stay within 47G for 12 CPUs on `nextgen`.** The partition allows 4027 MB
  per CPU; asking for more memory is billed as extra CPUs.
- **Keep `telescope assign --ncpu 1`.** Multi-CPU mode is broken; parallelism
  comes from running several samples at once.

More detail, including troubleshooting failed jobs, is in
[`skills/telescope/SKILL.md`](skills/telescope/SKILL.md).

## Known issues

- The step-2 memory coefficients come from two pilot samples and are
  generous on purpose; refit them from `telescope_max_rss_gb` once more
  manifests have run.
- Pandas reads the stranded mode `None` in `logs_<N>.tsv` as missing (NaN). Use
  `keep_default_na=False` when reading those logs.

## Credits

Original pipeline code by Rosario Distefano (2022).
