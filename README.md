# telescope-pipeline

Quantifies HERV and LINE retrotransposon expression from paired-end RNA-seq
with [Telescope](https://github.com/mlbendall/telescope), on the OSC Slurm
cluster.

```
FASTQ pairs ──► step 1: BOWTIE ─────────────────► tmp_<N>/<sample>.bam ──► step 2: TELESCOPE ──► <sample>-TE_counts.tsv
                trim (Trim Galore)                                          telescope assign
                strandedness (RSeQC)                                        (samples run concurrently)
                bowtie2 -k 100 | samtools view -b
```

Samples are processed in manifests of 4. Each manifest `N` gets one step-1
job and one step-2 job; `tmp_<N>/` carries the BAMs from step 1 to step 2 and
is deleted only after every sample in that manifest has a counts file.

## Repository layout

| Path | Contents |
|---|---|
| `code/run.py` | Pipeline entry point, run inside each Slurm job |
| `code/lib/worker_patch.py` | The pipeline itself (manifest-driven, starts from FASTQ) |
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
| `REF/gencode.v39.annotation.gtf`, `REF/gencode.v39.transcripts.fa` | Strandedness check |
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

**0. Configure.** Create a new project directory and a config in it:

```bash
cp skills/telescope/config.example.json /fs/scratch/PDE0005/projects/MYPROJECT/telescope.config.json
```

Set `project`, `samples_dir`, `code_dir`, `ref_dir`, `index_dir` and
`job_prefix`. `code_dir` must point to an up-to-date copy of this repo's
`code/` — each job copies `run.py` and `lib/` from there when it starts.

**1. Check FASTQ integrity** (`gzip -t` on every file; results are cached):

```bash
python check_fastq.py --slurm --submit     # full cohort, as a job array
python check_fastq.py -j 8                 # a few files, locally
```

**2. Build manifests** (`manifests/manifest_<N>.tsv`, 4 samples each):

```bash
python make_manifests.py --require-qc
```

Manifest numbers are permanent. To add samples later use `--append`; never
renumber a project that has submitted jobs.

**3. Generate job scripts** (`jobs_step1/job_<N>.sh`, `jobs_step2/job_<N>.sh`):

```bash
python make_jobs.py --step both
```

**4. Submit**, naming the manifests explicitly:

```bash
python submit.py --step 1 --ids 0-9 --dry-run
python submit.py --step 1 --ids 0-9
# after step 1 finishes for those manifests:
python submit.py --step 2 --ids 0-9
```

`submit.py` checks inputs, queue state and resources before submitting and
aborts the batch if anything is wrong.

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

## Credits

Original pipeline code by Rosario Distefano (2022).
