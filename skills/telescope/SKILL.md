---
name: telescope
description: Run the two-step Bowtie2 -> Telescope HERV/LINE quantification pipeline on an OSC Slurm cluster. Use when asked to validate FASTQ files, build manifests, generate step-1 (BOWTIE) or step-2 (TELESCOPE) job scripts, or submit/resubmit those jobs. Covers gzip integrity checks, 4-sample manifests, SLURM script generation, and guarded batch submission.
---

# Telescope pipeline

Two steps, one manifest at a time:

| step | workflow | reads | writes | typical cost |
|------|----------|-------|--------|--------------|
| 1 | `BOWTIE` | `<sample>_1.fastq.gz` / `_2.fastq.gz` | `tmp_<N>/<sample>.bam` + `.bai` | 24 cpu, 64 G, ~12–36 h |
| 2 | `TELESCOPE` | `tmp_<N>/*.bam` | `output/TELESCOPE/<sample>/<sample>-TE_counts.tsv` | 12 cpu, 47 G, 1–3 days |

Step 2 deletes `tmp_<N>` only after every sample in that manifest has a counts
file, so `tmp_<N>` is the handoff between the steps and the rerun safety net.

## Before anything else: read the config

Every script takes `--config <path>` (default `./telescope.config.json`).
If the project has no config yet, copy `config.example.json` and edit it.
**Never hardcode paths into a command — put them in the config**, so the four
scripts agree with each other and with whatever ran last month.

```bash
cp skills/telescope/config.example.json telescope.config.json
```

## 1. Validate the FASTQ files

```bash
python skills/telescope/scripts/check_fastq.py --config telescope.config.json
```

Runs `gzip -t` on every `*.fastq.gz` in `samples_dir` and writes
`qc/fastq_check.tsv` (`path  size  mtime  status  detail`).

`gzip -t` decompresses the whole file. At ~10 GB per file and a few thousand
files this is tens of TB of I/O — **do not run it unqualified on a login node.**
Pick one:

- `--slurm` — generate and submit a job array, one task per chunk of files.
  This is the right default for a full cohort. Add `--array-chunk N` to set
  files per task (default 20) and `--submit` to actually submit.
- `-j N` — local parallel check. Fine for a pilot or a handful of re-downloads.
- `--sample-only SL123,SL456` — check named samples only.

Results are cached by `(path, size, mtime)`, so re-running skips files that
already passed. Use `--recheck` to force.

Read the summary it prints. `status=BAD` means the transfer is corrupt and the
sample must be re-downloaded — it will not fail until hours into bowtie2
otherwise.

## 2. Generate manifests

```bash
python skills/telescope/scripts/make_manifests.py --config telescope.config.json
```

Pairs `_1`/`_2` (or `_R1`/`_R2`) by prefix, sorts, and writes
`manifests/manifest_<N>.tsv` with `samples_per_manifest` (default 4) samples
each. Samples with no mate go to `manifests/UNPAIRED.txt` and are excluded.

**Manifest indices are permanent.** Job scripts, `tmp_<N>` directories, log
filenames and `logs_<N>.tsv` are all keyed to them, so renumbering silently
divorces finished work from its manifest. The script therefore refuses to touch
existing manifests unless told:

- default — writes only if `manifests/` has none; otherwise reports and exits
- `--append` — keeps every existing manifest byte-for-byte and adds only
  newly-discovered samples, numbered from the current maximum + 1
- `--force` — rewrites everything. Only for a project with no submitted jobs.

Add `--require-qc` to include only samples that passed step 1's `gzip -t`.

## 3. Generate job scripts

```bash
python skills/telescope/scripts/make_jobs.py --config telescope.config.json --step both
```

Renders `templates/step1.sh.tmpl` and `templates/step2.sh.tmpl` into
`jobs_step1/job_<N>.sh` and `jobs_step2/job_<N>.sh`, one per manifest,
`chmod 750`.

Like manifests, existing scripts are not overwritten by default:

- default — writes only the scripts that do not exist yet
- `--force` — rewrite all. Safe only for jobs that have not been submitted;
  a queued job reads its script at start time, so rewriting under a pending
  job changes what runs.
- `--ids 233-332` — restrict to a range or comma list
- `--dry-run` — print what would be written

Run this again after editing a template, with `--ids` limited to unsubmitted
manifests.

## 4. Submit

```bash
python skills/telescope/scripts/submit.py --config telescope.config.json \
    --step 2 --ids 233-332
```

**Only submit when the user has explicitly asked for it, and only the range
they named.** These jobs run for days and consume real allocation; never
submit "the rest" on your own initiative.

Preflight runs first and aborts the whole batch if anything fails:

| check | step 1 | step 2 |
|-------|--------|--------|
| job script exists and `bash -n` passes | ✓ | ✓ |
| manifest exists | ✓ | ✓ |
| every input present and non-empty | FASTQ pairs | `tmp_<N>/<s>.bam` + `.bai` |
| output not already complete | — | `-TE_counts.tsv` absent |
| not already in the queue | ✓ | ✓ |
| resource directives uniform across the batch | ✓ | ✓ |

Useful flags: `--dry-run` (preflight only), `--allow-partial` (submit the
manifests that pass instead of aborting), `--max N` (cap the batch).

On success it prints `submitted=N failed=0 range=<first>..<last>`. Record that
range — it is how the jobs are found later with `sacct`.

## Cluster facts that are easy to get wrong

- **`nextgen` enforces `MaxMemPerCPU=4027 MB`.** 12 cpus × 4027 MB = 47.2 G, so
  `--mem=47G` is the ceiling for a 12-cpu step-2 job. Ask for more and Slurm
  allocates extra cpus to cover the RAM and bills for them — one job asking
  256 G was billed 66 cores. The `partition=batch` on older scripts is INACTIVE
  and jobs sit forever.
- **`--gres=pfsdir:scratch`** gives a job-private `$PFSDIR`. Put the SAM and
  trimmed FASTQ there; only the BAM belongs in `tmp_<N>`. `pfsdir:ess` is no
  longer advertised.
- **`telescope assign --ncpu > 1` is broken.** Parallelism is per-sample:
  `--telescope_n_samples` runs the manifest's samples as concurrent processes,
  each single-threaded.
- **Telescope's C extension needs the Intel runtime.** The step-2 template
  `LD_PRELOAD`s `libimf/libsvml/libirc` and verifies with a `ctypes.CDLL` probe
  before doing any work. Do not strip that probe; without it the job dies
  hours in.
- Step 2 captures `run.py`'s exit status into `$rc` instead of letting
  `set -e` fire, so results are always copied back before the job reports
  failure. A one-sample failure previously discarded a whole manifest's worth
  of finished work. Keep that ordering in any template edit.

## Known upstream defect: fragment bundling

`alignment.fetch_bundle` groups **consecutive** records by `query_name`, but
step 1 pipes bowtie2 through `samtools sort` (coordinate order). Records for
one read end up scattered, so nearly every alignment record is treated as its
own fragment.

Symptom: `total_fragments / nmap_idx` in `<sample>-run_stats.tsv` sits near
1.0 instead of the expected 0.01–0.05.

```bash
find output/TELESCOPE -maxdepth 2 -name '*-run_stats.tsv' -exec head -1 {} \; \
| grep -oP '(?<=\t)nmap_idx:\K[0-9]+|(?<=\t)total_fragments:\K[0-9]+' | paste - - \
| awk '$1>0{r=$2/$1; n++; s+=r} END{printf "n=%d mean=%.3f\n", n, s/n}'
```

Note the `(?<=\t)` — without it the pattern also matches inside `nunmap_idx`
and silently misaligns every pair.

Fixing this means `samtools collate` or `sort -n` in step 1 and redoing step 2
for the whole cohort. **Raise it with the user; do not switch sort order on
your own** — mixing sort orders within a cohort makes samples incomparable,
which is worse than a consistent known bias.

## Diagnosing a failed job

1. `sacct -j <jobid> -o JobID,State,Elapsed,MaxRSS,ExitCode -S <date> -E now`
2. Read `logs/step<N>/job_<M>-<jobid>.err` — the real traceback lives here.
3. For step 2, look for `output/TELESCOPE/<sample>/<sample>_telescope.FAILED.log`,
   which holds Telescope's own stderr for that sample.
4. `tmp_<M>` is kept whenever a job fails, so a rerun is just resubmitting
   `jobs_step2/job_<M>.sh`.

Two cautions from this pipeline's history:

- `sacct`'s **`MaxRSS` is the peak of a single process, not the sum** across
  the concurrent samples. It does not prove or disprove memory exhaustion on
  its own.
- Failures that correlate with sample size are not automatically caused by
  sample size. A batch of large-sample jobs once failed together and every one
  succeeded unchanged on resubmission; the size ordering was a coincidence of
  which jobs were running during a filesystem event. **Resubmit once before
  theorising.**

## Files

```
skills/telescope/
  SKILL.md                    this file
  config.example.json         copy to <project>/telescope.config.json and edit
  scripts/
    _common.py                config loading, id ranges, manifest parsing
    check_fastq.py            gzip -t sweep, cached, local or Slurm array
    make_manifests.py         paired-FASTQ discovery -> manifest_<N>.tsv
    make_jobs.py              templates -> jobs_step{1,2}/job_<N>.sh
    submit.py                 preflight + sbatch a named range
  templates/
    step1.sh.tmpl             BOWTIE
    step2.sh.tmpl             TELESCOPE
```

Every script accepts `--dry-run`. Use it.
