---
name: telescope
description: Run the two-step Bowtie2 -> Telescope HERV/LINE quantification pipeline on an OSC Slurm cluster, from paired FASTQ or aligned BAM (e.g. GDC) input. Use when asked to validate FASTQ files, determine library strandedness / Telescope --stranded_mode, build manifests, generate step-1 (BOWTIE) or step-2 (TELESCOPE) job scripts, or submit/resubmit those jobs. Covers gzip integrity checks, strandedness calls, BAM->FASTQ conversion, 4-sample manifests, SLURM script generation, and guarded batch submission.
---

# Telescope pipeline

Two steps, one manifest at a time:

| step | workflow | reads | writes | typical cost |
|------|----------|-------|--------|--------------|
| 1 | `BOWTIE` | a FASTQ pair **or** a BAM per sample (BAM → FASTQ inside the job) | `tmp_<N>/<sample>.bam` (aligner order, no index) + `tmp_<N>/<sample>.stranded_mode` | 24 cpu, 64 G, ~12–36 h |
| 2 | `TELESCOPE` | `tmp_<N>/*.bam` + `.stranded_mode` | `output/TELESCOPE/<sample>/<sample>-TE_counts.tsv` | 12 cpu, 47 G, 1–3 days |

Step 2 deletes `tmp_<N>` only after every sample in that manifest has a counts
file, so `tmp_<N>` is the handoff between the steps and the rerun safety net.

Order of work: **check_fastq (FASTQ only) → check_strandedness → make_manifests
--require-strand → make_jobs → submit step 1 → submit step 2.**

## Input: FASTQ or BAM

Set `input_type` in the config:

- `"fastq"` (default) — pairs `<sample>_1.fastq.gz` / `_2.fastq.gz` in `samples_dir`.
- `"bam"` — one BAM per row of `bam_table`, a TSV with a sample column
  (`bam_sample_col`, default `id`) and a path column (`bam_path_col`, default
  `path`). A gdc-download `batches/batch_NNN/files.tsv` works as is: samples
  are named by GDC file id, so `done_marker: ".../results/{sample}.done"` is
  exactly the marker gdc-download's cleanup looks for.

Step 1 converts BAM rows with `samtools collate | samtools fastq` into the
job's `$PFSDIR` (`code/lib/bam_input.py`). It keeps unmapped reads, drops
secondary/supplementary and vendor-QC-fail records (`bam2fq_exclude_flags`,
default `0xB00`), writes reverse-strand reads back in sequencing orientation,
and **fails the sample** unless R1 = R2 > 0 and 2×pairs + singletons equals the
BAM's primary reads — samtools itself exits 0 even when it writes nothing.
Paired-end only; single-end BAMs are refused.

## Before anything else: read the config

Every script takes `--config <path>` (default `./telescope.config.json`).
If the project has no config yet, copy `config.example.json` and edit it.
**Never hardcode paths into a command — put them in the config**, so the four
scripts agree with each other and with whatever ran last month.

```bash
cp skills/telescope/config.example.json telescope.config.json
```

**Snapshot the code into the project** and point `code_dir` at the copy, not
at this checkout. Jobs read `code_dir` when they start, so a branch switch or
pull while jobs are queued silently changes what they run:

```bash
cp -r code <project>/code && git rev-parse --short HEAD > <project>/code/SNAPSHOT_OF
```

**Python.** Every script runs under the env's Python 3.6 and the system 3.9.
`check_strandedness.py` needs the env active for its tools (samtools, RSeQC,
kallisto). When editing, stay 3.6-compatible: no `capture_output=`/`text=`
in `subprocess.run` (use `stdout=PIPE, stderr=PIPE, universal_newlines=True`),
no walrus operator, no `dict | dict`.

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

## 2. Decide each sample's Telescope `--stranded_mode`

```bash
conda activate immuno_telescope
python skills/telescope/scripts/check_strandedness.py -j 8               # BAMs: ~20 s each, fine here for a batch
python skills/telescope/scripts/check_strandedness.py --slurm --submit   # large FASTQ cohorts, as an array
python skills/telescope/scripts/check_strandedness.py --merge            # after the array finishes
```

Writes `qc/strandedness.tsv` (cached by input path, size, mtime) with, per
sample: the fraction of reads explained by RF (`1+-,1-+,2++,2--`) and FR
(`1++,1--,2+-,2-+`), the call, and — for BAMs — the sequencing centers in the
`@RG` header.

| input | method | time |
|---|---|---|
| BAM | RSeQC `infer_experiment.py` on its own alignments | ~20 s |
| FASTQ | `check_strandedness` (kallisto pseudo-alignment + `infer_experiment.py`); the kallisto index is built **once** at `kallisto_index` | ~1–2 min (+10 min once for the index) |

Calls (`code/lib/strandedness.py`, the same code step 1 uses): RF share of the
informative reads ≥ 0.8 → `RF`; ≤ 0.2 → `FR`; 0.35–0.65 → `None`; otherwise
`UNCLEAR`, which Telescope runs as `None` (counts both strands; never wrong about
orientation). UNCLEAR usually means libraries from different protocols merged
in one BAM (TARGET `PARLVL-09A`: BCGSC stranded + HAIB unstranded).

**Show the user the report before building manifests.** Strandedness follows the
sequencing protocol, not the project: in TARGET batch 1, 959 BAMs were RF and 41
were unstranded — all 31 OS samples (NCI-MELTZER), 5 HAIB, and 4 BCCAGSC/BCGSC
AML. Never set one value for a cohort.

To force a value, add `sample<TAB>stranded_mode<TAB>reason` to
`qc/strandedness_overrides.tsv`; overrides win over measured calls.

## 3. Generate manifests

```bash
python skills/telescope/scripts/make_manifests.py --config telescope.config.json
```

Pairs `_1`/`_2` (or `_R1`/`_R2`) by prefix — or reads `bam_table` — sorts,
and writes `manifests/manifest_<N>.tsv` with `samples_per_manifest` (default 4)
samples each, columns `sample fq1 fq2 bam stranded_mode` (a row uses fq1+fq2
or bam). Samples with no mate go to `manifests/UNPAIRED.txt` and are excluded.

`stranded_mode` comes from `qc/strandedness.tsv` / the overrides. With
`--require-strand` (recommended) samples without a current call are dropped and
listed; without it their mode is left blank and **measured inside step 1**.
Old three-column manifests still run; their modes are measured in step 1.

**Manifest indices are permanent.** Job scripts, `tmp_<N>` directories, log
filenames and `logs_<N>.tsv` are all keyed to them, so renumbering silently
divorces finished work from its manifest. The script therefore refuses to touch
existing manifests unless told:

- default — writes only if `manifests/` has none; otherwise reports and exits
- `--append` — keeps every existing manifest byte-for-byte and adds only
  newly-discovered samples, numbered from the current maximum + 1
- `--force` — rewrites everything. Only for a project with no submitted jobs.

Add `--require-qc` to include only samples that passed step 1's `gzip -t`.

## 4. Generate job scripts

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

## 5. Submit

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
| every input present and non-empty | FASTQ pair or BAM | `tmp_<N>/<s>.bam` |
| valid `stranded_mode` | manifest value (blank = measured in job) | `tmp_<N>/<s>.stranded_mode` |
| output not already complete | — | `-TE_counts.tsv` absent |
| not already in the queue | ✓ | ✓ |
| resource directives uniform across the batch | ✓ | ✓ |

Useful flags: `--dry-run` (preflight only), `--allow-partial` (submit the
manifests that pass instead of aborting), `--max N` (cap the batch),
`--no-auto-size` (step 2: keep the script's `--mem`/`--cpus`).

**Jobs are submitted as Unix group `submit_group`** (default: `account`), via
`sg <group> -c "sbatch ..."`. A job's `$PFSDIR` — where step 1 stages 75–230 GB
of FASTQ and step 2 copies its BAMs — is created under
`/fs/scratch/<job's group>/`, which is the user's *primary* group (here
`PCON0100`) unless the group is set, whatever `--account` says.

**Step 2 is sized per manifest.** Telescope's peak memory grows with sample
size (pilot: ~68 GB alone for a 285 M-pair sample; four ~80 M-pair samples fit
in 47 G). `submit.py --step 2` reads each remaining sample's read pairs from
`output/BOWTIE/<s>/<s>_bowties2.log`, estimates its peak as
`s2_mem_base_gb + s2_mem_gb_per_mpairs × M pairs`, and requests the sum over the
samples that run at once × `s2_mem_margin` — as cpus × `max_mem_per_cpu_mb`,
because nextgen sells memory only with cpus. If even `max_cpus_per_job` is too
little it lowers concurrency (`S2_CONCURRENT_OVERRIDE`). The plan is printed in
`--dry-run` and recorded in `logs/step2_resources.tsv`. Step 2 therefore cannot
be pre-queued with `afterok`: submit it once step 1 has finished.

The default coefficients are deliberately generous (two pilot data points).
Each sample's real peak is now logged as `telescope_max_rss_gb` (from
`/usr/bin/time -v`) in `output/TELESCOPE/logs_<N>.tsv` and
`<sample>_telescope_summary.tsv`; refit the two coefficients from those after
the first few dozen manifests.

On success it prints `submitted=N failed=0 range=<first>..<last>`. Record that
range — it is how the jobs are found later with `sacct`.

**Submit step 2 after step 1 has finished**, with `submit.py --step 2 --ids …`;
its sizing needs step 1's read counts. Keep a small table of manifest →
step-1 job → step-2 job in the project; job ids are how cost and failures are
traced later.

**A killed step-2 job keeps its finished samples.** Each sample is copied to
`output/TELESCOPE/<sample>/` (via `<sample>.partial`, then renamed) the moment
it finishes, and a re-submitted step 2 skips samples whose counts are already
there. Resubmitting after a failure therefore only redoes what was lost.

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

## Stranded mode: step 2 never guesses

Before this change the pipeline ran `check_strandedness` in step 1 and threw
the result away (its mapping only knew `unstranded`, and step 2 never read it),
so **every sample ran Telescope's default `--stranded_mode None`**. In
`telescope_TCC00573_redownloaded`, 1459 of 1503 samples were measured
`RF/fr-firststrand` and still counted unstranded. Counts from before this change
are therefore unstranded counts; do not mix them with stranded ones.

Now step 1 writes the decided mode to `tmp_<N>/<sample>.stranded_mode`, step 2
passes it to `telescope assign --stranded_mode`, and a sample without a valid
mode fails in step 2 (and in `submit.py --step 2` preflight) instead of running
as `None`. The mode used is recorded in `output/TELESCOPE/logs_<N>.tsv`.

## BAM order: never coordinate-sort the step-1 BAM

Telescope's `alignment.fetch_bundle` groups **consecutive** records by
`query_name`, so all alignments of a read pair must be contiguous. Bowtie2's
own output already is (also with `-p > 1`), so step 1 writes it straight to BAM
with `samtools view -b` — no sort, no `.bai` (such a BAM cannot be indexed, and
Telescope reads it front to back). `samtools sort -n` gives identical counts
but costs a full sort of a `-k 100` BAM; a **coordinate** sort scatters each
read's alignments so nearly every record becomes its own fragment and
multi-mappers are counted as unique. `telescope()` refuses a BAM whose header
says `SO:coordinate`.

Cohorts quantified before this change (e.g. `telescope_TCC00573_redownloaded`,
1503 samples) used coordinate-sorted BAMs. Their counts are biased and must not
be mixed with counts from the current pipeline — redo steps 1 and 2 for the
whole cohort.

Health check — with correct grouping, `total_fragments` in
`<sample>-run_stats.tsv` equals the read-pair count on the first line of
`output/BOWTIE/<sample>/<sample>_bowties2.log`. A coordinate-sorted run shows
several times more fragments than pairs.

```bash
for f in output/TELESCOPE/*/*-run_stats.tsv; do
  s=$(basename "$(dirname "$f")")
  frag=$(head -1 "$f" | grep -oP '(?<=\t)total_fragments:\K[0-9]+')
  pairs=$(head -1 "output/BOWTIE/$s/${s}_bowties2.log" 2>/dev/null | awk '{print $1}')
  echo -e "$s\t$frag\t$pairs"
done | awk -F'\t' '$3>0{r=$2/$3; n++; if(r>1.01||r<0.99)bad++} END{printf "n=%d off=%d\n", n, bad+0}'
```

## Known issues

- Pandas reads the stranded mode `None` in `logs_<N>.tsv` as NaN; pass
  `keep_default_na=False` when reading the logs.

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
    check_strandedness.py     stranded_mode per sample -> qc/strandedness.tsv
    make_manifests.py         FASTQ pairs or bam_table -> manifest_<N>.tsv
    make_jobs.py              templates -> jobs_step{1,2}/job_<N>.sh
    submit.py                 preflight + sbatch a named range
  templates/
    step1.sh.tmpl             BOWTIE
    step2.sh.tmpl             TELESCOPE
```

Every script accepts `--dry-run`. Use it.
