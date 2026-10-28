#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Patched worker.py (manifest-first workflow)
-----------------------------------------
This version supports running the Bowtie2 → Telescope pipeline starting from
**paired-end FASTQs or aligned BAMs** listed in a manifest TSV with columns:

    sample\tfq1\tfq2\tbam\tstranded_mode

A row uses either fq1+fq2 or bam. `bam` rows are converted to FASTQ inside
step 1 (lib/bam_input.py). `stranded_mode` is Telescope's --stranded_mode
(None, RF, FR, R, F); when it is blank, step 1 measures it (lib/strandedness.py:
RSeQC on the BAM for BAM rows, kallisto on the FASTQ for FASTQ rows). Old
three-column manifests (sample, fq1, fq2) still work: their mode is measured.

Step 1 hands the mode to step 2 next to the BAM, as tmp_<N>/<sample>.stranded_mode,
and step 2 passes it to `telescope assign --stranded_mode`. Step 2 refuses a
sample without it rather than silently running Telescope's default (None).

It keeps the original class shape (Worker) and depends on a CmdParser-like
object that provides CLI-style attributes. Newly required attributes:

- cmd_obj.manifest: path to the manifest TSV (e.g., manifest_0.tsv)
- cmd_obj.samples_dir: directory containing FASTQ files referenced by manifest

Existing attributes used:
- gtf, genome, transcript, herv_gtf, bowtiew2_idx
- n_cores, trimgalore_n_cores, sample_check_n_cores, telescope_n_cores
- out_dir, scratch_dir
- workflows (e.g., ["BOWTIE", "TELESCOPE"]) 

Major changes from the original:
- Remove GDC download + BAM→FASTQ path; start from FASTQs.
- Add `prepare_samples_info_manifest()` that reads the simple manifest.
- Add `stage_import_fastq()` to lay out FASTQs under FASTQ/<sample>/
  as <sample>_R1.fq.gz and <sample>_R2.fq.gz (via symlink or copy).
- Strandedness decides Telescope's --stranded_mode (it used to be measured
  and then ignored, so every sample ran as None).
- BAM input (converted to FASTQ in step 1).
- Add `prepare_telescope_samples_info_from_scratch()` to discover BAMs
  produced by Bowtie2 in scratch for the Telescope stage.

Author of original: Rosario Distefano
Patched by: (assistant)
"""

from concurrent.futures import ProcessPoolExecutor
import errno
from functools import partial
import glob
import gzip
from hashlib import md5
import json
import math
import os
import pandas as pd
import pysam
from shutil import move, rmtree, copyfileobj, copyfile, copytree
import subprocess
import sys
import time
from typing import Dict, List, Any
import traceback

try:
    from . import strandedness as strand_lib  # type: ignore
    from .bam_input import bam_to_fastq  # type: ignore
except ImportError:  # run as a loose module
    import strandedness as strand_lib  # type: ignore
    from bam_input import bam_to_fastq  # type: ignore

try:
    # keep relative import if this file is part of a package
    from .cmd import CmdParser  # type: ignore
except Exception:  # pragma: no cover
    CmdParser = object  # fallback for static analysis / loose usage


class PatchedWorker:
    
    def __init__(self, cmd_obj):
        """
        Worker class (FASTQ manifest-driven)

        Arguments:
        -----------
        - cmd_obj: A CmdParser-like object exposing required attributes.
        """
        self.cmd_obj: CmdParser = cmd_obj  # type: ignore
        self.cohort_ids: List[str] = []
        self.sample_check_n_cores: int = getattr(self.cmd_obj, 'sample_check_n_cores', 2)
        self.trimgalore_n_cores: int = getattr(self.cmd_obj, 'trimgalore_n_cores', 4)
        self.gtf_fpath: str = self.cmd_obj.gtf
        self.genome_dna_fpath: str = getattr(self.cmd_obj, 'genome_dna', '')
        self.genome_fpath: str = self.cmd_obj.genome
        self.seed: int = getattr(self.cmd_obj, 'seed', 1)
        self.transcript_fpath: str = self.cmd_obj.transcript
        self.bowtiew2_genome_dpath: str = self.cmd_obj.bowtiew2_idx
        self.n_cores: int = self.cmd_obj.n_cores
        self.herv_gtf_fpath: str = self.cmd_obj.herv_gtf
        self.telescope_n_cores: int = self.cmd_obj.telescope_n_cores
        self.workflows: List[str] = list(self.cmd_obj.workflows)
        self.fastq_samples_limit: int = 6
        self.bam_samples_limit: int = max(1, getattr(self.cmd_obj, 'telescope_n_samples', 1))

        # NEW: manifest + samples_dir
        self.manifest_tsv: str = getattr(self.cmd_obj, 'manifest', '')
        self.samples_dir: str = getattr(self.cmd_obj, 'samples_dir', '')

        # Strandedness -> telescope assign --stranded_mode
        self.strand_n_reads: int = getattr(self.cmd_obj, 'strand_n_reads', 1000000)
        self.kallisto_index: str = getattr(self.cmd_obj, 'kallisto_index', '') or 'kallisto_gencode.idx'
        self.bam2fq_exclude_flags: str = getattr(self.cmd_obj, 'bam2fq_exclude_flags', '0xB00')
        self.gene_bed_fpath: str = 'genes.protein_coding.bed12'   # derived from --gtf when needed

    # ------------------------------------------------------------------
    # Manifest + sample info helpers
    # ------------------------------------------------------------------
    def prepare_samples_info_manifest(self) -> Dict[str, Any]:
        """Read the manifest and build the per-sample info map.

        Columns are matched by name (case-insensitive): sample, and either
        fq1+fq2 or bam; stranded_mode is optional. Relative paths resolve
        against `self.samples_dir`.
        """
        if not self.manifest_tsv:
            raise ValueError("'manifest' argument is required for the BOWTIE workflow")

        df = pd.read_csv(self.manifest_tsv, sep='\t', dtype=str, keep_default_na=False)
        df = df.rename(columns={c: c.strip().lower() for c in df.columns})
        if 'sample' not in df.columns:
            raise ValueError("Manifest %s has no 'sample' column" % self.manifest_tsv)
        for c in ('fq1', 'fq2', 'bam', 'stranded_mode'):
            if c not in df.columns:
                df[c] = ''

        def _abs(p):
            p = str(p).strip()
            if p and not os.path.isabs(p) and self.samples_dir:
                p = os.path.join(self.samples_dir, p)
            return p

        sample_info: Dict[str, Any] = {}
        for _, row in df.iterrows():
            s = str(row['sample']).strip()
            if not s:
                continue
            fq1, fq2, bam = _abs(row['fq1']), _abs(row['fq2']), _abs(row['bam'])
            mode = str(row['stranded_mode']).strip()
            if bam and (fq1 or fq2):
                raise ValueError("%s: give either bam or fq1/fq2, not both" % s)
            if not bam and not (fq1 and fq2):
                raise ValueError("%s: needs fq1 and fq2, or bam" % s)
            if mode and mode not in strand_lib.TELESCOPE_MODES:
                raise ValueError("%s: stranded_mode %r is not one of %s"
                                 % (s, mode, ', '.join(strand_lib.TELESCOPE_MODES)))
            sample_info[s] = {
                'sample_id': s,
                'input_type': 'bam' if bam else 'fastq',
                'fq1': fq1, 'fq2': fq2, 'bam': bam,
                'cohort_id': '',
                'stranded_mode': mode,                        # '' -> measured in step 1
                'stranded_mode_source': 'manifest' if mode else '',
            }
        return sample_info

    def stage_import_bam(self, fastq_dpath: str, n_threads: int, sample_info: Dict[str, Any], sample: str) -> Dict[str, Any]:
        """BAM -> FASTQ/<sample>/<sample>_R1.fq.gz / _R2.fq.gz (see lib/bam_input.py)."""
        s_time = time.time()
        bam = sample_info[sample]['bam']
        try:
            r = bam_to_fastq(bam, os.path.join(fastq_dpath, sample), sample,
                             threads=n_threads, exclude_flags=self.bam2fq_exclude_flags)
            d = {'status': 'Ok', 'note': '', 'sample_id': sample, 'R1_src': bam, 'R2_src': bam,
                 'bam2fq_pairs': r['pairs'], 'bam2fq_singletons_dropped': r['singletons_dropped']}
        except Exception as e:
            d = {'status': 'Err', 'note': 'BAM->FASTQ: %s' % e, 'sample_id': sample,
                 'R1_src': bam, 'R2_src': bam, 'bam2fq_pairs': 0, 'bam2fq_singletons_dropped': 0}
        d['fastq_import_time_in_sec'] = round(time.time() - s_time, 2)
        print('\t'.join(str(d[c]) for c in ('status', 'note', 'sample_id', 'bam2fq_pairs',
                                              'bam2fq_singletons_dropped', 'fastq_import_time_in_sec')))
        return d

    def stage_import_fastq(self, fastq_dpath: str, sample: str, fq1: str, fq2: str) -> Dict[str, Any]:
        """
        Prepare FASTQ layout expected by downstream functions:
        FASTQ/<sample>/<sample>_R1.fq.gz and <sample>_R2.fq.gz
        Uses symlinks when possible; falls back to copy.
        """
        s_time = time.time()
        cols = ['status', 'note', 'sample_id', 'R1_src', 'R2_src', 'fastq_import_time_in_sec']

        o_dpath = os.path.join(fastq_dpath, sample)
        if os.path.isdir(o_dpath):
            self.delete_file_or_dir(o_dpath)
        self.create_dir(o_dpath)

        r1_dst = os.path.join(o_dpath, f'{sample}_R1.fq.gz')
        r2_dst = os.path.join(o_dpath, f'{sample}_R2.fq.gz')

        def _place(src: str, dst: str):
            if not os.path.isfile(src):
                raise FileNotFoundError(f"Missing FASTQ: {src}")
            try:
                os.symlink(os.path.abspath(src), dst)
            except Exception:
                copyfile(src, dst)

        try:
            _place(fq1, r1_dst)
            _place(fq2, r2_dst)
            d = {
                'status': 'Ok', 'note': '', 'sample_id': sample,
                'R1_src': fq1, 'R2_src': fq2,
                'fastq_import_time_in_sec': round(time.time() - s_time, 2)
            }
        except Exception as e:
            d = {
                'status': 'Err', 'note': str(e), 'sample_id': sample,
                'R1_src': fq1, 'R2_src': fq2,
                'fastq_import_time_in_sec': round(time.time() - s_time, 2)
            }
        print('\t'.join([str(d[c]) for c in cols]))
        return d

    # ------------------------------------------------------------------
    # Existing helpers (kept mostly intact)
    # ------------------------------------------------------------------
    def get_fastqgz_nrows(self, fpath: str) -> int:
        if not os.path.isfile(fpath):
            return 0
        i = -1
        with gzip.open(fpath, 'rb') as f:
            for i, _ in enumerate(f):
                pass
        return i + 1

    def check_fastq_gz(self, fpath: str) -> bool:
        if not os.path.isfile(fpath):
            return False
        try:
            cmd = f"gzip -vt {fpath}"
            log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).decode('utf-8')
            msg = log.replace(f'{fpath}:', '').strip()
            return True if msg == 'OK' else False
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"command '{e.cmd}' return with error (code {e.returncode}): {e.output}")

    @staticmethod
    def _detect_input(o_dpath: str, sample: str, read_no: int) -> str:
        suffixes = [
            f"_{read_no}.fq.gz",
            f"_R{read_no}.fq.gz",
            f"_{read_no}.fastq.gz",
            f"_R{read_no}.fastq.gz",
        ]
        for suffix in suffixes:
            candidate = os.path.join(o_dpath, f"{sample}{suffix}")
            if os.path.isfile(candidate):
                return candidate
        raise FileNotFoundError(
            f"Missing FASTQ for sample '{sample}' read {read_no} in {o_dpath}"
        )

    def fastq_trimming(self, sample_info: Dict[str, Any], fastq_dpath: str, n_cores: int, sample: str) -> Dict[str, Any]:
        """ Trimgalore trimming """
        s_time = time.time()
        cols = ['status', 'note', 'sample_id', 'cohort_id', 'R1_trim_size', 'R2_trim_size', 'fastqtrim_time_in_sec']

        o_dpath = os.path.join(fastq_dpath, sample)

        def _trimmed_output_name(input_path: str, read_no: int) -> str:
            base = os.path.basename(input_path)
            if base.endswith('.fastq.gz'):
                stem = base[:-len('.fastq.gz')]
                ext = '.fastq.gz'
            elif base.endswith('.fq.gz'):
                stem = base[:-len('.fq.gz')]
                ext = '.fq.gz'
            else:
                raise ValueError(f"Unexpected FASTQ extension: {base}")
            return os.path.join(o_dpath, f"{stem}_val_{read_no}{ext}")

        r1_input = self._detect_input(o_dpath, sample, 1)
        r2_input = self._detect_input(o_dpath, sample, 2)
        cmd = (
            "trim_galore "
            "--paired "
            "--no_report_file "
            "--gzip "
            f"--cores {n_cores} "
            f"-o {o_dpath} "
            f"{r1_input} "
            f"{r2_input}"
        )
        try:
            log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
            with open(os.path.join(o_dpath, f'{sample}_trimgalore.log'), 'w') as f:
                f.write(log.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"command '{e.cmd}' return with error (code {e.returncode}): {e.output}")

        # Normalize output filenames back to <sample>_R1/_R2
        r1_trimmed = _trimmed_output_name(r1_input, 1)
        r2_trimmed = _trimmed_output_name(r2_input, 2)

        if not os.path.isfile(r1_trimmed) or not os.path.isfile(r2_trimmed):
            raise FileNotFoundError(
                f"Trimmed FASTQ not found for sample '{sample}' in {o_dpath}"
            )

        r1_dst = os.path.join(o_dpath, f'{sample}_R1.fq.gz')
        r2_dst = os.path.join(o_dpath, f'{sample}_R2.fq.gz')

        if os.path.isfile(r1_dst):
            os.remove(r1_dst)
        if os.path.isfile(r2_dst):
            os.remove(r2_dst)

        os.rename(r1_trimmed, r1_dst)
        os.rename(r2_trimmed, r2_dst)

        r1 = os.path.getsize(os.path.join(o_dpath, f'{sample}_R1.fq.gz')) if os.path.isfile(os.path.join(o_dpath, f'{sample}_R1.fq.gz')) else 0
        r2 = os.path.getsize(os.path.join(o_dpath, f'{sample}_R2.fq.gz')) if os.path.isfile(os.path.join(o_dpath, f'{sample}_R2.fq.gz')) else 0

        d = {
            'status': 'Ok' if r1 and r2 else 'Err',
            'note': '' if r1 and r2 else 'Trimgalore error',
            'sample_id': sample,
            'cohort_id': sample_info[sample].get('cohort_id', ''),
            'R1_trim_size': r1,
            'R2_trim_size': r2,
            'fastqtrim_time_in_sec': round(time.time() - s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d

    def resolve_strandedness(self, sample_info: Dict[str, Any], fastq_dpath: str, sample: str) -> Dict[str, Any]:
        """Decide this sample's Telescope --stranded_mode.

        The manifest value wins (check_strandedness.py measured it before the
        job, or someone set it by hand). Otherwise measure it here: RSeQC on
        the source BAM for BAM rows, kallisto on the staged (untrimmed) FASTQ
        for FASTQ rows. UNCLEAR calls run as None.
        """
        s_time = time.time()
        info = sample_info[sample]
        d = {'status': 'Ok', 'note': '', 'sample_id': sample,
             'stranded_mode': info.get('stranded_mode', ''),
             'stranded_mode_source': info.get('stranded_mode_source', ''),
             'strand_call': '', 'strand_frac_rf': '', 'strand_frac_fr': '', 'strand_note': ''}
        if not d['stranded_mode']:
            try:
                if info['input_type'] == 'bam':
                    r = strand_lib.infer_bam(info['bam'], self.gene_bed_fpath, self.strand_n_reads)
                else:
                    o_dpath = os.path.join(fastq_dpath, sample)
                    r = strand_lib.infer_fastq(
                        self._detect_input(o_dpath, sample, 1), self._detect_input(o_dpath, sample, 2),
                        self.gtf_fpath, self.transcript_fpath, self.kallisto_index, self.strand_n_reads)
                d.update({'stranded_mode': r['telescope_mode'], 'stranded_mode_source': r['method'],
                          'strand_call': r['stranded_mode'], 'strand_frac_rf': r['frac_rf'],
                          'strand_frac_fr': r['frac_fr'], 'strand_note': r['note']})
            except Exception as e:
                d.update({'status': 'Err', 'note': 'strandedness: %s' % e})
        info['stranded_mode'] = d['stranded_mode']
        d['strandedness_time_in_sec'] = round(time.time() - s_time, 2)
        print('\t'.join(str(d[c]) for c in ('status', 'note', 'sample_id', 'stranded_mode',
                                              'stranded_mode_source', 'strand_note', 'strandedness_time_in_sec')))
        return d

    def telescope_input(self, sample_info: Dict[str, Any], fastq_dpath: str, output_dpath: str, n_cores: int, genome_fpath: str, scratch_dpath: str, sample: str) -> Dict[str, Any]:
        """ Run Bowtie2 and produce sorted BAM + index; move to scratch. """
        s_time = time.time()
        cols = [
            'status', 'note', 'sample_id', 'cohort_id', 'bowties2',
            'bowties2_sam_size', 'bowties2_time_in_sec', 'samtools_sam_to_bam',
            'samtools_bam_size', 'samtools_sam_to_bam_time_in_sec',
            'samtools_bai', 'samtools_bai_time_in_sec'
        ]

        telescope_dpath = os.path.join(output_dpath, sample)
        if os.path.isdir(telescope_dpath):
            self.delete_file_or_dir(telescope_dpath)
        self.create_dir(telescope_dpath)

        # Align and compress to BAM in a single pass, in bowtie2's own output
        # order. Do NOT coordinate-sort this BAM: Telescope groups a fragment's
        # alignments by taking CONSECUTIVE records with the same query name
        # (telescope.utils.alignment.fetch_bundle). bowtie2 already emits every
        # alignment of a read pair as one contiguous block (also with -p > 1),
        # which is exactly what Telescope needs; `samtools sort -n` gives
        # identical counts at the cost of a full sort. A coordinate sort
        # scatters each read's alignments, so nearly every record became its
        # own fragment and multi-mappers were counted as unique.
        #
        # Two details matter here:
        #   * bowtie2 writes its summary to stderr and the SAM to stdout, so
        #     the two must be kept apart - merging them corrupts the stream.
        #   * pipefail is required because a shell pipeline reports only the
        #     LAST command's exit status. Without it, a bowtie2 failure part
        #     way through leaves samtools view writing a valid but truncated
        #     BAM and exiting 0, i.e. silently incomplete counts.
        view_n_cores = max(2, n_cores // 6)
        bt2_log_fpath = os.path.join(telescope_dpath, f'{sample}_bowties2.log')
        bam_fpath = os.path.join(telescope_dpath, f'{sample}.bam')
        cmd = (
            "set -o pipefail; "
            "bowtie2 "
            "-k 100 "
            f"-p {n_cores} "
            f"-x {genome_fpath} "
            f"-1 {os.path.join(fastq_dpath, sample, f'{sample}_R1.fq.gz')} "
            f"-2 {os.path.join(fastq_dpath, sample, f'{sample}_R2.fq.gz')} "
            f"--seed {self.seed} "
            f"2> {bt2_log_fpath} "
            "| "
            "samtools view "
            "-b "
            f"-@ {view_n_cores} "
            f"-o {bam_fpath} "
            "-"
        )
        print(cmd)
        log = subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')
        with open(os.path.join(telescope_dpath, f'{sample}_sam_to_bam.log'), 'w') as f:
            f.write(log.decode('utf-8'))

        # No SAM is materialised any more. The column is retained at 0 so the
        # logs_*.tsv schema stays identical to the previously profiled samples.
        # Align+compress is one fused step, so its wall time is reported under
        # bowties2_time_in_sec and samtools_sam_to_bam_time_in_sec is ~0.
        bowties2_sam_size = 0
        bowties2_sam_time = round(time.time() - s_time, 2)
        s_time = time.time()

        bam_fpath = os.path.join(telescope_dpath, f'{sample}.bam')
        if not os.path.isfile(bam_fpath):
            d = {
                'status': 'Err', 'note': 'No BAM file found',
                'sample_id': sample, 'cohort_id': sample_info[sample].get('cohort_id', ''),
                'bowties2': 'Ok', 'bowties2_sam_size': bowties2_sam_size,
                'bowties2_time_in_sec': bowties2_sam_time,
                'samtools_sam_to_bam': 'Err', 'samtools_bam_size': 0,
                'samtools_sam_to_bam_time_in_sec': round(time.time() - s_time, 2),
                'samtools_bai': '', 'samtools_bai_time_in_sec': .0
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d

        samtools_bam_size = os.path.getsize(bam_fpath)
        samtools_bam_time = round(time.time() - s_time, 2)

        # No index: a BAM in aligner order cannot be indexed, and Telescope
        # reads it front to back (until_eof) so it does not need one. The
        # samtools_bai columns are kept as 'Skipped' for schema stability.

        # Move into scratch for Telescope stage
        src_fpath = os.path.join(telescope_dpath, f'{sample}.bam')
        dst_fpath = os.path.join(scratch_dpath, f'{sample}.bam')
        move(src_fpath, dst_fpath)

        # Clean remaining duplicates if any
        if os.path.exists(src_fpath):
            self.delete_file_or_dir(src_fpath)

        d = {
            'status': 'Ok', 'note': '', 'sample_id': sample,
            'cohort_id': sample_info[sample].get('cohort_id', ''),
            'bowties2': 'Ok', 'bowties2_sam_size': bowties2_sam_size,
            'bowties2_time_in_sec': bowties2_sam_time,
            'samtools_sam_to_bam': 'Ok', 'samtools_bam_size': samtools_bam_size,
            'samtools_sam_to_bam_time_in_sec': samtools_bam_time,
            'samtools_bai': 'Skipped', 'samtools_bai_time_in_sec': .0
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d

    def telescope(self, sample_info: Dict[str, Any], output_dpath: str, n_cores: int, herv_gtf_fpath: str, sample: str) -> Dict[str, Any]:
        """ Run Telescope assign on a prepared BAM in the per-sample outdir. """
        s_time = time.time()
        cols = ['status', 'note', 'sample_id', 'cohort_id', 'stranded_mode', 'telescope',
                'telescope_time_in_sec', 'telescope_max_rss_gb']

        # Never fall back to Telescope's default here: a sample whose mode was
        # lost would be counted unstranded without anyone noticing.
        mode = sample_info.get(sample, {}).get('stranded_mode', '')
        if mode not in strand_lib.TELESCOPE_MODES:
            raise RuntimeError("%s: no valid stranded_mode (%r); re-run step 1 for this sample" % (sample, mode))

        telescope_dpath = os.path.join(output_dpath, sample)
        tmp_dpath = os.path.join(telescope_dpath, 'tmp')
        if os.path.isdir(tmp_dpath):
            self.delete_file_or_dir(tmp_dpath)
        self.create_dir(tmp_dpath)

        # Refuse a coordinate-sorted BAM (e.g. one left in tmp_<N> by an older
        # step 1): Telescope would run to completion on it and silently treat
        # each alignment record as its own fragment.
        bam_fpath = os.path.join(telescope_dpath, f'{sample}.bam')
        with pysam.AlignmentFile(bam_fpath, check_sq=False) as bam:
            sort_order = bam.header.to_dict().get('HD', {}).get('SO', '')
        if sort_order == 'coordinate':
            raise RuntimeError(
                f"{sample}: {bam_fpath} is coordinate-sorted; Telescope needs each "
                "read's alignments contiguous. Re-run step 1 for this sample.")

        # /usr/bin/time -v records this sample's own peak memory. sacct's MaxRSS
        # cannot: it is one number per job, and samples run concurrently.
        time_fpath = os.path.join(telescope_dpath, f'{sample}_telescope.time')
        timer = f"/usr/bin/time -v -o {time_fpath} " if os.path.exists('/usr/bin/time') else ''
        cmd = (
            f"{timer}telescope assign "
            f"{os.path.join(telescope_dpath, f'{sample}.bam')} "
            f"{herv_gtf_fpath} "
            f"--ncpu {n_cores} "
            f"--tempdir {tmp_dpath} "
            f"--outdir {telescope_dpath} "
            f"--exp_tag {sample} "
            f"--stranded_mode {mode}"
        )
        log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
        with open(os.path.join(telescope_dpath, f'{sample}_telescope.log'), 'w') as f:
            f.write(log.decode('utf-8'))

        # Cleanup tmp + intermediate BAMs written by telescope (if any)
        for fp in (
            os.path.join(telescope_dpath, f'{sample}.bam'),
            os.path.join(telescope_dpath, f'{sample}.bam.bai'),
            os.path.join(telescope_dpath, f'{sample}-checkpoint.npz'),
        ):
            if os.path.exists(fp):
                self.delete_file_or_dir(fp)
        if os.path.isdir(tmp_dpath):
            self.delete_file_or_dir(tmp_dpath)

        d = {
            'status': 'Ok', 'note': '', 'sample_id': sample,
            'cohort_id': sample_info.get(sample, {}).get('cohort_id', ''),
            'stranded_mode': mode,
            'telescope': 'Ok', 'telescope_time_in_sec': round(time.time() - s_time, 2),
            'telescope_max_rss_gb': self.max_rss_gb(time_fpath),
        }
        # The sample's own log row, next to its counts, so it survives even if
        # the job dies before logs.tsv is copied back.
        with open(os.path.join(telescope_dpath, f'{sample}_telescope_summary.tsv'), 'w') as f:
            f.write('\t'.join(cols) + '\n' + '\t'.join(str(d[c]) for c in cols) + '\n')
        self.publish_sample(telescope_dpath, sample)
        print('\t'.join([str(d[c]) for c in cols]))
        return d

    @staticmethod
    def max_rss_gb(time_fpath: str):
        """Peak RSS in GB from a /usr/bin/time -v report ('' if unavailable)."""
        try:
            with open(time_fpath) as f:
                for line in f:
                    if 'Maximum resident set size' in line:
                        return round(int(line.rsplit(':', 1)[1]) / 1024 ** 2, 2)
        except (OSError, ValueError):
            pass
        return ''

    def publish_sample(self, telescope_dpath: str, sample: str) -> None:
        """Copy one finished sample's outputs to --final_dir right away.

        Results used to be copied back only at the end of the job, so a job
        killed late (e.g. out of memory on its largest sample) lost every
        sample, finished or not. Copy to <sample>.partial and rename, so a
        half-copied directory is never mistaken for a finished one.
        """
        final_root = getattr(self.cmd_obj, 'final_dir', '')
        if not final_root:
            return
        dst = os.path.join(final_root, sample)
        tmp = dst + '.partial'
        for p in (tmp, dst):
            if os.path.exists(p):
                self.delete_file_or_dir(p)
        self.create_dir(final_root)
        copytree(telescope_dpath, tmp)
        os.rename(tmp, dst)

    # ------------------------------------------------------------------
    # Small utilities
    # ------------------------------------------------------------------
    def samples_chunking(self, samples, chunk_size):
        return [samples[i:i + chunk_size] for i in range(0, len(samples), chunk_size)]

    def prepare_telescope_samples_info_from_scratch(self, scratch_path: str) -> Dict[str, Any]:
        """Discover samples that have a BAM in scratch (output of Bowtie2 stage),
        with the --stranded_mode step 1 wrote next to it ('' if missing)."""
        info: Dict[str, Any] = {}
        for bam in glob.glob(os.path.join(scratch_path, '*.bam')):
            sample = os.path.basename(bam)[:-len('.bam')]
            if os.path.isfile(bam):
                info[sample] = {'sample_id': sample, 'cohort_id': '',
                                'stranded_mode': self.read_stranded_mode(scratch_path, sample)}
        return info

    @staticmethod
    def stranded_mode_fpath(scratch_path: str, sample: str) -> str:
        return os.path.join(scratch_path, f'{sample}.stranded_mode')

    def write_stranded_mode(self, scratch_path: str, sample: str, mode: str) -> None:
        fp = self.stranded_mode_fpath(scratch_path, sample)
        with open(fp + '.tmp', 'w') as f:
            f.write(mode + '\n')
        os.replace(fp + '.tmp', fp)

    def read_stranded_mode(self, scratch_path: str, sample: str) -> str:
        try:
            with open(self.stranded_mode_fpath(scratch_path, sample)) as f:
                return f.read().strip()
        except OSError:
            return ''

    # ------------------------------------------------------------------
    # Orchestration
    # ------------------------------------------------------------------
    def run(self) -> None:
        if 'BOWTIE' not in self.workflows and 'TELESCOPE' not in self.workflows:
            return

        # Initialization
        o_path: str = self.cmd_obj.out_dir
        scratch_path: str = self.cmd_obj.scratch_dir

        if os.path.isdir(o_path):
            self.delete_file_or_dir(o_path)
        self.create_dir(o_path)

        for workflow in self.workflows:
            w_path = os.path.join(o_path, workflow)
            if os.path.isdir(w_path):
                self.delete_file_or_dir(w_path)
            self.create_dir(w_path)

        logs_map: Dict[str, Any] = {}

        # =====================================================================
        # BOWTIE (BAM→FASTQ | FASTQ → strandedness → trim → Bowtie2 → BAM in scratch)
        # =====================================================================
        if 'BOWTIE' in self.workflows:
            sample_info = self.prepare_samples_info_manifest()

            fastq_dpath = 'FASTQ'
            if os.path.isdir(fastq_dpath):
                self.delete_file_or_dir(fastq_dpath)
            self.create_dir(fastq_dpath)

            # -----------------------------------------------------------------
            # Stage inputs as FASTQ/<sample>/<sample>_R{1,2}.fq.gz:
            # FASTQ rows are symlinked, BAM rows are converted (and verified).
            # -----------------------------------------------------------------
            s_time = time.time()
            for sample, info in sample_info.items():
                logs_map[sample] = {'status': 'Ok', 'note': '', 'sample_id': sample,
                                    'cohort_id': info.get('cohort_id', ''), 'input_type': info['input_type']}
            fq_samples = [s for s, i in sample_info.items() if i['input_type'] == 'fastq']
            bam_samples = [s for s, i in sample_info.items() if i['input_type'] == 'bam']
            for sample in fq_samples:
                info = sample_info[sample]
                result = self.stage_import_fastq(fastq_dpath, sample, info['fq1'], info['fq2'])
                logs_map[sample].update(result)
            if bam_samples:
                n_workers = max(1, min(len(bam_samples), self.n_cores // 6))
                func = partial(self.stage_import_bam, fastq_dpath, max(2, self.n_cores // n_workers), sample_info)
                with ProcessPoolExecutor(max_workers=n_workers) as executor:
                    for result in executor.map(func, bam_samples, chunksize=1):
                        logs_map[result['sample_id']].update(result)
            print("Input staging (wall-time): " f"{round(time.time() - s_time, 2)} (s)")
            sample_to_process: List[str] = [s for s in sample_info if logs_map[s]['status'] == 'Ok']

            if not sample_to_process:
                self.save_logs(logs_map, o_path)
                return

            # -----------------------------------------------------------------
            # Strandedness -> telescope --stranded_mode. Taken from the manifest
            # when set; otherwise measured here. A sample whose mode cannot be
            # decided is dropped: Telescope must not guess.
            # -----------------------------------------------------------------
            s_time = time.time()
            need = [s for s in sample_to_process if not sample_info[s]['stranded_mode']]
            try:
                if any(sample_info[s]['input_type'] == 'bam' for s in need):
                    strand_lib.bed_from_gtf(self.gtf_fpath, self.gene_bed_fpath)
                if any(sample_info[s]['input_type'] == 'fastq' for s in need):
                    # build once, before the parallel calls, not once per sample
                    strand_lib.ensure_kallisto_index(self.transcript_fpath, self.kallisto_index)
            except Exception as e:
                traceback.print_exc()
                for s in need:
                    logs_map[s].update({'status': 'Err', 'note': 'strandedness setup: %s' % e})
            func = partial(self.resolve_strandedness, sample_info, fastq_dpath)
            todo = [s for s in sample_to_process if logs_map[s]['status'] == 'Ok']
            with ProcessPoolExecutor(max_workers=max(1, min(4, len(todo)))) as executor:
                for result in executor.map(func, todo, chunksize=1):
                    logs_map[result['sample_id']].update(result)
                    sample_info[result['sample_id']]['stranded_mode'] = result['stranded_mode']
            print("Strandedness (wall-time): " f"{round(time.time() - s_time, 2)} (s)")
            sample_to_process = [s for s in sample_to_process if logs_map[s]['status'] == 'Ok']

            if not sample_to_process:
                self.save_logs(logs_map, o_path)
                return

            # -----------------------------------------------------------------
            # FASTQ trimming
            # -----------------------------------------------------------------
            sample_to_process2: List[str] = []
            try:
                n_workers = max(1, int(math.floor(self.n_cores / self.trimgalore_n_cores)))
                func = partial(self.fastq_trimming, sample_info, fastq_dpath, self.trimgalore_n_cores)
                s_time = time.time()
                with ProcessPoolExecutor(max_workers=n_workers) as executor:
                    for result in executor.map(func, sample_to_process, chunksize=1):
                        d = {**logs_map[result['sample_id']], **result}
                        logs_map[result['sample_id']] = d
                print("FastQ samples trimming (wall-time): " f"{round(time.time() - s_time, 2)} (s)")
                sample_to_process2 = [s for s in sample_to_process if logs_map[s]['status'] == 'Ok']
            except Exception:
                traceback.print_exc()

            if not sample_to_process2:
                self.save_logs(logs_map, o_path)
                return

            # -----------------------------------------------------------------
            # Bowtie2 → BAM in scratch
            # -----------------------------------------------------------------
            try:
                s_time = time.time()
                for sample in sample_to_process2:
                    result = self.telescope_input(
                        sample_info,
                        fastq_dpath,
                        os.path.join(o_path, 'BOWTIE'),
                        self.n_cores,
                        self.bowtiew2_genome_dpath,
                        scratch_path,
                        sample,
                    )
                    d = {**logs_map[result['sample_id']], **result}
                    logs_map[result['sample_id']] = d
                    if result['status'] == 'Ok':
                        # handoff to step 2, next to the BAM it belongs to
                        self.write_stranded_mode(scratch_path, sample, sample_info[sample]['stranded_mode'])
                print("Bowtie2 processing (wall-time): " f"{round(time.time() - s_time, 2)} (s)")
            except Exception:
                traceback.print_exc()
                self.save_logs(logs_map, o_path)

        # =====================================================================
        # TELESCOPE (BAM from scratch → Telescope assign)
        # =====================================================================
        if 'TELESCOPE' in self.workflows:
            sample_info_tele = self.prepare_telescope_samples_info_from_scratch(scratch_path)
            init_samples: List[str] = []

            for sample, info in sample_info_tele.items():
                bam_fpath = os.path.join(scratch_path, f'{sample}.bam')
                final_root = getattr(self.cmd_obj, 'final_dir', '')
                done = os.path.join(final_root, sample, f'{sample}-TE_counts.tsv') if final_root else ''
                if done and os.path.isfile(done) and os.path.getsize(done) > 0:
                    # finished by an earlier (e.g. killed) run of this job; its
                    # outputs were already published to --final_dir
                    logs_map[sample] = {'status': 'Ok', 'note': 'already done (kept from an earlier run)',
                                        'sample_id': sample, 'cohort_id': info.get('cohort_id', ''),
                                        'stranded_mode': info['stranded_mode'], 'telescope': 'Skipped'}
                    print('\t'.join(['Ok', 'already done', sample]))
                    continue
                if info['stranded_mode'] not in strand_lib.TELESCOPE_MODES:
                    logs_map[sample] = {
                        'status': 'Err', 'sample_id': sample, 'cohort_id': info.get('cohort_id', ''),
                        'stranded_mode': info['stranded_mode'],
                        'note': 'no valid %s next to the BAM; re-run step 1 for this sample'
                                % os.path.basename(self.stranded_mode_fpath(scratch_path, sample))}
                    print('\t'.join(['Err', logs_map[sample]['note'], sample]))
                    continue
                if os.path.isfile(bam_fpath):
                    tele_dpath = os.path.join(o_path, 'TELESCOPE', sample)
                    if os.path.isdir(tele_dpath):
                        self.delete_file_or_dir(tele_dpath)
                    self.create_dir(tele_dpath)
                    # Link, don't copy: the BAM is already on the same scratch
                    # filesystem, and copying a manifest's -k 100 BAMs (100-270 GB)
                    # into $PFSDIR only doubled step 2's disk use. telescope()
                    # removes the link afterwards; the BAM itself stays in tmp_<N>
                    # until the job template reclaims it.
                    os.symlink(os.path.abspath(bam_fpath), os.path.join(tele_dpath, f'{sample}.bam'))
                    init_samples.append(sample)
                    logs_map.setdefault(sample, {
                        'status': 'Ok', 'note': '', 'sample_id': sample, 'cohort_id': info.get('cohort_id', '')
                    })

            if not init_samples:
                self.save_logs(logs_map, o_path)
                return

            s_time = time.time()
            for samples_chunk in self.samples_chunking(init_samples, self.bam_samples_limit):
                try:
                    func = partial(self.telescope, sample_info_tele, os.path.join(o_path, 'TELESCOPE'), self.telescope_n_cores, self.herv_gtf_fpath)
                    with ProcessPoolExecutor(max_workers=self.bam_samples_limit) as executor:
                        for result in executor.map(func, samples_chunk, chunksize=1):
                            d = {**logs_map[result['sample_id']], **result}
                            logs_map[result['sample_id']] = d
                except Exception:
                    traceback.print_exc()
                    self.save_logs(logs_map, o_path)
            print("TELESCOPE (wall-time): " f"{round(time.time() - s_time, 2)} (s)")

        self.save_logs(logs_map, o_path)

    # ------------------------------------------------------------------
    # Filesystem helpers
    # ------------------------------------------------------------------
    def save_logs(self, runs_log, dpath):
        df = pd.DataFrame([v for _, v in runs_log.items()])
        df.to_csv(os.path.join(dpath, 'logs.tsv'), sep='\t', encoding='utf-8', index=False)

    def delete_file_or_dir(self, path: str) -> None:
        if os.path.exists(path):
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                try:
                    rmtree(path)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
            else:
                raise ValueError(f'{path} is not a file or directory.')

    def create_dir(self, path: str):
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def print_error(self, *args, **kwargs):
        print(*args, file=sys.stderr, **kwargs)

    def print_error_and_abort(self, *args, **kwargs):
        print(*args, file=sys.stderr, **kwargs)
        sys.exit(1)
