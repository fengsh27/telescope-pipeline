#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Patched worker.py (manifest-first workflow)
-----------------------------------------
This version supports running the Bowtie2 → Telescope pipeline starting
**from paired-end FASTQs** listed in a simple manifest TSV with columns:

    sample\tfq1\tfq2

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
- Enable strandedness mapping for RF/FR so samples are not dropped.
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
from shutil import move, rmtree, copyfileobj, copyfile
import subprocess
import sys
import time
from typing import Dict, List, Any
import traceback

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
        self.check_strand_n_workers: int = 2
        self.fastq_samples_limit: int = 6
        self.bam_samples_limit: int = 1

        # NEW: manifest + samples_dir
        self.manifest_tsv: str = getattr(self.cmd_obj, 'manifest', '')
        self.samples_dir: str = getattr(self.cmd_obj, 'samples_dir', '')

        # Keep bowtie strandedness mapping here if ever needed (not used below)
        # self.bowtie2_strandedness = {...}

        # Allow all three strandedness outcomes for Telescope
        self.telescope_strandedness: Dict[str, str] = {
            'unstranded': 'None', # We are only interested in unstranded samples
            # 'RF/fr-firststrand': 'RF',
            # 'FR/fr-secondstrand': 'FR',
        }

    # ------------------------------------------------------------------
    # Manifest + sample info helpers
    # ------------------------------------------------------------------
    def prepare_samples_info_manifest(self) -> Dict[str, Any]:
        """Read the simple manifest (sample\tfq1\tfq2) and build info map.
        FASTQ paths are resolved relative to `self.samples_dir` if not absolute.
        """
        if not self.manifest_tsv:
            raise ValueError("'manifest' argument is required for FASTQ mode")

        # Try reading as TSV; fall back to any whitespace in case
        try:
            df = pd.read_csv(self.manifest_tsv, sep='\t', dtype=str)
        except Exception:
            df = pd.read_csv(self.manifest_tsv, sep='\s+', engine='python', dtype=str)

        required_cols = {'sample', 'fq1', 'fq2'}
        if not required_cols.issubset(df.columns.str.lower()):
            # try case-insensitive mapping
            cols = {c.lower(): c for c in df.columns}
            missing = required_cols - set(cols.keys())
            if missing:
                raise ValueError(f"Manifest missing columns: {missing}")
            df = df.rename(columns={cols['sample']: 'sample', cols['fq1']: 'fq1', cols['fq2']: 'fq2'})
        else:
            # normalize exactly
            df = df.rename(columns={c: c.lower() for c in df.columns})

        sample_info: Dict[str, Any] = {}
        for _, row in df.iterrows():
            s = str(row['sample']).strip()
            fq1 = str(row['fq1']).strip()
            fq2 = str(row['fq2']).strip()

            if not os.path.isabs(fq1):
                fq1 = os.path.join(self.samples_dir, fq1) if self.samples_dir else fq1
            if not os.path.isabs(fq2):
                fq2 = os.path.join(self.samples_dir, fq2) if self.samples_dir else fq2

            sample_info[s] = {
                'sample_id': s,
                'fq1': fq1,
                'fq2': fq2,
                'cohort_id': '',
                'telescope_strandedness': ''
            }
        return sample_info

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

    def fastq_trimming(self, sample_info: Dict[str, Any], fastq_dpath: str, n_cores: int, sample: str) -> Dict[str, Any]:
        """ Trimgalore trimming """
        s_time = time.time()
        cols = ['status', 'note', 'sample_id', 'cohort_id', 'R1_trim_size', 'R2_trim_size', 'fastqtrim_time_in_sec']

        o_dpath = os.path.join(fastq_dpath, sample)
        cmd = (
            "trim_galore "
            "--paired "
            "--no_report_file "
            "--gzip "
            f"--cores {n_cores} "
            f"-o {o_dpath} "
            f"{os.path.join(o_dpath, f'{sample}_R1.fq.gz')} "
            f"{os.path.join(o_dpath, f'{sample}_R2.fq.gz')}"
        )
        try:
            log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
            with open(os.path.join(o_dpath, f'{sample}_trimgalore.log'), 'w') as f:
                f.write(log.decode('utf-8'))
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"command '{e.cmd}' return with error (code {e.returncode}): {e.output}")

        # Normalize output filenames back to <sample>_R1/_R2
        for file in os.listdir(o_dpath):
            if file.endswith(f'{sample}_R1_val_1.fq.gz'):
                os.rename(os.path.join(o_dpath, file), os.path.join(o_dpath, f'{sample}_R1.fq.gz'))
            elif file.endswith(f'{sample}_R2_val_2.fq.gz'):
                os.rename(os.path.join(o_dpath, file), os.path.join(o_dpath, f'{sample}_R2.fq.gz'))

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

    def check_strandedness(self, sample_info: Dict[str, Any], fastq_dpath: str, gtf_fpath: str, transcript_fpath: str, sample: str) -> Dict[str, Any]:
        """ Strandedness """
        s_time = time.time()
        cols = ['status', 'note', 'sample_id', 'cohort_id', 'strandedness', 'strandedness_time_in_sec']

        o_dpath = os.path.join(fastq_dpath, sample)
        strandedness = ''

        cmd = (
            "check_strandedness "
            f"--gtf {gtf_fpath} "
            f"--transcripts {transcript_fpath} "
            "-n 1000000 "
            f"--reads_1 {os.path.join(o_dpath, f'{sample}_R1.fq.gz')} "
            f"--reads_2 {os.path.join(o_dpath, f'{sample}_R2.fq.gz')} "
            f"--kallisto_index {sample}.kallisto_idx"
        )
        try:
            log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
            with open(os.path.join(o_dpath, f'{sample}_strandedness.log'), 'w') as f:
                f.write(log.decode('utf-8'))
                rows = log.decode('utf-8').strip().split('\n')
                if rows:
                    strandedness = rows[-1].strip().replace('Data is likely ', '')
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"command '{e.cmd}' return with error (code {e.returncode}): {e.output}")

        # Cleanup kallisto artifacts
        kallisto_idx = f'{sample}.kallisto_idx'
        if os.path.isfile(kallisto_idx):
            self.delete_file_or_dir(kallisto_idx)
        stranded_test = f'stranded_test_{sample}_R1'
        if os.path.isdir(stranded_test):
            self.delete_file_or_dir(stranded_test)

        d = {
            'status': 'Ok', 'note': '', 'sample_id': sample,
            'cohort_id': sample_info[sample].get('cohort_id', ''),
            'strandedness': strandedness,
            'strandedness_time_in_sec': round(time.time() - s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
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

        cmd = (
            "bowtie2 "
            "-k 100 "
            f"-p {n_cores} "
            f"-x {genome_fpath} "
            f"-1 {os.path.join(fastq_dpath, sample, f'{sample}_R1.fq.gz')} "
            f"-2 {os.path.join(fastq_dpath, sample, f'{sample}_R2.fq.gz')} "
            f"-S {os.path.join(telescope_dpath, f'{sample}.sam')} "
            f"--seed {self.seed}"
        )
        print(cmd)
        log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
        with open(os.path.join(telescope_dpath, f'{sample}_bowties2.log'), 'w') as f:
            f.write(log.decode('utf-8'))

        sam_fpath = os.path.join(telescope_dpath, f'{sample}.sam')
        if not os.path.isfile(sam_fpath):
            d = {
                'status': 'Err', 'note': 'No SAM file found',
                'sample_id': sample, 'cohort_id': sample_info[sample].get('cohort_id', ''),
                'bowties2': 'Err', 'bowties2_sam_size': 0,
                'bowties2_time_in_sec': round(time.time() - s_time, 2),
                'samtools_sam_to_bam': '', 'samtools_bam_size': 0,
                'samtools_sam_to_bam_time_in_sec': .0,
                'samtools_bai': '', 'samtools_bai_time_in_sec': .0
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d

        bowties2_sam_size = os.path.getsize(sam_fpath)
        bowties2_sam_time = round(time.time() - s_time, 2)

        # Sort BAM
        s_time = time.time()
        cmd = (
            "samtools view "
            f"-@ {n_cores} "
            f"-bS {os.path.join(telescope_dpath, f'{sample}.sam')} "
            "| "
            "samtools sort -o "
            f"{os.path.join(telescope_dpath, f'{sample}.bam')}"
        )
        log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
        with open(os.path.join(telescope_dpath, f'{sample}_sam_to_sorted_bam.log'), 'w') as f:
            f.write(log.decode('utf-8'))

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

        # Remove SAM to save space
        if os.path.isfile(sam_fpath):
            self.delete_file_or_dir(sam_fpath)

        # Index BAM
        s_time = time.time()
        cmd = (
            "samtools index "
            f"-@ {n_cores} "
            f"{os.path.join(telescope_dpath, f'{sample}.bam')} "
            f"{os.path.join(telescope_dpath, f'{sample}.bam.bai')}"
        )
        log = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
        with open(os.path.join(telescope_dpath, f'{sample}_samtools_bai.log'), 'w') as f:
            f.write(log.decode('utf-8'))

        bai_fpath = os.path.join(telescope_dpath, f'{sample}.bam.bai')
        if not os.path.isfile(bai_fpath):
            d = {
                'status': 'Err', 'note': 'No BAM.BAI file found',
                'sample_id': sample, 'cohort_id': sample_info[sample].get('cohort_id', ''),
                'bowties2': 'Ok', 'bowties2_sam_size': bowties2_sam_size,
                'bowties2_time_in_sec': bowties2_sam_time,
                'samtools_sam_to_bam': 'Ok', 'samtools_bam_size': samtools_bam_size,
                'samtools_sam_to_bam_time_in_sec': samtools_bam_time,
                'samtools_bai': 'Err', 'samtools_bai_time_in_sec': round(time.time() - s_time, 2)
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d

        # Move into scratch for Telescope stage
        src_fpath = os.path.join(telescope_dpath, f'{sample}.bam')
        dst_fpath = os.path.join(scratch_dpath, f'{sample}.bam')
        move(src_fpath, dst_fpath)
        src_fpath = os.path.join(telescope_dpath, f'{sample}.bam.bai')
        dst_fpath = os.path.join(scratch_dpath, f'{sample}.bam.bai')
        move(src_fpath, dst_fpath)

        # Clean remaining duplicates if any
        for fp in (os.path.join(telescope_dpath, f'{sample}.bam'), os.path.join(telescope_dpath, f'{sample}.bam.bai')):
            if os.path.exists(fp):
                self.delete_file_or_dir(fp)

        d = {
            'status': 'Ok', 'note': '', 'sample_id': sample,
            'cohort_id': sample_info[sample].get('cohort_id', ''),
            'bowties2': 'Ok', 'bowties2_sam_size': bowties2_sam_size,
            'bowties2_time_in_sec': bowties2_sam_time,
            'samtools_sam_to_bam': 'Ok', 'samtools_bam_size': samtools_bam_size,
            'samtools_sam_to_bam_time_in_sec': samtools_bam_time,
            'samtools_bai': 'Ok', 'samtools_bai_time_in_sec': round(time.time() - s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d

    def telescope(self, sample_info: Dict[str, Any], output_dpath: str, n_cores: int, herv_gtf_fpath: str, sample: str) -> Dict[str, Any]:
        """ Run Telescope assign on a prepared BAM in the per-sample outdir. """
        s_time = time.time()
        cols = ['status', 'note', 'sample_id', 'cohort_id', 'telescope', 'telescope_time_in_sec']

        telescope_dpath = os.path.join(output_dpath, sample)
        tmp_dpath = os.path.join(telescope_dpath, 'tmp')
        if os.path.isdir(tmp_dpath):
            self.delete_file_or_dir(tmp_dpath)
        self.create_dir(tmp_dpath)

        cmd = (
            "telescope assign "
            f"{os.path.join(telescope_dpath, f'{sample}.bam')} "
            f"{herv_gtf_fpath} "
            f"--ncpu {n_cores} "
            f"--tempdir {tmp_dpath} "
            f"--outdir {telescope_dpath} "
            f"--exp_tag {sample}"
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
            'telescope': 'Ok', 'telescope_time_in_sec': round(time.time() - s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d

    # ------------------------------------------------------------------
    # Small utilities
    # ------------------------------------------------------------------
    def samples_chunking(self, samples, chunk_size):
        return [samples[i:i + chunk_size] for i in range(0, len(samples), chunk_size)]

    def prepare_telescope_samples_info_from_scratch(self, scratch_path: str) -> Dict[str, Any]:
        """Discover samples that have BAM+BAI in scratch (output of Bowtie2 stage)."""
        info: Dict[str, Any] = {}
        for bam in glob.glob(os.path.join(scratch_path, '*.bam')):
            sample = os.path.basename(bam).replace('.bam', '')
            bai = os.path.join(scratch_path, f'{sample}.bam.bai')
            if os.path.isfile(bam) and os.path.isfile(bai):
                info[sample] = {'sample_id': sample, 'cohort_id': ''}
        return info

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
        # BOWTIE (FASTQ → trim → strandedness → Bowtie2 → BAM in scratch)
        # =====================================================================
        if 'BOWTIE' in self.workflows:
            sample_info = self.prepare_samples_info_manifest()

            fastq_dpath = 'FASTQ'
            if os.path.isdir(fastq_dpath):
                self.delete_file_or_dir(fastq_dpath)
            self.create_dir(fastq_dpath)

            # -----------------------------------------------------------------
            # Import/link FASTQs into expected layout
            # -----------------------------------------------------------------
            sample_to_process: List[str] = []
            for sample, info in sample_info.items():
                result = self.stage_import_fastq(fastq_dpath, sample, info['fq1'], info['fq2'])
                logs_map[sample] = {
                    'status': result['status'],
                    'note': result['note'],
                    'sample_id': sample,
                    'cohort_id': info.get('cohort_id', ''),
                }
                if result['status'] == 'Ok':
                    sample_to_process.append(sample)

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
                sample_to_process2 = [v['sample_id'] for _, v in logs_map.items() if v['status'] == 'Ok']
            except Exception:
                traceback.print_exc()

            if not sample_to_process2:
                self.save_logs(logs_map, o_path)
                return

            # -----------------------------------------------------------------
            # Strandedness (optional info; we do not drop samples if not matched)
            # -----------------------------------------------------------------
            try:
                func = partial(self.check_strandedness, sample_info, fastq_dpath, self.gtf_fpath, self.transcript_fpath)
                s_time = time.time()
                with ProcessPoolExecutor(max_workers=self.check_strand_n_workers) as executor:
                    for result in executor.map(func, sample_to_process2, chunksize=1):
                        d = {**logs_map[result['sample_id']], **result}
                        logs_map[result['sample_id']] = d
                        # record strandedness mapping if recognized
                        str_mode = result.get('strandedness', '')
                        if str_mode in self.telescope_strandedness:
                            sample_info[result['sample_id']]['telescope_strandedness'] = self.telescope_strandedness[str_mode]
                print("Samples strandedness check (wall-time): " f"{round(time.time() - s_time, 2)} (s)")
            except Exception:
                traceback.print_exc()

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
                bai_fpath = os.path.join(scratch_path, f'{sample}.bam.bai')
                if os.path.isfile(bam_fpath) and os.path.isfile(bai_fpath):
                    tele_dpath = os.path.join(o_path, 'TELESCOPE', sample)
                    if os.path.isdir(tele_dpath):
                        self.delete_file_or_dir(tele_dpath)
                    self.create_dir(tele_dpath)
                    copyfile(bam_fpath, os.path.join(tele_dpath, f'{sample}.bam'))
                    copyfile(bai_fpath, os.path.join(tele_dpath, f'{sample}.bam.bai'))
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
