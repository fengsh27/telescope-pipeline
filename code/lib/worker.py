#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  24 10:29:10 2022

@author: Rosario Distefano
"""

from concurrent.futures import ProcessPoolExecutor
import csv
import errno
from functools import partial, reduce
import glob
import gzip
from hashlib import md5
from itertools import chain
import json
import math
import os
import pandas as pd
import re
from shutil import move, rmtree, copyfileobj, copyfile
import subprocess
import sys
import time
from typing import Dict, List, Tuple, Union, Any, Optional
import traceback

from .cmd import CmdParser

class Worker():
        
    def __init__(self, cmd_obj):
        """
        Worker class

        Arguments:
        -----------
        - cmd_obj: A CmdParser class object.
        """
        self.cmd_obj: CmdParser = cmd_obj
        self.cohort_ids: List[str] = []
        self.sample_check_n_cores: int = self.cmd_obj.sample_check_n_cores
        self.trimgalore_n_cores: int = self.cmd_obj.trimgalore_n_cores
        self.gtf_fpath: str = self.cmd_obj.gtf
        self.genome_dna_fpath: str = self.cmd_obj.genome_dna
        self.genome_fpath: str = self.cmd_obj.genome
        self.seed: int = self.cmd_obj.seed
        self.transcript_fpath: str = self.cmd_obj.transcript
        self.bowtiew2_genome_dpath: str = self.cmd_obj.bowtiew2_idx
        self.n_cores: int = self.cmd_obj.n_cores
        self.herv_gtf_fpath: str = self.cmd_obj.herv_gtf
        self.telescope_n_cores: int = self.cmd_obj.telescope_n_cores
        self.workflows: List[str] = cmd_obj.workflows
        self.check_strand_n_workers: int = 2
        self.fastq_samples_limit: int = 6
        self.bam_samples_limit: int = 1
        """
        self.bowtie2_strandedness: Dict[str,str] = {
            'unstranded': '--fr',
            'RF/fr-firststrand': '--fr --nofw',
            'FR/fr-secondstrand': '--fr --norc'
        }

        self.telescope_strandedness: Dict[str,str] = {
            'unstranded': 'None',
            'RF/fr-firststrand': 'RF',
            'FR/fr-secondstrand': 'FR'
        }
        """
        self.telescope_strandedness: Dict[str,str] = {
            'unstranded': 'None'
        }
        

    def prepare_samples_info(self) -> Dict[str,Any]:
        with open(self.cmd_obj.gdc_mt, 'r') as f:
            mdata: List[Dict[str,Any]] = json.load(f)
                
        manifest_df: pd.DataFrame = pd.read_csv(
            self.cmd_obj.gdc_mf,
            delimiter='\t',
            encoding='utf-8'
        )

        manifest_df.set_index('id', inplace=True)
        manifest: Dict[str,Any] = manifest_df.to_dict('index')

        sample_info: Dict[str,Any] = {}

        for sample in mdata:
            if sample['file_id'] in manifest:
                sample_ids: List[str] = sorted(
                    [s['entity_submitter_id']
                    for s in sample['associated_entities']]
                )

                sample_info[sample_ids[0]] = {
                    'sample_id': sample_ids[0],
                    'file_id': sample['file_id'],
                    'file_name': sample['file_name'],
                    'cohort_id': manifest[sample['file_id']]['Project ID'],
                    'md5': manifest[sample['file_id']]['md5'],
                    'telescope_strandedness': ''
                }
                self.cohort_ids.append(
                    manifest[sample['file_id']]['Project ID'])

        self.cohort_ids = sorted(list(set(self.cohort_ids)))
        del manifest_df, manifest, mdata
        return sample_info


    def prepare_telescope_samples_info(self) -> Dict[str,Any]:
        manifest_df: pd.DataFrame = pd.read_csv(
            self.cmd_obj.gdc_mf,
            delimiter='\t',
            encoding='utf-8'
        )

        # Keep only unstranded samples
        manifest_df = manifest_df[
            (manifest_df.strandedness=='unstranded')
        ].reset_index(drop=True)

        manifest_df.set_index('sample_id', inplace=True)
        return manifest_df.to_dict('index')


    def sample_check(self, sample_info: Dict[str,Any], bam_dpath: str,
        n_cores: int, sample: str) -> Dict[str,Any]:
        """ Initial samples check """

        s_time: float = time.time()

        cols: List[str] = [
            'status', 'note', 'sample_id', 'cohort_id', 'bam_init_size',
            'is_paired', 'bamcheck_time_in_sec'
        ]

        bam_fpath: str = os.path.join(bam_dpath, f'{sample}.bam')
        
        if not os.path.isfile(bam_fpath):
            self.delete_file_or_dir(bam_fpath)
            d: Dict[str,Any] = {
                'status': 'Err',
                'note': 'No BAM found',
                'sample_id': sample,
                'cohort_id': '',
                'bam_init_size': '',
                'is_paired': '',
                'bamcheck_time_in_sec': [time.time()-s_time]
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d
        
        if sample not in sample_info:
            self.delete_file_or_dir(bam_fpath)
            d: Dict[str,Any] = {
                'status': 'Err',
                'note': 'No sample ID in metadata',
                'sample_id': sample,
                'cohort_id': '',
                'bam_init_size': '',
                'is_paired': '',
                'bamcheck_time_in_sec': [time.time()-s_time]
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d

        finfo: os.stat_result = os.stat(bam_fpath)
        init_bam_size: int = finfo.st_size
        del finfo
        
        checksum: str = md5(open(bam_fpath,'rb').read()).hexdigest()
        if checksum != sample_info[sample]['md5']:
            self.delete_file_or_dir(bam_fpath)
            d: Dict[str,Any] = {
                'status': 'Err',
                'note': 'Invalid checksum',
                'sample_id': sample,
                'cohort_id': sample_info[sample]['cohort_id'],
                'bam_init_size': init_bam_size,
                'is_paired': '',
                'bamcheck_time_in_sec': time.time()-s_time
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d
        
        # Check if paired end (n_rows must be > than 0)
        cmd: str = f"samtools " \
            " view " \
            f"-@ {n_cores} " \
            "-c " \
            "-f " \
            f"1 {bam_fpath}"

        n_rows = int(subprocess.check_output(cmd, shell=True
            ).decode('utf-8').strip())

        if n_rows == 0:
            self.delete_file_or_dir(bam_fpath)
            d: Dict[str,Any] = {
                'status': 'Err',
                'note': 'Not Paired-End',
                'sample_id': sample,
                'cohort_id': sample_info[sample]['cohort_id'],
                'bam_init_size': init_bam_size,
                'is_paired': 'n',
                'bamcheck_time_in_sec': time.time()-s_time
            }                
            print('\t'.join([str(d[c]) for c in cols]))
            return d
        
        d: Dict[str,Any] = {
            'status': 'Ok',
            'note': '',
            'sample_id': sample,
            'cohort_id': sample_info[sample]['cohort_id'],
            'bam_init_size': init_bam_size,
            'is_paired': 'y',
            'bamcheck_time_in_sec': round(time.time()-s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d


    def get_fastqgz_nrows(self, fpath: str) -> int:
        if not os.path.isfile(fpath):
            return 0

        with gzip.open(fpath, 'rb') as f:
            for i,_ in enumerate(f):
                pass
        return i+1


    def check_fastq_gz(self, fpath: str) -> bool:
        if not os.path.isfile(fpath):
            return False

        try:
            cmd: str = f"gzip -vt {fpath}"

            log: str = subprocess.check_output(
                cmd, stderr=subprocess.STDOUT, shell=True
                ).decode('utf-8')

            msg: str = log.replace(f'{fpath}:','').strip()
            return True if msg == 'OK' else False
            del log
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"command '{e.cmd}' return with error " \
                f"(code {e.returncode}): {e.output}")

        
    def bam_to_fastq(
        self,
        sample_info: Dict[str,Any],
        bam_dpath: str,
        fastq_dpath: str,
        sample: str
        ) -> Dict[str,Any]:
        
        s_time: float = time.time()

        cols: List[str] = [
            'status', 'note', 'sample_id', 'cohort_id', 'R1_init_size',
            'R2_init_size', 'bamtofastq_time_in_sec'
        ]

        o_dpath: str = os.path.join(fastq_dpath, sample)
        if os.path.isdir(o_dpath):
            self.delete_file_or_dir(o_dpath)        
        self.create_dir(o_dpath)
        
        cmd: str = "bamtofastq exclude=QCFAIL,SECONDARY,SUPPLEMENTARY " \
            "collate=1 " \
            f"outputperreadgroup=1 outputperreadgroupsuffixF={sample}_R1.fq " \
            f"outputperreadgroupsuffixF2={sample}_R2.fq " \
            f"outputdir={o_dpath} " \
            f"filename={os.path.join(bam_dpath, f'{sample}.bam')} " \
            "tryoq=1"

        try:
            log: str = subprocess.check_output(
                cmd, stderr=subprocess.STDOUT, shell=True)
            del log
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"command '{e.cmd}' return with error " \
                f"(code {e.returncode}): {e.output}")
            
        # Check how many FastQ bamtofastq generated
        n_files: int = len(glob.glob1(o_dpath, '*_R1.fq'))

        if n_files == 1:
            for file in os.listdir(o_dpath):
                if file.endswith(f'{sample}_R1.fq'):
                    os.rename(
                        os.path.join(o_dpath, file),
                        os.path.join(o_dpath, f'{sample}_R1.fq')
                    )  
                elif file.endswith(f'{sample}_R2.fq'):
                    os.rename(
                        os.path.join(o_dpath, file),
                        os.path.join(o_dpath, f'{sample}_R2.fq')
                    )
        elif n_files >1:
            fq_prefixes: List[str] = sorted([
                f.replace(f'{sample}_R1.fq', '')
                for f in os.listdir(o_dpath)
                if f.endswith(f'{sample}_R1.fq')
            ])

            print(f"Multiple FastQ files:\t{sample}\t{n_files}")

            r1_fpath: str = os.path.join(
                o_dpath, f'{fq_prefixes[0]}{sample}_R1.fq')
            r2_fpath: str = os.path.join(
                o_dpath, f'{fq_prefixes[0]}{sample}_R2.fq')

            os.rename(
                r1_fpath,
                os.path.join(o_dpath, f'{sample}_R1.fq')
            )  
            os.rename(
                r2_fpath,
                os.path.join(o_dpath, f'{sample}_R2.fq')
            )

            with open(os.path.join(o_dpath, f'{sample}_R1.fq'),'ab+') as wfd:
                for prefix in fq_prefixes[1:]:
                    f_path: str = os.path.join(
                        o_dpath, f'{prefix}{sample}_R1.fq')
                    with open(f_path,'rb') as fd:
                        copyfileobj(fd, wfd)
                        self.delete_file_or_dir(f_path)

            with open(os.path.join(o_dpath, f'{sample}_R2.fq'),'ab+') as wfd:
                for prefix in fq_prefixes[1:]:
                    f_path: str = os.path.join(
                        o_dpath, f'{prefix}{sample}_R2.fq')
                    with open(f_path,'rb') as fd:
                        copyfileobj(fd, wfd)
                        self.delete_file_or_dir(f_path)

        r1: int = 0
        r1_fpath: str = os.path.join(o_dpath, f'{sample}_R1.fq')
        if os.path.isfile(r1_fpath):
            finfo: os.stat_result = os.stat(r1_fpath)
            r1: int = finfo.st_size
            del finfo

        try:
            cmd: str = "gzip " \
                "-c " \
                f"{r1_fpath} " \
                "> " \
                f"{r1_fpath}.gz"

            log: str = subprocess.check_output(
                cmd, stderr=subprocess.STDOUT, shell=True)

            if not self.check_fastq_gz(
                os.path.join(o_dpath, f'{sample}_R1.fq.gz')):
                d: Dict[str,Any] = {
                    'status': 'Err',
                    'note': 'Invalid R1.fq.gz file',
                    'sample_id': sample,
                    'cohort_id': sample_info[sample]['cohort_id'],
                    'multiple_files': n_files,
                    'R1_init_size': r1,
                    'R2_init_size': 0,
                    'bamtofastq_time_in_sec': round(time.time()-s_time, 2)
                }
                print('\t'.join([str(d[c]) for c in cols]))
                return d

            self.delete_file_or_dir(os.path.join(o_dpath, f'{sample}_R1.fq'))
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"command '{e.cmd}' return with error " \
                f"(code {e.returncode}): {e.output}")

        r2: int = 0
        r2_fpath: str = os.path.join(o_dpath, f'{sample}_R2.fq')
        if os.path.isfile(r2_fpath):
            finfo: os.stat_result = os.stat(r2_fpath)
            r2: int = finfo.st_size
            del finfo

        try:
            cmd: str = "gzip " \
                "-c " \
                f"{r2_fpath} " \
                "> " \
                f"{r2_fpath}.gz"

            log: str = subprocess.check_output(
                cmd, stderr=subprocess.STDOUT, shell=True)
            
            if not self.check_fastq_gz(
                os.path.join(o_dpath, f'{sample}_R2.fq.gz')):
                d: Dict[str,Any] = {
                    'status': 'Err',
                    'note': 'Invalid R2.fq.gz file',
                    'sample_id': sample,
                    'cohort_id': sample_info[sample]['cohort_id'],
                    'multiple_files': n_files,
                    'R1_init_size': r1,
                    'R2_init_size': 0,
                    'bamtofastq_time_in_sec': round(time.time()-s_time, 2)
                }
                print('\t'.join([str(d[c]) for c in cols]))
                return d

            self.delete_file_or_dir(os.path.join(o_dpath, f'{sample}_R2.fq'))
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"command '{e.cmd}' return with error " \
                f"(code {e.returncode}): {e.output}")

        if r1 and r2:
            self.delete_file_or_dir(os.path.join(bam_dpath, f'{sample}.bam'))

        d: Dict[str,Any] = {
            'status': 'Ok' if r1 and r2 else 'Err',
            'note': '' if r1 and r2 else 'Invalid R1/R2 FASTQ files',
            'sample_id': sample,
            'cohort_id': sample_info[sample]['cohort_id'],
            'multiple_files': n_files,
            'R1_init_size': r1,
            'R2_init_size': r2,
            'bamtofastq_time_in_sec': round(time.time()-s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d


    def fastq_trimming(
        self,
        sample_info: Dict[str,Any],
        fastq_dpath: str,
        n_cores: int,
        sample: str
        ) -> Dict[str,Any]:
        """ Trimgalore trimming """

        s_time: float = time.time()

        cols: List[str] = [
            'status', 'note', 'sample_id', 'cohort_id', 'R1_trim_size',
            'R2_trim_size', 'fastqtrim_time_in_sec'
        ]
        
        o_dpath: str = os.path.join(fastq_dpath, sample)
        
        cmd: str = "trim_galore " \
            "--paired " \
            "--no_report_file " \
            "--gzip " \
            f"--cores {n_cores} " \
            f"-o {o_dpath} " \
            f"{os.path.join(o_dpath, f'{sample}_R1.fq.gz')} " \
            f"{os.path.join(o_dpath, f'{sample}_R2.fq.gz')}"

        try:
            log: str = subprocess.check_output(
                cmd, stderr=subprocess.STDOUT, shell=True)
            with open(
                os.path.join(o_dpath, f'{sample}_trimgalore.log'), 'w') as f:
                f.write(log.decode('utf-8'))
            del log
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"command '{e.cmd}' return with error " \
                f"(code {e.returncode}): {e.output}")

        for file in os.listdir(o_dpath):
            if file.endswith(f'{sample}_R1_val_1.fq.gz'):
                os.rename(
                    os.path.join(o_dpath, file),
                    os.path.join(o_dpath, f'{sample}_R1.fq.gz')
                )  
            elif file.endswith(f'{sample}_R2_val_2.fq.gz'):
                os.rename(
                    os.path.join(o_dpath, file),
                    os.path.join(o_dpath, f'{sample}_R2.fq.gz')
                )

        r1: int = 0
        r2: int = 0
        r1_fpath: str = os.path.join(o_dpath, f'{sample}_R1.fq.gz')
        if os.path.isfile(r1_fpath):
            finfo: os.stat_result = os.stat(r1_fpath)
            r1: int = finfo.st_size
            del finfo

        r2_fpath: str = os.path.join(o_dpath, f'{sample}_R2.fq.gz')
        if os.path.isfile(r2_fpath):
            finfo: os.stat_result = os.stat(r2_fpath)
            r2: int = finfo.st_size
            del finfo
            
        d: Dict[str,Any] = {
            'status': 'Ok' if r1 and r2 else 'Err',
            'note': '' if r1 and r2 else 'Trimgalore error',
            'sample_id': sample,
            'cohort_id': sample_info[sample]['cohort_id'],
            'R1_trim_size': r1,
            'R2_trim_size': r2,
            'fastqtrim_time_in_sec': round(time.time()-s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d


    def check_strandedness(
        self,
        sample_info: Dict[str,Any],
        fastq_dpath: str,
        gtf_fpath: str,
        transcript_fpath: str,
        sample: str
        ) -> Dict[str,Any]:
        """ Strandedness """

        s_time: float = time.time()

        cols: List[str] = [
            'status', 'note', 'sample_id', 'cohort_id', 'strandedness',
            'strandedness_time_in_sec'
        ]
        
        o_dpath: str = os.path.join(fastq_dpath, sample)
        strandedness: str = ''

        cmd: str = "check_strandedness " \
            f"--gtf {gtf_fpath} " \
            f"--transcripts {transcript_fpath} " \
            "-n 1000000 " \
            f"--reads_1 {os.path.join(o_dpath, f'{sample}_R1.fq.gz')} " \
            f"--reads_2 {os.path.join(o_dpath, f'{sample}_R2.fq.gz')} " \
            f"--kallisto_index {sample}.kallisto_idx"

        try:
            log: str = subprocess.check_output(
                cmd, stderr=subprocess.STDOUT, shell=True)
            with open(
                os.path.join(o_dpath, f'{sample}_strandedness.log'), 'w') as f:
                f.write(log.decode('utf-8'))

                rows: List[str] = log.decode('utf-8').strip().split('\n')

                if rows:
                    strandedness = \
                        rows[-1].strip().replace('Data is likely ', '')
                del rows
            del log
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"command '{e.cmd}' return with error " \
                f"(code {e.returncode}): {e.output}")

        kallisto_idx: str = f'{sample}.kallisto_idx'
        if os.path.isfile(kallisto_idx):
            self.delete_file_or_dir(kallisto_idx)

        stranded_test: str = f'stranded_test_{sample}_R1'
        if os.path.isdir(stranded_test):
            self.delete_file_or_dir(stranded_test)

        d: Dict[str,Any] = {
            'status': 'Ok',
            'note': '',
            'sample_id': sample,
            'cohort_id': sample_info[sample]['cohort_id'],
            'strandedness': strandedness,
            'strandedness_time_in_sec': round(time.time()-s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d


    def telescope_input(
        self,
        sample_info: Dict[str,Any],
        fastq_dpath: str,
        output_dpath: str,
        n_cores: int,
        genome_fpath: str,
        scracth_dpath: str,
        sample: str
        ) -> pd.DataFrame:
        """ Input preparation for TELESCOPE """

        s_time: float = time.time()

        cols: List[str] = [
            'status', 'note', 'sample_id', 'cohort_id', 'bowties2',
            'bowties2_sam_size', 'bowties2_time_in_sec', 'samtools_sam_to_bam',
            'samtools_bam_size', 'samtools_sam_to_bam_time_in_sec',
            'samtools_bai', 'samtools_bai_time_in_sec'
        ]
            
        telescope_dpath: str = os.path.join(output_dpath, sample)
        if os.path.isdir(telescope_dpath):
            self.delete_file_or_dir(telescope_dpath)        
        self.create_dir(telescope_dpath)
        
        cmd: str = "bowtie2 " \
            "-k 100 " \
            f"-p {n_cores} " \
            f"-x {genome_fpath} " \
            f"-1 {os.path.join(fastq_dpath, sample, f'{sample}_R1.fq.gz')} " \
            f"-2 {os.path.join(fastq_dpath, sample, f'{sample}_R2.fq.gz')} " \
            f"-S {os.path.join(telescope_dpath, f'{sample}.sam')}"
        
        log: str = subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, shell=True)
        
        with open(os.path.join(
            telescope_dpath, f'{sample}_bowties2.log'), 'w') as f:
            f.write(log.decode('utf-8'))
        del log
        
        sam_fpath: str = os.path.join(telescope_dpath, f'{sample}.sam')
        if not os.path.isfile(sam_fpath):
            d: Dict[str,Any] = {
                'status': 'Err',
                'note': 'No SAM file found',
                'sample_id': sample,
                'cohort_id': sample_info[sample]['cohort_id'],
                'bowties2': 'Err',
                'bowties2_sam_size': 0,
                'bowties2_time_in_sec': round(time.time()-s_time, 2),
                'samtools_sam_to_bam': '',
                'samtools_bam_size': 0,
                'samtools_sam_to_bam_time_in_sec': .0,
                'samtools_bai': '',
                'samtools_bai_time_in_sec': .0
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d

        finfo: os.stat_result = os.stat(sam_fpath)
        bowties2_sam_size: int = finfo.st_size
        bowties2_sam_time: float = round(time.time()-s_time, 2)
        del finfo
        
        r1_fpath: str = os.path.join(fastq_dpath,  sample, f'{sample}_R1.fq.gz')  
        if os.path.isfile(r1_fpath):
            self.delete_file_or_dir(r1_fpath)
            
        r2_fpath: str = os.path.join(fastq_dpath, sample, f'{sample}_R2.fq.gz')  
        if os.path.isfile(r2_fpath):
            self.delete_file_or_dir(r2_fpath)
        
        s_time: float = time.time()
        
        cmd: str = "samtools " \
            "view " \
            f"-@ {n_cores} " \
            f"-bS {os.path.join(telescope_dpath, f'{sample}.sam')} " \
            "| " \
            "samtools " \
            "sort -o " \
            f"{os.path.join(telescope_dpath, f'{sample}.bam')}"
        
        log: str = subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, shell=True)
        
        with open(os.path.join(
            telescope_dpath, f'{sample}_sam_to_sorted_bam.log'), 'w') as f:
            f.write(log.decode('utf-8'))
        del log
        
        bam_fpath: str = os.path.join(telescope_dpath, f'{sample}.bam')
        if not os.path.isfile(bam_fpath):
            d: Dict[str,Any] = {
                'status': 'Err',
                'note': 'No BAM file found',
                'sample_id': sample,
                'cohort_id': sample_info[sample]['cohort_id'],
                'bowties2': 'Ok',
                'bowties2_sam_size': bowties2_sam_size,
                'bowties2_time_in_sec': bowties2_sam_time,
                'samtools_sam_to_bam': 'Err',
                'samtools_bam_size': 0,
                'samtools_sam_to_bam_time_in_sec': \
                    round(time.time()-s_time, 2),
                'samtools_bai': '',
                'samtools_bai_time_in_sec': .0
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d
        
        finfo: os.stat_result = os.stat(bam_fpath)
        samtools_bam_size: int = finfo.st_size
        samtools_bam_time: float = round(time.time()-s_time, 2)
        del finfo
        
        sam_fpath: str = os.path.join(telescope_dpath, f'{sample}.sam') 
        if os.path.isfile(sam_fpath):
            self.delete_file_or_dir(sam_fpath)
        
        s_time: float = time.time()
        
        cmd: str = "samtools " \
            "index " \
            f"-@ {n_cores} " \
            f"{os.path.join(telescope_dpath, f'{sample}.bam')} " \
            f"{os.path.join(telescope_dpath, f'{sample}.bam.bai')}"
        
        log: str = subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, shell=True)
        
        with open(os.path.join(
            telescope_dpath, f'{sample}_samtools_bai.log'), 'w') as f:
            f.write(log.decode('utf-8'))
        del log
        
        bai_fpath: str = os.path.join(telescope_dpath, f'{sample}.bam.bai')
        if not os.path.isfile(bai_fpath):
            d: Dict[str,Any] = {
                'status': 'Err',
                'note': 'No BAM.BAI file found',
                'sample_id': sample,
                'cohort_id': sample_info[sample]['cohort_id'],
                'bowties2': 'Ok',
                'bowties2_sam_size': bowties2_sam_size,
                'bowties2_time_in_sec': bowties2_sam_time,
                'samtools_sam_to_bam': 'Ok',
                'samtools_bam_size': samtools_bam_size,
                'samtools_sam_to_bam_time_in_sec': samtools_bam_time,
                'samtools_bai': 'Err',
                'samtools_bai_time_in_sec': round(time.time()-s_time, 2)
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d

        src_fpath: str = os.path.join(telescope_dpath, f'{sample}.bam')
        dst_fpath: str = os.path.join(scracth_dpath, f'{sample}.bam')
        move(src_fpath, dst_fpath)

        src_fpath: str = os.path.join(telescope_dpath, f'{sample}.bam.bai')
        dst_fpath: str = os.path.join(scracth_dpath, f'{sample}.bam.bai')
        move(src_fpath, dst_fpath)

        src_fpath: str = os.path.join(telescope_dpath, f'{sample}.bam')
        if os.path.isfile(src_fpath):
            self.delete_file_or_dir(src_fpath)
        
        src_fpath: str = os.path.join(telescope_dpath, f'{sample}.bam.bai')
        if os.path.isfile(src_fpath):
            self.delete_file_or_dir(src_fpath)

        d: Dict[str,Any] = {
            'status': 'Ok',
            'note': '',
            'sample_id': sample,
            'cohort_id': sample_info[sample]['cohort_id'],
            'bowties2': 'Ok',
            'bowties2_sam_size': bowties2_sam_size,
            'bowties2_time_in_sec': bowties2_sam_time,
            'samtools_sam_to_bam': 'Ok',
            'samtools_bam_size': samtools_bam_size,
            'samtools_sam_to_bam_time_in_sec': samtools_bam_time,
            'samtools_bai': 'Ok',
            'samtools_bai_time_in_sec': round(time.time()-s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d


    def telescope(
        self,
        sample_info: Dict[str,Any],
        output_dpath: str,
        n_cores: int,
        herv_gtf_fpath: str,
        sample: str
        ) -> pd.DataFrame:
        """ TELESCOPE alignment """
        
        s_time: float = time.time()

        cols: List[str] = [
            'status', 'note', 'sample_id', 'cohort_id', 'telescope',
            'telescope_time_in_sec'
        ]

        """
        if sample_info[sample]['strandedness'] == '' or \
            sample_info[sample]['strandedness'] not in \
                self.telescope_strandedness.keys():
            d: Dict[str,Any] = {
                'status': 'Err',
                'note': "Wrong strandedness: " \
                    f"{sample_info[sample]['strandedness']}",
                'sample_id': sample,
                'cohort_id': sample_info[sample]['cohort_id'],
                'telescope': '',
                'telescope_time_in_sec': round(time.time()-s_time, 2)
            }
            print('\t'.join([str(d[c]) for c in cols]))
            return d
        """

        telescope_dpath: str = os.path.join(output_dpath, sample)
        tmp_dpath: str = os.path.join(telescope_dpath, 'tmp')
        if os.path.isdir(tmp_dpath):
            self.delete_file_or_dir(tmp_dpath)        
        self.create_dir(tmp_dpath)
        
        """
        strandedness: str = self.telescope_strandedness[
            sample_info[sample]['strandedness']]

        cmd: str = "telescope " \
            "bulk assign " \
            f"--ncpu {n_cores} " \
            f"--tempdir {tmp_dpath} " \
            f"--outdir {telescope_dpath} " \
            f"--exp_tag {sample} " \
            f"--stranded_mode {strandedness} " \
            f"{os.path.join(telescope_dpath, f'{sample}.bam')} " \
            f"{herv_gtf_fpath}"
        """

        cmd: str = "telescope " \
            "assign " \
            f"{os.path.join(telescope_dpath, f'{sample}.bam')} " \
            f"{herv_gtf_fpath} " \
            f"--ncpu {n_cores} " \
            f"--tempdir {tmp_dpath} " \
            f"--outdir {telescope_dpath} " \
            f"--exp_tag {sample}"         

        log: str = subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, shell=True)

        with open(os.path.join(
            telescope_dpath, f'{sample}_telescope.log'), 'w') as f:
            f.write(log.decode('utf-8'))
        del log
        
        tbam_fpath: str = os.path.join(telescope_dpath, f'{sample}.bam')  
        if os.path.isfile(tbam_fpath):
            self.delete_file_or_dir(tbam_fpath)
        
        tbai_fpath: str = os.path.join(telescope_dpath, f'{sample}.bam.bai')
        if os.path.isfile(tbai_fpath):
            self.delete_file_or_dir(tbai_fpath)
        
        tnpz_fpath: str = os.path.join(
            telescope_dpath, f'{sample}-checkpoint.npz')
        if os.path.isfile(tnpz_fpath):
            self.delete_file_or_dir(tnpz_fpath)
        
        ttmp_dpath: str = os.path.join(telescope_dpath, 'tmp')
        if os.path.isdir(ttmp_dpath):
            self.delete_file_or_dir(ttmp_dpath)

        d: Dict[str,Any] = {
            'status': 'Ok',
            'note': '',
            'sample_id': sample,
            'cohort_id': sample_info[sample]['cohort_id'],
            'telescope': 'Ok',
            'telescope_time_in_sec': round(time.time()-s_time, 2)
        }
        print('\t'.join([str(d[c]) for c in cols]))
        return d
        

    def samples_chunking(self, samples, chunk_size):
        return [
            samples[i:i+chunk_size] for i in range(0, len(samples), chunk_size)]


    def run(self) -> None:
        if 'BOWTIE' not in self.workflows and \
            'TELESCOPE' not in self.workflows:
            return

        # Initialization
        o_path: str = self.cmd_obj.out_dir
        scratch_path: str = self.cmd_obj.scratch_dir
        
        if os.path.isdir(o_path):
            self.delete_file_or_dir(o_path)
        self.create_dir(o_path)

        for workflow in self.workflows:
            w_path: str = os.path.join(o_path, workflow)
            if os.path.isdir(w_path):
                self.delete_file_or_dir(w_path)
            self.create_dir(w_path)

        logs_map: Dict[str,Any] = {}

        # =====================================================================
        # BOWTIE
        # =====================================================================
        if 'BOWTIE' in self.workflows:
            # Download samples from GDC
            sample_info: Dict[str,Any] = self.prepare_samples_info()

            bam_dpath: str = 'BAM'
            if os.path.isdir(bam_dpath):
                self.delete_file_or_dir(bam_dpath)
            self.create_dir(bam_dpath)
            fastq_dpath: str = 'FASTQ'
            if os.path.isdir(fastq_dpath):
                self.delete_file_or_dir(fastq_dpath)
            self.create_dir(fastq_dpath)

            # =================================================================
            # Samples download
            # =================================================================
            downloaded_samples: List[str] = []
            
            for sample,info in sample_info.items():
                cmd: str = f"{self.cmd_obj.gdc_cl} " \
                    "download " \
                    "--no-related-files " \
                    "--no-annotations " \
                    f"-t {self.cmd_obj.gdc_tk} " \
                    f"-d {bam_dpath} " \
                    f"{info['file_id']}"

                log: str = subprocess.check_output(
                    cmd, stderr=subprocess.STDOUT, shell=True
                ).decode('utf-8')
                
                del log
                
                src_fpath: str = os.path.join(
                    bam_dpath, info['file_id'], info['file_name']
                )
                dst_fpath: str = os.path.join(bam_dpath, f'{sample}.bam')

                status: str = 'Err'
                if os.path.isfile(src_fpath):
                    move(src_fpath, dst_fpath)
                    self.delete_file_or_dir(
                        os.path.join(bam_dpath, info['file_id']))
                    downloaded_samples.append(sample)
                    status = 'Ok'

                if sample not in logs_map:
                    logs_map[sample] = {
                        'status': status,
                        'note': \
                            '' if status=='Ok' else 'No BAM file downloaded',
                        'sample_id': sample,
                        'cohort_id': info['cohort_id']
                    }

            """
            # !!! Only for local mapping !!!
            for sample,info in sample_info.items():
                dst_fpath: str = os.path.join(bam_dpath, f'{sample}.bam')
                status: str = 'Err'
                if os.path.isfile(dst_fpath):
                    downloaded_samples.append(sample)
                    status = 'Ok'
                
                if sample not in logs_map:
                    logs_map[sample] = {
                        'status': status,
                        'note': \
                            '' if status=='Ok' else 'No BAM file downloaded',
                        'sample_id': sample,
                        'cohort_id': info['cohort_id']
                    }

            """

            init_samples: List[str] = sorted(list(
                set(sample_info.keys()).intersection(set(downloaded_samples))
            ))
            del downloaded_samples

            if not init_samples:
                self.save_logs(logs_map, o_path)
                del logs_map
                return

            # =================================================================
            # Samples check
            # =================================================================
            sample_to_process: List[str] = []
            try:
                n_workers: int = int(math.floor(
                    self.n_cores/self.sample_check_n_cores))

                func = partial(
                    self.sample_check,
                    sample_info,
                    bam_dpath,
                    self.sample_check_n_cores
                )

                s_time: float = time.time()
                    
                with ProcessPoolExecutor(max_workers=n_workers) as executor:
                    for result in executor.map(
                        func, init_samples, chunksize=1):
                        d: Dict[str,Any] = \
                            {**logs_map[result['sample_id']], **result}
                        logs_map[result['sample_id']] = d
                        del d
                        
                print("BAM samples validation (wall-time): " \
                    f"{round(time.time()-s_time, 2)} (s)")

                del s_time, init_samples

                sample_to_process = [
                    v['sample_id'] for _,v in logs_map.items()
                    if v['status']=='Ok'
                ]
            except Exception:
                traceback.print_exc()
            
            if not sample_to_process:
                self.save_logs(logs_map, o_path)
                del logs_map
                return

            # =================================================================
            # BAM to FASTQ
            # =================================================================
            sample_to_process1: List[str] = []

            for samples_chunk_samples in self.samples_chunking(
                sample_to_process, self.fastq_samples_limit):

                try:
                    n_workers: int = self.n_cores
                    func = partial(
                        self.bam_to_fastq,
                        sample_info,
                        bam_dpath,
                        fastq_dpath
                    )
                    s_time: float = time.time()
                        
                    with ProcessPoolExecutor(
                        max_workers=n_workers) as executor:
                        for result in executor.map(
                            func, samples_chunk_samples, chunksize=1):
                            d: Dict[str,Any] = \
                                {**logs_map[result['sample_id']], **result}
                            logs_map[result['sample_id']] = d
                            del d
                            
                    print("BAM to FastQ (wall-time): " \
                        f"{round(time.time()-s_time, 2)} (s)")

                    del s_time

                    sample_to_process1.extend([v['sample_id'] 
                        for _,v in logs_map.items() if v['status']=='Ok'])
                except Exception:
                    traceback.print_exc()

            sample_to_process1 = list(set(sample_to_process1))

            if not sample_to_process1:
                self.save_logs(logs_map, o_path)
                del logs_map
                return

            # =================================================================
            # FASTQ trimming
            # =================================================================
            sample_to_process2: List[str] = []
            try:           
                n_workers: int = int(math.floor(
                    self.n_cores/self.trimgalore_n_cores))

                func = partial(
                    self.fastq_trimming,
                    sample_info,
                    fastq_dpath,
                    self.trimgalore_n_cores
                )
                s_time: float = time.time()
                    
                with ProcessPoolExecutor(max_workers=n_workers) as executor:
                    for result in executor.map(
                        func, sample_to_process1, chunksize=1):
                        d: Dict[str,Any] = \
                            {**logs_map[result['sample_id']], **result}
                        logs_map[result['sample_id']] = d
                        del d
                
                print("FastQ samples trimming (wall-time): " \
                    f"{round(time.time()-s_time, 2)} (s)")

                del s_time, sample_to_process1

                sample_to_process2 = [v['sample_id'] 
                    for _,v in logs_map.items() if v['status']=='Ok']
        
            except Exception:
                traceback.print_exc()
                self.save_logs(logs_map, o_path)

            if not sample_to_process2:
                self.save_logs(logs_map, o_path)
                del logs_map
                return

            # =================================================================
            # Strandedness
            # =================================================================
            sample_to_process3: List[str] = []
            try:
                func = partial(
                    self.check_strandedness,
                    sample_info,
                    fastq_dpath,
                    self.gtf_fpath,
                    self.transcript_fpath
                )
                s_time: float = time.time()
                    
                with ProcessPoolExecutor(
                    max_workers=self.check_strand_n_workers) as executor:
                    for result in executor.map(
                        func, sample_to_process2, chunksize=1):
                        d: Dict[str,Any] = \
                            {**logs_map[result['sample_id']], **result}
                        logs_map[result['sample_id']] = d

                        if result['strandedness'] in self.telescope_strandedness:
                            sample_to_process3.append(result['sample_id'])

                            sample_info[result['sample_id']]['telescope_strandedness'] = \
                                self.telescope_strandedness[result['strandedness']]
                        else:
                            print("Strandedness not supported for sample:", sample)
                            # Remove sample files
                            b_fpath: str = os.path.join(bam_dpath, f'{sample}.bam')
                            if os.path.isfile(b_fpath):
                                self.delete_file_or_dir(b_fpath)

                            r1_fpath: str = os.path.join(
                                fastq_dpath,  sample, f'{sample}_R1.fq.gz')  
                            if os.path.isfile(r1_fpath):
                                self.delete_file_or_dir(r1_fpath)
                                
                            r2_fpath: str = os.path.join(
                                fastq_dpath, sample, f'{sample}_R2.fq.gz')  
                            if os.path.isfile(r2_fpath):
                                self.delete_file_or_dir(r2_fpath)
                        del d
                
                print("Samples strandedness check (wall-time): " \
                    f"{round(time.time()-s_time, 2)} (s)")

                del s_time
            except Exception:
                traceback.print_exc()
                self.save_logs(logs_map, o_path)

            if not sample_to_process3:
                self.save_logs(logs_map, o_path)
                del logs_map
                return

            # BOWTIE
            try:
                s_time: float = time.time()

                for sample in sample_to_process3:
                    result: Dict[str,Any] = self.telescope_input(
                        sample_info,
                        fastq_dpath,
                        os.path.join(o_path, 'BOWTIE'),
                        self.n_cores,
                        self.bowtiew2_genome_dpath,
                        scratch_path,
                        sample
                    )
                    d: Dict[str,Any] = \
                        {**logs_map[result['sample_id']], **result}
                    logs_map[result['sample_id']] = d
                    del d
                    
                print("Bowtie2 processing (wall-time): " \
                    f"{round(time.time()-s_time, 2)} (s)")
                del s_time

            except Exception:
                traceback.print_exc()
                self.save_logs(logs_map, o_path)
                pass

            del sample_info

        # =====================================================================
        # TELESCOPE
        # =====================================================================
        if 'TELESCOPE' in self.workflows:

            sample_info: Dict[str,Any] = self.prepare_telescope_samples_info()
            init_samples: List[str] = []

            for sample,info in sample_info.items():
                bam_fpath: str = os.path.join(scratch_path, f'{sample}.bam')
                bai_fpath: str = \
                    os.path.join(scratch_path, f'{sample}.bam.bai')
                
                if os.path.isfile(bam_fpath) and \
                    os.path.isfile(bai_fpath):

                    tele_dpath: str = os.path.join(o_path, 'TELESCOPE', sample)
                    if os.path.isdir(tele_dpath):
                        self.delete_file_or_dir(tele_dpath)        
                    self.create_dir(tele_dpath)

                    src_fpath: str = \
                        os.path.join(scratch_path, f'{sample}.bam')
                    dst_fpath: str = \
                        os.path.join(tele_dpath, f'{sample}.bam')
                    #move(src_fpath, dst_fpath)
                    copyfile(src_fpath, dst_fpath)

                    src_fpath: str = \
                        os.path.join(scratch_path, f'{sample}.bam.bai')
                    dst_fpath: str = \
                        os.path.join(tele_dpath, f'{sample}.bam.bai')
                    #move(src_fpath, dst_fpath)
                    copyfile(src_fpath, dst_fpath)
                    
                    init_samples.append(sample)

                    logs_map[sample] = {
                        'status': 'Ok',
                        'note': '',
                        'sample_id': sample,
                        'cohort_id': info['cohort_id']
                    }

            if not init_samples:
                self.save_logs(logs_map, o_path)
                del logs_map
                return

            s_time: float = time.time()

            for samples_chunk_samples in self.samples_chunking(
                init_samples, self.bam_samples_limit):

                try:
                    func = partial(
                        self.telescope,
                        sample_info,
                        os.path.join(o_path, 'TELESCOPE'),
                        self.telescope_n_cores,
                        self.herv_gtf_fpath
                    )

                    with ProcessPoolExecutor(
                        max_workers=self.bam_samples_limit) as executor:
                        for result in executor.map(func, samples_chunk_samples, 
                            chunksize=1):
                            d: Dict[str,Any] = \
                                {**logs_map[result['sample_id']], **result}
                            logs_map[result['sample_id']] = d
                            del d
                    
                except Exception:
                    traceback.print_exc()
                    self.save_logs(logs_map, o_path)
                    pass

            print("TELESCOPE (wall-time): " \
                f"{round(time.time()-s_time, 2)} (s)")

            del s_time, init_samples, sample_info

        self.save_logs(logs_map, o_path)
        del logs_map


    def save_logs(self, runs_log, dpath):
        df: pd.DataFrame = pd.DataFrame([
            v for _,v in runs_log.items()
        ])
        df.to_csv(
            os.path.join(dpath, 'logs.tsv'),
            sep='\t',
            encoding='utf-8',
            index=False
        )


    def delete_file_or_dir(self, path: str) -> None:
        """
        Deletes a file or directory
        
        Arguments:
        -----------
        - path: A file or directory path. It could either be relative or absolute.
        """
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
        """
        Create a directory
        
        Arguments:
        -----------
        - path: A directory path. It could either be relative or absolute.
        """
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
                
    
    def print_error(self, *args, **kwargs):
        print(*args, file = sys.stderr, **kwargs)
        
    def print_error_and_abort(self, *args, **kwargs):
        print(*args, file = sys.stderr, **kwargs)
        sys.exit(1)
