#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  24 10:29:10 2022

@author: Rosario Distefano
"""

from __future__ import print_function
from argparse import ArgumentParser
from ast import arg
import os.path
import sys
from typing import List
import json

class CmdParser():
        
    def __init__(self):
        """ Command line argument class"""
        
        # Define the AgumentParser object
        cmd_args: ArgumentParser = ArgumentParser(description="Run run.py")
        
        cmd_args.add_argument(
            '--gdc_cl', 
            type=str, 
            default='', 
            required=False,
            help='GDC client path'
        )

        cmd_args.add_argument(
            '--gdc_tk', 
            type=str, 
            default='', 
            required=False,
            help='GDC tocken file path (TXT)'
        )

        cmd_args.add_argument(
            '--gdc_mf', 
            type=str, 
            default='', 
            required=False,
            help='GDC manifest file path (TXT)'
        )

        cmd_args.add_argument(
            '--gdc_mt', 
            type=str, 
            default='', 
            required=False,
            help='GDC metadata file path (JSON)'
        )
        
        cmd_args.add_argument(
            '--gtf', 
            type=str, 
            default='', 
            required=True,
            help='gencode.v39.annotation.gtf path'
        )

        cmd_args.add_argument(
            '--genome_dna', 
            type=str, 
            default='', 
            required=False,
            help='GRCh38.primary_assembly.genome.fa path'
        )

        cmd_args.add_argument(
            '--genome', 
            type=str, 
            default='', 
            required=True,
            help='GRCh38.p13.genome.fa.gz path'
        )

        cmd_args.add_argument(
            '--transcript', 
            type=str, 
            default='', 
            required=True,
            help='gencode.v39.transcripts.fa path'
        )

        cmd_args.add_argument(
            '--herv_gtf', 
            type=str, 
            default='', 
            required=True,
            help='HG38_HERV_all_families_telescope_ann.gtf path'
        )

        cmd_args.add_argument(
            '--bowtiew2_idx', 
            type=str, 
            default='', 
            required=True,
            help='Bowtie2 genome index directory path'
        )

        cmd_args.add_argument(
            '--workflows', 
            type=str, 
            default='', 
            required=True,
            help='Workflows to perform (i.e., BOWTIE,TELESCOPE)'
        )

        cmd_args.add_argument(
            '--out_dir', 
            type=str, 
            default='', 
            required=True,
            help='Output directory path'
        )

        cmd_args.add_argument(
            '--scratch_dir', 
            type=str, 
            default='', 
            required=True,
            help='Scratch directory path'
        )

        cmd_args.add_argument(
            '--n_cores', 
            type=int, 
            default=1, 
            required=True,
            help='Number of CPU cores'
        )

        cmd_args.add_argument(
            '--sample_check_n_cores', 
            type=int, 
            default=1, 
            required=False,
            help='Number of CPU cores for sample check'
        )

        cmd_args.add_argument(
            '--trimgalore_n_cores', 
            type=int, 
            default=1, 
            required=True,
            help='Number of CPU cores for trimgalore'
        )

        cmd_args.add_argument(
            '--telescope_n_cores', 
            type=int, 
            default=1, 
            required=True,
            help='Number of CPU cores for telescope'
        )

        cmd_args.add_argument(
            '--seed', 
            type=int, 
            default=123456, 
            required=False,
            help='Seed (default: 123456)'
        )

        cmd_args.add_argument(
            '--manifest', 
            type=str, 
            default='', 
            required=False,
            help='Manifest file path'
        )

        cmd_args.add_argument(
            '--samples_dir', 
            type=str, 
            default='', 
            required=False,
            help='Samples directory path'
        )

        cmd_args.add_argument(
            '--fastq_mode',
            action='store_true',
            help='Fastq mode, read local fastq files'
        )
              
        args = cmd_args.parse_args()      
        dict_args = vars(args)
        print(json.dumps(dict_args))

        allowed_workflows: List[str] = ['BOWTIE','TELESCOPE']
        self.gdc_cl: str = args.gdc_cl
        self.gdc_tk: str = args.gdc_tk
        self.gdc_mf: str = args.gdc_mf
        self.gdc_mt: str= args.gdc_mt
        self.gtf: str= args.gtf
        self.genome_dna: str = args.genome_dna
        self.genome: str = args.genome
        self.transcript: str = args.transcript
        self.herv_gtf: str= args.herv_gtf
        self.bowtiew2_idx: str = args.bowtiew2_idx.rstrip('/')
        self.out_dir: str = args.out_dir.rstrip('/')
        self.scratch_dir: str = args.scratch_dir.rstrip('/')
        self.workflows: List[str] = [
            wk.strip() for wk in args.workflows.split(',')
            if wk.strip() in allowed_workflows
        ]
        self.n_cores: int = args.n_cores
        self.sample_check_n_cores: int = args.sample_check_n_cores
        self.trimgalore_n_cores: int = args.trimgalore_n_cores
        self.telescope_n_cores: int = args.telescope_n_cores
        self.seed: int = args.seed
        self.manifest: str = args.manifest
        self.samples_dir: str = args.samples_dir
        self.fastq_mode: bool = args.fastq_mode

        if not self.fastq_mode and not os.path.isfile(self.gdc_cl):
            self.print_error_and_abort('ERROR. No GDC client found')

        if not self.fastq_mode and (not os.path.isfile(self.gdc_tk) or 
            not self.gdc_tk.lower().endswith('.txt')):
            self.print_error_and_abort('ERROR. No GDC token found')
            
        if not self.fastq_mode and not (os.path.isfile(self.gdc_mf) or \
            not self.gdc_mf.lower().endswith('.txt')):
            self.print_error_and_abort('ERROR. No GDC manifest found')
            
        if not self.fastq_mode and (not os.path.isfile(self.gdc_mt) or \
            not self.gdc_mt.lower().endswith('.json')):
            self.print_error_and_abort('ERROR. No GDC metadata found')

        if self.fastq_mode and ("BOWTIE" in self.workflows) and not os.path.isfile(self.manifest):
            self.print_error_and_abort('ERROR. No manifest found')

        if self.fastq_mode and ("BOWTIE" in self.workflows) and not os.path.isdir(self.samples_dir):
            self.print_error_and_abort('ERROR, No samples directory found ')
            
        if not os.path.isdir(self.out_dir):
            self.print_error_and_abort('ERROR. No output directory found')

        if not os.path.isdir(self.scratch_dir):
            self.print_error_and_abort('ERROR. No scratch directory found')

        if self.n_cores <= 0:
            self.print_error_and_abort('ERROR. The analysis need at least 10 CPU cores')
    
        
    def print_error_and_abort(self, *args, **kwargs):
        print(*args, file = sys.stderr, **kwargs)
        sys.exit(1)                
    
