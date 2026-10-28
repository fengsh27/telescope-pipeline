#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Aligned BAM -> paired FASTQ.gz, so BAM samples can enter the FASTQ pipeline.

    samtools collate -O -u in.bam | samtools fastq -F <exclude> -1 R1 -2 R2 -0 other -s singletons

Why each part matters for Telescope:

* collate   -- samtools fastq pairs mates by adjacency; a coordinate-sorted BAM
               would put most mates far apart and turn them into singletons.
* -F 0x900  -- drop secondary and supplementary records, or every multi-mapped
               read is written once per alignment (duplicated input to bowtie2).
               0x200 (vendor QC-fail) is dropped too by default; those reads are
               normally not delivered in FASTQ files.
* unmapped reads are KEPT -- reads the original aligner could not place are
               exactly the ones a repeat-aware aligner may still place.
* samtools fastq writes reverse-strand reads reverse-complemented, and keeps
  the READ1/READ2 assignment, so the FASTQ reproduces the sequenced reads and
  the library strandedness measured on the BAM still applies.

The conversion is then checked by counting, because samtools exits 0 even when
it wrote nothing:

    R1 == R2 > 0,   2*R1 + singletons + other == primary reads kept by -F

Runs on Python 3.6 (the immuno_telescope env).
"""

import os
import subprocess


def _sh(cmd):
    p = subprocess.Popen(cmd, shell=True, executable='/bin/bash',
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    out, err = p.communicate()
    if p.returncode != 0:
        raise RuntimeError("command failed (rc=%d): %s\n%s" % (p.returncode, cmd, err[-1000:]))
    return out, err


def count_fastq(path, threads=4):
    out, _ = _sh("set -o pipefail; pigz -p %d -dc '%s' | wc -l" % (threads, path))
    n = int(out.strip())
    if n % 4:
        raise RuntimeError('%s has %d lines, not a multiple of 4 -- truncated FASTQ' % (path, n))
    return n // 4


def is_paired(bam, samtools='samtools', probe=10000):
    """True if the first primary records carry the paired flag."""
    # no pipefail here: head closing the pipe early is expected, not an error
    out, _ = _sh("%s view -F 0x900 '%s' | head -%d | awk '{print $2}'" % (samtools, bam, probe))
    flags = [int(x) for x in out.split()]
    if not flags:
        raise RuntimeError('%s contains no primary reads' % bam)
    return sum(1 for f in flags if f & 1) > 0.5 * len(flags)


def bam_to_fastq(bam, out_dir, sample, threads=4, exclude_flags='0xB00',
                 samtools='samtools', max_singleton_frac=0.05):
    """-> dict with R1/R2 paths and the counts that prove nothing was lost.

    Raises on single-end input, empty output, count mismatches, or more than
    max_singleton_frac of reads left without a mate.
    """
    if not os.path.isfile(bam) or os.path.getsize(bam) == 0:
        raise FileNotFoundError('missing or empty BAM: %s' % bam)
    if not is_paired(bam, samtools):
        raise ValueError('%s is single-end; the pipeline needs paired-end reads' % bam)
    os.makedirs(out_dir, exist_ok=True)
    r1 = os.path.join(out_dir, '%s_R1.fq.gz' % sample)
    r2 = os.path.join(out_dir, '%s_R2.fq.gz' % sample)
    other = os.path.join(out_dir, '%s_other.fq.gz' % sample)
    single = os.path.join(out_dir, '%s_singletons.fq.gz' % sample)
    tmp_prefix = os.path.join(out_dir, '%s.collate' % sample)
    t = max(1, threads // 2)
    _, err = _sh(
        "set -o pipefail; "
        "%s collate -@ %d -O -u '%s' '%s' | "
        "%s fastq -@ %d -F %s -1 '%s' -2 '%s' -0 '%s' -s '%s' -"
        % (samtools, t, bam, tmp_prefix, samtools, t, exclude_flags, r1, r2, other, single))

    out, _ = _sh("%s view -c -@ %d -F %s '%s'" % (samtools, threads, exclude_flags, bam))
    primary = int(out.strip())
    n1, n2 = count_fastq(r1, threads), count_fastq(r2, threads)
    n0, ns = count_fastq(other, threads), count_fastq(single, threads)
    for p in (other, single):
        os.remove(p)

    problems = []
    if n1 == 0:
        problems.append('no read pairs written')
    if n1 != n2:
        problems.append('R1 has %d reads but R2 has %d' % (n1, n2))
    if 2 * n1 + n0 + ns != primary:
        problems.append('wrote %d reads (2x%d pairs + %d singletons + %d other) but the BAM has %d '
                        'primary reads after -F %s' % (2 * n1 + n0 + ns, n1, ns, n0, primary, exclude_flags))
    if primary and (ns + n0) > max_singleton_frac * primary:
        problems.append('%d of %d reads (%.1f%%) have no mate' % (ns + n0, primary, 100.0 * (ns + n0) / primary))
    if problems:
        for p in (r1, r2):
            if os.path.exists(p):
                os.remove(p)
        raise RuntimeError('%s: BAM->FASTQ check failed: %s' % (sample, '; '.join(problems)))
    return {'R1': r1, 'R2': r2, 'pairs': n1, 'singletons_dropped': ns + n0,
            'primary_reads': primary, 'exclude_flags': exclude_flags}
