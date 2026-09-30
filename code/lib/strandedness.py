#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Library strandedness -> Telescope --stranded_mode.

Two ways to measure it, one way to call it:

* BAM input  : RSeQC infer_experiment.py on the sample's own alignments
               (seconds; no index needed).
* FASTQ input: how_are_we_stranded_here's check_strandedness, which pseudo-
               aligns a read sample with kallisto and runs the same
               infer_experiment.py on the result. The kallisto index is built
               once and reused (it is the slow part: ~8 min per build).

Both end in the same two fractions -- reads explained by "1++,1--,2+-,2-+"
(FR, read 1 sense) and by "1+-,1-+,2++,2--" (RF, read 1 antisense) -- and
call_mode() turns them into a Telescope value with one set of thresholds.

    >= 0.8 of informative reads RF  -> RF       (single-end: R)
    >= 0.8 of informative reads FR  -> FR       (single-end: F)
    0.35-0.65                       -> None     (unstranded)
    anything else                   -> UNCLEAR  (mixed/weak; Telescope gets None)

UNCLEAR typically means several libraries with different protocols merged in
one BAM, or a weakly stranded library. None is the safe Telescope value for it:
it counts both strands and is never wrong about orientation.

Runs on Python 3.6 (the immuno_telescope env). Also usable as a CLI:

    python strandedness.py bam   --bam X.bam --gtf gencode.gtf [--bed genes.bed12]
    python strandedness.py fastq --fq1 R1.fq.gz --fq2 R2.fq.gz --gtf gencode.gtf \
                                 --transcripts gencode.transcripts.fa --kallisto-index idx
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

TELESCOPE_MODES = ("None", "RF", "FR", "R", "F")

FR_PATTERNS = (r'"1\+\+,1--,2\+-,2-\+"', r'"\+\+,--"')
RF_PATTERNS = (r'"1\+-,1-\+,2\+\+,2--"', r'"\+-,-\+"')


def _run(cmd, cwd=None):
    """-> (returncode, combined stdout+stderr) ; Python 3.6 compatible."""
    p = subprocess.Popen(cmd, shell=isinstance(cmd, str), cwd=cwd,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         universal_newlines=True, executable='/bin/bash' if isinstance(cmd, str) else None)
    out, _ = p.communicate()
    return p.returncode, out


# ---------------------------------------------------------------------------
# annotation
# ---------------------------------------------------------------------------

def bed_from_gtf(gtf, out_bed):
    """Protein-coding gene spans from a GENCODE GTF as BED12 (what RSeQC reads).

    Written atomically, and reused if it already exists and is newer than the
    GTF. Only chr1-22/X/Y are kept: the strand of a read on a scaffold tells us
    nothing extra and those contigs are not in the bowtie2 index anyway.
    """
    if os.path.isfile(out_bed) and os.path.getmtime(out_bed) >= os.path.getmtime(gtf):
        return out_bed
    rx_id = re.compile(r'gene_id "([^"]+)"')
    rows = []
    with open(gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'gene' or 'gene_type "protein_coding"' not in p[8]:
                continue
            if not re.match(r'^chr([0-9]+|X|Y)$', p[0]):
                continue
            m = rx_id.search(p[8])
            s, e = int(p[3]) - 1, int(p[4])
            rows.append((p[0], s, e, m.group(1) if m else '.', p[6]))
    if not rows:
        raise ValueError('no protein_coding genes on chr1-22/X/Y in %s -- not a GENCODE GTF?' % gtf)
    rows.sort(key=lambda r: (r[0], r[1]))
    os.makedirs(os.path.dirname(os.path.abspath(out_bed)), exist_ok=True)
    tmp = out_bed + '.tmp%d' % os.getpid()
    with open(tmp, 'w') as f:
        for c, s, e, gid, strand in rows:
            f.write('%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0\t1\t%d,\t0,\n' % (c, s, e, gid, strand, s, e, e - s))
    os.replace(tmp, out_bed)
    return out_bed


def bed_contigs(bed):
    with open(bed) as f:
        return {ln.split('\t', 1)[0] for ln in f if ln.strip()}


# ---------------------------------------------------------------------------
# parsing + calling
# ---------------------------------------------------------------------------

def parse_infer_experiment(text):
    """-> dict(layout, frac_fr, frac_rf, frac_undetermined, n_reads) from
    infer_experiment.py output (also found inside check_strandedness output)."""
    def frac(patterns):
        for pat in patterns:
            m = re.search(pat + r'\s*:\s*([0-9.]+)', text)
            if m:
                return float(m.group(1))
        return None
    und = re.search(r'failed to determine\s*:\s*([0-9.]+)', text)
    n = re.search(r'Total (\d+) usable reads', text)
    layout = 'PE' if 'PairEnd' in text else 'SE' if 'SingleEnd' in text else ''
    return {
        'layout': layout,
        'frac_fr': frac(FR_PATTERNS),
        'frac_rf': frac(RF_PATTERNS),
        'frac_undetermined': float(und.group(1)) if und else None,
        'n_reads': int(n.group(1)) if n else None,
    }


def call_mode(layout, frac_fr, frac_rf):
    """-> (stranded_mode, telescope_mode, note)."""
    if frac_fr is None or frac_rf is None:
        return 'UNCLEAR', 'None', 'could not parse strand fractions'
    tot = frac_fr + frac_rf
    if tot < 0.5:
        return 'UNCLEAR', 'None', 'too few informative reads (%.2f)' % tot
    s = frac_fr / tot
    pe = layout != 'SE'
    if s <= 0.2:
        mode, note = ('RF' if pe else 'R'), 'stranded, read 1 antisense'
    elif s >= 0.8:
        mode, note = ('FR' if pe else 'F'), 'stranded, read 1 sense'
    elif 0.35 <= s <= 0.65:
        mode, note = 'None', 'unstranded'
    else:
        return 'UNCLEAR', 'None', 'partly stranded (FR share %.2f); Telescope gets None' % s
    return mode, mode, note


# ---------------------------------------------------------------------------
# measuring
# ---------------------------------------------------------------------------

def read_groups(bam, samtools='samtools'):
    """-> [(center, library), ...] from @RG lines (empty if none)."""
    rc, out = _run([samtools, 'view', '-H', bam])
    if rc != 0:
        raise RuntimeError('samtools view -H failed on %s: %s' % (bam, out[-400:]))
    rgs = []
    for line in out.splitlines():
        if not line.startswith('@RG'):
            continue
        tags = dict(t.split(':', 1) for t in line.split('\t')[1:] if ':' in t)
        rgs.append((tags.get('CN', ''), tags.get('LB', '')))
    return rgs


def infer_bam(bam, bed, n_reads=1000000, infer_exe='infer_experiment.py', samtools='samtools'):
    """Strandedness of an aligned BAM. Contig names must match the BED."""
    rc, hdr = _run([samtools, 'view', '-H', bam])
    if rc != 0:
        raise RuntimeError('cannot read BAM header of %s: %s' % (bam, hdr[-400:]))
    sq = set(re.findall(r'^@SQ\tSN:(\S+)', hdr, re.M))
    shared = sq & bed_contigs(bed)
    if not shared:
        raise ValueError('%s shares no contig names with %s (BAM has e.g. %s) -- use an '
                         'annotation with the same naming (chr1 vs 1)' % (bam, bed, sorted(sq)[:3]))
    rc, out = _run([infer_exe, '-r', bed, '-i', bam, '-s', str(n_reads)])
    if rc != 0 or 'Fraction' not in out:
        raise RuntimeError('infer_experiment.py failed on %s: %s' % (bam, out[-600:]))
    res = parse_infer_experiment(out)
    res['stranded_mode'], res['telescope_mode'], res['note'] = call_mode(
        res['layout'], res['frac_fr'], res['frac_rf'])
    rgs = read_groups(bam, samtools)
    centers = sorted({c for c, _ in rgs if c})
    res['centers'] = ','.join(centers)
    res['n_libraries'] = len({l for _, l in rgs}) if rgs else 0
    if res['stranded_mode'] == 'UNCLEAR' and len(centers) > 1:
        res['note'] += '; libraries from several centers merged in one BAM'
    res['method'] = 'rseqc_bam'
    res['raw'] = out
    return res


def ensure_kallisto_index(transcripts, index_path, kallisto='kallisto'):
    """Build the kallisto index once; safe when many jobs ask at the same time.

    Same recipe as check_strandedness uses when it builds one itself: GENCODE
    headers ('ENST...|ENSG...|...') are split on '|' so the transcript name
    kallisto keeps is the bare ENST id that the GTF-derived BED uses.
    """
    import fcntl
    if os.path.isfile(index_path) and os.path.getsize(index_path) > 0:
        return index_path
    os.makedirs(os.path.dirname(os.path.abspath(index_path)) or '.', exist_ok=True)
    with open(index_path + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if os.path.isfile(index_path) and os.path.getsize(index_path) > 0:
            return index_path                      # another job built it while we waited
        tmp_fa = index_path + '.transcripts.fa.tmp%d' % os.getpid()
        tmp_idx = index_path + '.tmp%d' % os.getpid()
        try:
            rc, out = _run("sed 's/[|]/ /g' %s > %s" % (transcripts, tmp_fa))
            if rc != 0:
                raise RuntimeError('could not rewrite %s: %s' % (transcripts, out[-400:]))
            rc, out = _run([kallisto, 'index', '-i', tmp_idx, tmp_fa])
            if rc != 0 or not os.path.isfile(tmp_idx):
                raise RuntimeError('kallisto index failed: %s' % out[-800:])
            os.replace(tmp_idx, index_path)
        finally:
            for p in (tmp_fa, tmp_idx):
                if os.path.exists(p):
                    os.remove(p)
    return index_path


def infer_fastq(fq1, fq2, gtf, transcripts, kallisto_index, n_reads=1000000,
                exe='check_strandedness', workdir=None):
    """Strandedness of a FASTQ pair via kallisto pseudo-alignment.

    The kallisto index is built on first use (ensure_kallisto_index) and reused
    afterwards, instead of once per sample.
    """
    # check_strandedness runs in a scratch directory, so every path it gets must be absolute
    fq1, fq2, gtf, transcripts = (os.path.abspath(p) for p in (fq1, fq2, gtf, transcripts))
    ensure_kallisto_index(transcripts, kallisto_index)
    work = tempfile.mkdtemp(prefix='strand_', dir=workdir)
    try:
        cmd = [exe, '--gtf', gtf, '--transcripts', transcripts, '-n', str(n_reads),
               '--reads_1', fq1, '--reads_2', fq2, '--kallisto_index', os.path.abspath(kallisto_index)]
        rc, out = _run(cmd, cwd=work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if rc != 0 or 'Fraction' not in out:
        raise RuntimeError('check_strandedness failed on %s: %s' % (fq1, out[-800:]))
    res = parse_infer_experiment(out)
    res['layout'] = res['layout'] or 'PE'
    res['stranded_mode'], res['telescope_mode'], res['note'] = call_mode(
        res['layout'], res['frac_fr'], res['frac_rf'])
    res['centers'], res['n_libraries'] = '', 0
    res['method'] = 'kallisto_fastq'
    res['raw'] = out
    return res


# ---------------------------------------------------------------------------
# cli
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='kind')
    b = sp.add_parser('bam')
    b.add_argument('--bam', required=True)
    b.add_argument('--gtf', help='GENCODE GTF (a BED12 is derived and cached next to --bed)')
    b.add_argument('--bed', help='gene BED12; default <gtf>.protein_coding.genes.bed12 in $TMPDIR')
    f = sp.add_parser('fastq')
    f.add_argument('--fq1', required=True)
    f.add_argument('--fq2', required=True)
    f.add_argument('--gtf', required=True)
    f.add_argument('--transcripts', required=True)
    f.add_argument('--kallisto-index', required=True)
    for p in (b, f):
        p.add_argument('-n', '--n-reads', type=int, default=1000000)
    a = ap.parse_args(argv)
    if a.kind == 'bam':
        bed = a.bed
        if not bed:
            if not a.gtf:
                ap.error('bam needs --bed or --gtf')
            bed = os.path.join(tempfile.gettempdir(), os.path.basename(a.gtf) + '.protein_coding.genes.bed12')
        if a.gtf:
            bed_from_gtf(a.gtf, bed)
        res = infer_bam(a.bam, bed, a.n_reads)
    elif a.kind == 'fastq':
        res = infer_fastq(a.fq1, a.fq2, a.gtf, a.transcripts, a.kallisto_index, a.n_reads)
    else:
        ap.error('choose bam or fastq')
    res.pop('raw', None)
    print(json.dumps(res))
    return 0


if __name__ == '__main__':
    sys.exit(main())
