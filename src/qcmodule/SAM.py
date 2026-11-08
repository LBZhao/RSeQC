#!/usr/bin/env python

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

# -*- coding: utf-8 -*-
'''manipulate BAM/SAM file.'''
#import built-in modules
import sys
import re
import collections
import random

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.bitset_utils import *
import pysam
from qcmodule import bam_cigar
from qcmodule import readfilter
from qcmodule import sampling
from qcmodule import BED
from importlib.metadata import version

__author__ = "Liguo Wang"
__version__ = version("RSeQC")


def r_string(text):
    """Quote *text* as an R string literal (file paths may contain quotes)."""
    return '"' + str(text).replace('\\', '\\\\').replace('"', '\\"') + '"'


class NoDataFound(Exception):
    """The input has nothing to report (e.g. no paired reads, no junctions).

    Command-line tools print the message and exit with status 0, as in
    earlier releases.
    """


JUNCTION_CATEGORIES = ("annotated", "novel_combination", "partial_novel", "complete_novel")


def load_reference_introns(refgene):
    """Read introns of a BED12 gene model, keyed by upper-case chromosome.

    Returns ``{chrom: (pairs, starts, ends)}`` where *pairs* holds
    ``(start, end)`` of every annotated intron and *starts*/*ends* the
    annotated intron starts and ends (0-based, half-open).
    """
    reference = {}
    with open(refgene, 'r') as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.split()
            if len(fields) < 12:
                print("Invalid bed line (skipped):", line, end=' ', file=sys.stderr)
                continue
            chrom = fields[0].upper()
            tx_start = int(fields[1])
            sizes = [int(x) for x in fields[10].rstrip(',').split(',')]
            starts = [tx_start + int(x) for x in fields[11].rstrip(',').split(',')]
            exons = sorted(zip(starts, (st + size for st, size in zip(starts, sizes))))
            pairs, i_starts, i_ends = reference.setdefault(chrom, (set(), set(), set()))
            for (_, i_st), (i_end, _) in zip(exons[:-1], exons[1:]):
                pairs.add((i_st, i_end))
                i_starts.add(i_st)
                i_ends.add(i_end)
    return reference


def classify_junction(reference, chrom, start, end):
    """Classify an intron ``[start, end)`` against :func:`load_reference_introns`."""
    pairs, starts, ends = reference.get(chrom.upper(), ((), (), ()))
    if (start, end) in pairs:
        return "annotated"
    known_start = start in starts
    known_end = end in ends
    if known_start and known_end:
        return "novel_combination"
    if known_start or known_end:
        return "partial_novel"
    return "complete_novel"


def load_transcript_exons(refbed):
    """Read a BED12 file into ``[(chrom, tx_start, tx_end, [(start, end), ...]), ...]``."""
    transcripts = []
    with open(refbed, 'r') as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.split()
            try:
                tx_start = int(fields[1])
                sizes = [int(x) for x in fields[10].rstrip(',').split(',')]
                starts = [tx_start + int(x) for x in fields[11].rstrip(',').split(',')]
                transcripts.append((fields[0], tx_start, int(fields[2]),
                                    sorted((st, st + size) for st, size in zip(starts, sizes))))
            except (IndexError, ValueError):
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
    return transcripts


def sampling_percentages(sample_start, sample_end, sample_step):
    """Percentages at which saturation curves are evaluated (always ends at 100).

    ``range(sample_start, sample_end, sample_step)`` plus 100. Each step adds
    the reads between the previous percentage and this one, so every read is
    used exactly once by 100%. (Before RSeQC 5.1.0+lbzhao each step added the slice
    ``[pertl - sample_step, pertl)``, which skipped reads whenever
    ``sample_start`` differed from ``sample_step`` or ``sample_end`` was below
    100.)
    """
    percents = sorted(set(range(sample_start, sample_end, sample_step)) | {100})
    return [p for p in percents if 0 < p <= 100]


class ParseBAM(object):
    '''This class provides fuctions to parsing/processing/transforming SAM or BAM files. The input
    file could be either SAM or BAM format file'''
    
    multi_hit_tags=['H0','H1','H2','IH','NH']
    def __init__(self,inputFile):
        '''constructor. input could be bam or sam'''
        # pysam detects SAM/BAM/CRAM itself. Errors (missing file, not an
        # alignment file) propagate with their real message instead of a
        # second, misleading "could not open as SAM" error.
        self.samfile = pysam.AlignmentFile(inputFile)
        if len(self.samfile.header) == 0:
            raise ValueError(f"alignment file has no header section: {inputFile}")
        self.bam_format = self.samfile.is_bam or self.samfile.is_cram

    def stat (self,q_cut=30):
        '''Calculate mapping statistics'''
        R_total=0
        R_qc_fail=0
        R_duplicate=0
        R_nonprimary=0
        R_supplementary=0
        R_unmap =0
        
        R_multipleHit=0
        R_uniqHit=0 #all the following count should be sum to uniqHit
        
        R_read1=0
        R_read2=0
        R_reverse =0
        R_forward=0
        R_nonSplice=0
        R_splice=0
        R_properPair =0 
        R_pair_diff_chrom = 0
        R_mitochondria = 0

        if self.bam_format:
            print("Load BAM file ... ", end=' ', file=sys.stderr)
        else:
            print("Load SAM file ... ", end=' ', file=sys.stderr)

        try:
            while(1):
                #flag=0
                aligned_read = next(self.samfile)
                R_total +=1
                if aligned_read.is_qcfail:          #skip QC fail read
                    R_qc_fail +=1
                    continue
                if aligned_read.is_duplicate:       #skip duplicate read
                    R_duplicate +=1
                    continue
                if aligned_read.is_secondary:       #skip non primary hit
                    R_nonprimary +=1
                    continue
                if aligned_read.is_supplementary:   #skip supplementary (chimeric) hit
                    R_supplementary +=1
                    continue
                if aligned_read.is_unmapped:        #skip unmap read
                    R_unmap +=1
                    continue        
                    
                if aligned_read.mapq < q_cut:
                    R_multipleHit +=1
                    continue                        #skip multiple map read             
                if aligned_read.mapq >= q_cut:
                    R_uniqHit +=1
                    
                    if aligned_read.is_read1:
                        R_read1 +=1
                    if aligned_read.is_read2:
                        R_read2 +=1
                    if aligned_read.is_reverse:
                        R_reverse +=1
                    else:
                        R_forward +=1
                    introns = bam_cigar.fetch_intron('chr1', aligned_read.pos, aligned_read.cigar)
                    if len(introns)==0:
                        R_nonSplice +=1
                    if len(introns)>=1:
                        R_splice +=1
                    if aligned_read.is_proper_pair:
                        R_properPair +=1
                        R_read1_ref = self.samfile.getrname(aligned_read.tid)
                        R_read2_ref = self.samfile.getrname(aligned_read.rnext)
                        if R_read1_ref != R_read2_ref:
                            #print aligned_read.qname
                            R_pair_diff_chrom +=1
                
        except StopIteration:
            print("Done", file=sys.stderr)      
        #self.samfile.seek(current_pos)
                
        print("\n#==================================================", file=sys.stdout)
        print("#All numbers are READ count", file=sys.stdout)
        print("#==================================================\n", file=sys.stdout)
        rows = [
            ("total_records", "Total records:", R_total),
            None,
            ("qc_failed", "QC failed:", R_qc_fail),
            ("duplicates", "Optical/PCR duplicate:", R_duplicate),
            ("non_primary", "Non primary hits", R_nonprimary),
            ("supplementary", "Supplementary hits:", R_supplementary),
            ("unmapped", "Unmapped reads:", R_unmap),
            ("low_mapq", "mapq < mapq_cut (non-unique):", R_multipleHit),
            None,
            ("unique", "mapq >= mapq_cut (unique):", R_uniqHit),
            ("read_1", "Read-1:", R_read1),
            ("read_2", "Read-2:", R_read2),
            ("forward_strand", "Reads map to '+':", R_forward),
            ("reverse_strand", "Reads map to '-':", R_reverse),
            ("non_spliced", "Non-splice reads:", R_nonSplice),
            ("spliced", "Splice reads:", R_splice),
            ("proper_pairs", "Reads mapped in proper pairs:", R_properPair),
            ("proper_pairs_diff_chrom", "Proper-paired reads map to different chrom:", R_pair_diff_chrom),
        ]
        for row in rows:
            if row is None:
                print("\n", end='', file=sys.stdout)
            else:
                print("%-40s%d" % (row[1], row[2]), file=sys.stdout)
        return [(key, value) for key, _label, value in (r for r in rows if r is not None)]
    
    def configure_experiment(self,refbed,sample_size, q_cut = 30, seed = 0):
        """Guess the RNA-seq library layout and strandedness.

        Reads overlapping annotated transcripts are sampled across the genome
        (see :func:`qcmodule.sampling.sample_alignments`) and their mapping
        strand is compared with the strand of the overlapping gene(s).

        Returns ``[protocol, fraction_1, fraction_2, fraction_other]`` where
        *protocol* is ``"PairEnd"`` or ``"SingleEnd"`` (the layout of most
        sampled reads; a warning is printed if both occur) or ``"Mixture"``
        if no usable reads were found.
        """
        gene_ranges = {}
        regions = []
        print("Reading reference gene model " + refbed + ' ...', end=' ', file=sys.stderr)
        with open(refbed, 'r') as handle:
            for line in handle:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                fields = line.split()
                try:
                    chrom = fields[0]
                    tx_start = int(fields[1])
                    tx_end = int(fields[2])
                    strand = fields[5]
                except (IndexError, ValueError):
                    print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                    continue
                gene_ranges.setdefault(chrom, Intersecter()).insert(tx_start, tx_end, strand)
                regions.append((chrom, tx_start, tx_end))
        print("Done", file=sys.stderr)

        p_strandness = collections.defaultdict(int)
        s_strandness = collections.defaultdict(int)

        def gene_strand(read):
            """Strand(s) of genes overlapping the read, e.g. '+' or '+:-'."""
            chrom = read.reference_name
            if chrom not in gene_ranges or read.reference_end is None:
                return None
            strands = set(gene_ranges[chrom].find(read.reference_start, read.reference_end))
            return ':'.join(sorted(strands)) if strands else None

        def usable(read):
            return readfilter.passes(read, min_mapq=q_cut) and gene_strand(read) is not None

        print("Loading SAM/BAM file ... ", end=' ', file=sys.stderr)
        count = 0
        for read in sampling.sample_alignments(self.samfile, regions, sample_size, usable, seed=seed):
            map_strand = '-' if read.is_reverse else '+'
            strand_from_gene = gene_strand(read)
            if read.is_paired:
                read_id = '1' if read.is_read1 else '2'
                p_strandness[read_id + map_strand + strand_from_gene] += 1
            else:
                s_strandness[map_strand + strand_from_gene] += 1
            count += 1
        print("Finished", file=sys.stderr)
        print("Total " + str(count) + " usable reads were sampled", file=sys.stderr)

        n_paired = sum(p_strandness.values())
        n_single = sum(s_strandness.values())
        if n_paired and n_single:
            print(
                f"Warning: sample contains both paired ({n_paired}) and unpaired "
                f"({n_single}) reads; reporting the {'paired' if n_paired >= n_single else 'unpaired'} reads only.",
                file=sys.stderr,
            )
        if n_paired and n_paired >= n_single:
            total = float(n_paired)
            spec1 = (p_strandness['1++'] + p_strandness['1--'] + p_strandness['2+-'] + p_strandness['2-+']) / total
            spec2 = (p_strandness['1+-'] + p_strandness['1-+'] + p_strandness['2++'] + p_strandness['2--']) / total
            return ["PairEnd", spec1, spec2, 1 - spec1 - spec2]
        if n_single:
            total = float(n_single)
            spec1 = (s_strandness['++'] + s_strandness['--']) / total
            spec2 = (s_strandness['+-'] + s_strandness['-+']) / total
            return ["SingleEnd", spec1, spec2, 1 - spec1 - spec2]
        return ["Mixture", 0, 0, 0]

    def bamTowig(self,outfile,chrom_sizes, chrom_file,skip_multi=True,strand_rule=None,WigSumFactor=None,q_cut=30, window=10_000_000):
        """Write per-base read coverage as wiggle (variableStep) and BigWig.

        Each aligned base of each read adds 1 (soft clips, deletions and
        introns add nothing). With *strand_rule* (from infer_experiment, e.g.
        ``"1++,1--,2+-,2-+"``) coverage is split into ``.Forward.wig`` and
        ``.Reverse.wig`` (reverse values are negative); reads the rule does
        not cover (e.g. unpaired reads under a paired-end rule) are skipped
        and counted. With *WigSumFactor* every value is multiplied by it.
        Reads below *q_cut* are skipped when *skip_multi* is True, using the
        same rule as :meth:`calWigSum`.
        """
        import shutil
        import subprocess
        import numpy as np

        rule = readfilter.parse_strand_rule(strand_rule)
        min_mapq = q_cut if skip_multi else 0
        factor = 1.0 if WigSumFactor is None else float(WigSumFactor)
        if rule:
            paths = {'+': outfile + '.Forward.wig', '-': outfile + '.Reverse.wig'}
        else:
            paths = {None: outfile + '.wig'}
        handles = {key: open(path, 'w') for key, path in paths.items()}
        references = set(self.samfile.references)
        n_unassigned = 0

        try:
            for chr_name, chr_size in chrom_sizes.items():
                if chr_name not in references:
                    print("No alignments for " + chr_name + '. skipped', file=sys.stderr)
                    continue
                print("Processing " + chr_name + " ...", file=sys.stderr)
                for handle in handles.values():
                    handle.write('variableStep chrom=' + chr_name + '\n')

                for win_st in range(0, chr_size, window):
                    win_end = min(chr_size, win_st + window)
                    span = win_end - win_st
                    starts = {key: [] for key in handles}
                    ends = {key: [] for key in handles}
                    for aligned_read in self.samfile.fetch(chr_name, win_st, win_end):
                        if not readfilter.passes(aligned_read, min_mapq=min_mapq):
                            continue
                        if rule:
                            key = readfilter.gene_strand(aligned_read, rule)
                            if key is None:
                                if aligned_read.reference_start >= win_st:
                                    n_unassigned += 1
                                continue
                        else:
                            key = None
                        for block_st, block_end in aligned_read.get_blocks():
                            block_st = max(block_st, win_st)
                            block_end = min(block_end, win_end)
                            if block_st < block_end:
                                starts[key].append(block_st - win_st)
                                ends[key].append(block_end - win_st)
                    for key, handle in handles.items():
                        if not starts[key]:
                            continue
                        diff = np.bincount(starts[key], minlength=span + 1)[: span + 1]
                        diff -= np.bincount(ends[key], minlength=span + 1)[: span + 1]
                        coverage = np.cumsum(diff[:span])
                        positions = np.flatnonzero(coverage)
                        values = coverage[positions] * (-factor if key == '-' else factor)
                        np.savetxt(
                            handle,
                            np.column_stack((positions + win_st + 1, values)),
                            fmt=('%d', '%.2f'),
                            delimiter='\t',
                        )
        finally:
            for handle in handles.values():
                handle.close()

        if n_unassigned:
            print(f"Warning: {n_unassigned} reads not covered by strand rule {strand_rule!r} were skipped.", file=sys.stderr)

        executable = shutil.which("wigToBigWig")
        if executable is None:
            print("wigToBigWig was not found; only wiggle files were written.", file=sys.stderr)
            return
        for path in paths.values():
            bigwig = path[: -len('.wig')] + '.bw'
            print("Run wigToBigWig " + path + ' ' + chrom_file + ' ' + bigwig, file=sys.stderr)
            result = subprocess.run([executable, "-clip", path, chrom_file, bigwig])
            if result.returncode != 0:
                print('Failed to run "wigToBigWig" for ' + path, file=sys.stderr)

    def calWigSum(self,chrom_sizes, skip_multi=True, q_cut=30):
        """Total aligned bases, using the same read filter as :meth:`bamTowig`.

        (Before RSeQC 5.1.0+lbzhao multi-mapped reads were identified here by the NH
        tag but in :meth:`bamTowig` by MAPQ, so the normalization total and
        the signal could count different reads.)
        """
        print("Calcualte wigsum ... ", file=sys.stderr)
        min_mapq = q_cut if skip_multi else 0
        references = set(self.samfile.references)
        wigsum = 0.0
        for chr_name in chrom_sizes:
            if chr_name not in references:
                print("No alignments for " + chr_name + '. skipped', file=sys.stderr)
                continue
            print("Processing " + chr_name + " ...", file=sys.stderr)
            for aligned_read in self.samfile.fetch(chr_name):
                if not readfilter.passes(aligned_read, min_mapq=min_mapq):
                    continue
                wigsum += sum(end - st for st, end in aligned_read.get_blocks())
        return wigsum

    def bam2fq(self, prefix: str, paired: bool = True) -> None:
        """Convert BAM/SAM alignments into FASTQ files."""

        complement_table = str.maketrans("ACGTNX", "TGCANX")

        def fastq_values(aligned_read):
            """Return the FASTQ name, sequence, and quality string."""
            read_name = aligned_read.query_name
            read_seq = aligned_read.query_sequence
            read_qual = aligned_read.qual
            if read_seq is None:
                raise ValueError(
                    f"Read {read_name!r} does not contain sequence information"
                )

            read_seq = read_seq.upper()

            if read_qual is None:
                raise ValueError(
                    f"Read {read_name!r} does not contain quality information"
                )

            if aligned_read.is_reverse:
                read_seq = read_seq.translate(complement_table)[::-1]
                read_qual = read_qual[::-1]

            return read_name, read_seq, read_qual

        def write_fastq(handle, name: str, sequence: str, quality: str) -> None:
            """Write one FASTQ record."""
            handle.write(f"@{name}\n")
            handle.write(f"{sequence}\n")
            handle.write("+\n")
            handle.write(f"{quality}\n")

        print(
            "Convert BAM/SAM file into FASTQ format ...",
            file=sys.stderr,
        )

        if paired:
            read1_count = 0
            read2_count = 0

            read1_path = f"{prefix}.R1.fastq"
            read2_path = f"{prefix}.R2.fastq"

            with open(read1_path, "w", encoding="ascii") as out1, open(
                read2_path,
                "w",
                encoding="ascii",
            ) as out2:
                for aligned_read in self.samfile:
                    # Write each read once: secondary/supplementary records
                    # repeat (or hold only part of) the primary read.
                    if readfilter.is_non_primary(aligned_read):
                        continue
                    read_name, read_seq, read_qual = fastq_values(aligned_read)

                    if aligned_read.is_read1:
                        read1_count += 1
                        name = (
                            read_name
                            if read_name.endswith("/1")
                            else f"{read_name}/1"
                        )
                        write_fastq(out1, name, read_seq, read_qual)

                    elif aligned_read.is_read2:
                        read2_count += 1
                        name = (
                            read_name
                            if read_name.endswith("/2")
                            else f"{read_name}/2"
                        )
                        write_fastq(out2, name, read_seq, read_qual)

            print("Done", file=sys.stderr)
            print(f"read_1 count: {read1_count}", file=sys.stderr)
            print(f"read_2 count: {read2_count}", file=sys.stderr)

        else:
            read_count = 0
            output_path = f"{prefix}.fastq"

            with open(output_path, "w", encoding="ascii") as output:
                for aligned_read in self.samfile:
                    # Write each read once: secondary/supplementary records
                    # repeat (or hold only part of) the primary read.
                    if readfilter.is_non_primary(aligned_read):
                        continue
                    read_name, read_seq, read_qual = fastq_values(aligned_read)
                    write_fastq(output, read_name, read_seq, read_qual)
                    read_count += 1

            print("Done", file=sys.stderr)
            print(f"read count: {read_count}", file=sys.stderr)
                    
    def readsNVC(self,outfile=None,nx=True, q_cut = 30):
        '''for each read, calculate nucleotide frequency vs position'''
        if outfile is None:
            outfile1 = self.fileName + ".NVC.xls"
            outfile2 = self.fileName +".NVC_plot.r"
        else:
            outfile1 = outfile + ".NVC.xls"
            outfile2 = outfile +".NVC_plot.r"
        FO=open(outfile1,'w')
        RS=open(outfile2,'w')
        #PPcount=0
        
        transtab = str.maketrans("ACGTNX","TGCANX")
        base_freq=collections.defaultdict(int)
        a_count=[]
        c_count=[]
        g_count=[]
        t_count=[]
        n_count=[]
        x_count=[]
        if self.bam_format:print("Read BAM file ... ", end=' ', file=sys.stderr)
        else:print("Read SAM file ... ", end=' ', file=sys.stderr)

        max_len = 0     # longest read: positions run 0 .. max_len-1 (reads may be trimmed)
        try:
            while(1):
                aligned_read = next(self.samfile)
                if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                if aligned_read.mapq < q_cut: continue
                if not aligned_read.query_sequence: continue           #SEQ is '*'

                RNA_read = aligned_read.seq.upper()     
                max_len = max(max_len, len(RNA_read))
                if aligned_read.is_reverse:
                    RNA_read = RNA_read.translate(transtab)[::-1]
                for i,j in enumerate(RNA_read):
                    key = str(i) + j
                    base_freq[key] += 1
        except StopIteration:
            print("Done", file=sys.stderr)
        
        print("generating data matrix ...", file=sys.stderr)
        print("Position\tA\tC\tG\tT\tN\tX", file=FO)
        for i in range(max_len):
            print(str(i) + '\t', end=' ', file=FO)
            print(str(base_freq[str(i) + "A"]) + '\t', end=' ', file=FO)
            a_count.append(str(base_freq[str(i) + "A"]))
            print(str(base_freq[str(i) + "C"]) + '\t', end=' ', file=FO)
            c_count.append(str(base_freq[str(i) + "C"]))
            print(str(base_freq[str(i) + "G"]) + '\t', end=' ', file=FO)
            g_count.append(str(base_freq[str(i) + "G"]))
            print(str(base_freq[str(i) + "T"]) + '\t', end=' ', file=FO)
            t_count.append(str(base_freq[str(i) + "T"]))
            print(str(base_freq[str(i) + "N"]) + '\t', end=' ', file=FO)
            n_count.append(str(base_freq[str(i) + "N"]))
            print(str(base_freq[str(i) + "X"]) + '\t', file=FO)
            x_count.append(str(base_freq[str(i) + "X"]))
        FO.close()
        
        #generating R scripts
        print("generating R script  ...", file=sys.stderr)
        print("position=c(" + ','.join([str(i) for i in range(max_len)]) + ')', file=RS)
        print("A_count=c(" + ','.join(a_count) + ')', file=RS)
        print("C_count=c(" + ','.join(c_count) + ')', file=RS)
        print("G_count=c(" + ','.join(g_count) + ')', file=RS)
        print("T_count=c(" + ','.join(t_count) + ')', file=RS)
        print("N_count=c(" + ','.join(n_count) + ')', file=RS)
        print("X_count=c(" + ','.join(x_count) + ')', file=RS)
        
        if nx:
            print("total= A_count + C_count + G_count + T_count + N_count + X_count", file=RS)
            print("ym=max(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total) + 0.05", file=RS)
            print("yn=min(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total)", file=RS)
            
            print("pdf(%s)" % r_string(outfile +".NVC_plot.pdf"), file=RS)
            print('plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")', file=RS)
            print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
            print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
            print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
            print('lines(position,N_count/total,type="o",pch=20,col="black")', file=RS)     
            print('lines(position,X_count/total,type="o",pch=20,col="grey")', file=RS)  
            print('legend('+ str(max_len-10) + ',ym,legend=c("A","T","G","C","N","X"),col=c("dark green","red","blue","cyan","black","grey"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan","black","grey"))', file=RS)
            print("dev.off()", file=RS)
        else:
            print("total= A_count + C_count + G_count + T_count", file=RS)
            print("ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05", file=RS)
            print("yn=min(A_count/total,C_count/total,G_count/total,T_count/total)", file=RS)
        
            print("pdf(%s)" % r_string(outfile +".NVC_plot.pdf"), file=RS)
            print('plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")', file=RS)
            print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
            print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
            print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
            print('legend('+ str(max_len-10) + ',ym,legend=c("A","T","G","C"),col=c("dark green","red","blue","cyan"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan"))', file=RS)
            print("dev.off()", file=RS)
        
        RS.close()

        return {
            "counts": {base: [base_freq[str(i) + base] for i in range(max_len)] for base in "ACGTNX"},
            "include_nx": nx,
        }

    def readsQual_boxplot(self,outfile,shrink=1000, q_cut=30):
        '''calculate phred quality score for each base in read (5->3)'''

        output = outfile + ".qual.r"
        FO=open(output,'w')

        if self.bam_format:print("Read BAM file ... ", end=' ', file=sys.stderr)
        else:print("Read SAM file ... ", end=' ', file=sys.stderr)

        quality = collections.defaultdict(dict) #read_pos=>quality score=>count
        q_max = -1
        q_min = 10000
        q_list=[]
        i_box={}    #key is read postion,value is 
        read_len = 0
        try:
            while(1):
                aligned_read = next(self.samfile)
                if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                if aligned_read.mapq < q_cut: continue
                if not aligned_read.query_qualities: continue          #QUAL is '*'

                qual_str = aligned_read.qqual
                read_len = max(read_len, len(qual_str))    # longest read (reads may be trimmed)
                if aligned_read.is_reverse:
                    qual_str = qual_str[::-1]

                for i,j in enumerate(qual_str):
                    q=ord(j)-33
                    if q > q_max: q_max = q
                    if q < q_min: q_min = q
                    quality[i][q] = quality[i].get(q, 0) + 1
        except StopIteration:
            print("Done", file=sys.stderr)
        
        for p in range(0,read_len):
            #print str(p) + ':',
            val=[]
            occurrence=[]
            for q in range(q_min,q_max+1):
                if p in quality and q in quality[p]:
                    val.append(str(q))              
                    occurrence.append(str(quality[p][q]))   
                    q_list.append(str(quality[p][q]))
                else:
                    q_list.append(str(0))
            i_box[p] = 'rep(c(' + ','.join(val) + '),times=c(' + ','.join(occurrence) + ')/' + str(shrink)+ ')'
        
        
        #generate R script for boxplot
        print("pdf(%s)" % r_string(outfile + ".qual.boxplot.pdf"), file=FO)
        for i in sorted(i_box):
            print('p'+str(i) + '<-' + i_box[i], file=FO)
        print('boxplot(' + ','.join(['p'+str(i) for i in i_box]) + ',xlab=\"Position of Read(5\'->3\')\",ylab=\"Phred Quality Score\",outline=F' + ')', file=FO)
        print("dev.off()", file=FO)
        
        
        #generate R script for heatmap
        print('\n', file=FO)
        print("pdf(%s)" % r_string(outfile + ".qual.heatmap.pdf"), file=FO)
        print("qual=c(" + ','.join(q_list)  + ')', file=FO)
        print("mat=matrix(qual,ncol=%s,byrow=F)" % (read_len), file=FO)
        print('Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c(\'spline\'))', file=FO)
        print("heatmap(mat,Rowv=NA,Colv=NA,xlab=\"Position of Read\",ylab=\"Phred Quality Score\",labRow=seq(from=%s,to=%s),col = Lab.palette(256),scale=\"none\" )" % (q_min,q_max), file=FO)
        print('dev.off()', file=FO)
        FO.close()
        return {"quality": {pos: dict(quality[pos]) for pos in range(read_len)}}

    def readGC(self,outfile=None, q_cut=30):
        '''GC content distribution of reads'''
        if outfile is None:
            outfile1 = self.fileName + ".GC.xls"
            outfile2 = self.fileName +".GC_plot.r"
        else:
            outfile1 = outfile + ".GC.xls"
            outfile2 = outfile + ".GC_plot.r"
        FO=open(outfile1,'w')
        RS=open(outfile2,'w')
        
        gc_hist=collections.defaultdict(int)    #key is GC percent, value is count of reads

        if self.bam_format:print("Read BAM file ... ", end=' ', file=sys.stderr)
        else:print("Read SAM file ... ", end=' ', file=sys.stderr)

        try:
            while(1):
                aligned_read = next(self.samfile)
                if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                if aligned_read.is_unmapped:continue    #skip unmapped read
                if aligned_read.is_qcfail:continue      #skip low quality
                if aligned_read.mapq < q_cut: continue
                RNA_read = aligned_read.seq.upper()     
                gc_percent = "%4.2f" % ((RNA_read.count('C') + RNA_read.count('G'))/(len(RNA_read)+0.0)*100)
                #print gc_percent
                gc_hist[gc_percent] += 1
        except StopIteration:
            print("Done", file=sys.stderr)
        
        print("writing GC content ...", file=sys.stderr)    
        print("GC%\tread_count", file=FO)
        for i in list(gc_hist.keys()):
            print(i + '\t' + str(gc_hist[i]), file=FO)
            
        print("writing R script ...", file=sys.stderr)
        print("pdf(%s)" % r_string(outfile +  ".GC_plot.pdf"), file=RS)
        print('gc=rep(c(' + ','.join([i for i in list(gc_hist.keys())]) + '),' + 'times=c(' + ','.join([str(i) for i in list(gc_hist.values())]) + '))', file=RS)
        print('hist(gc,probability=T,breaks=%d,xlab="GC content (%%)",ylab="Density of Reads",border="blue",main="")' % 100, file=RS)
        #print >>RS, "lines(density(gc),col='red')"
        print("dev.off()", file=RS)     
        FO.close()
        RS.close()
        return {"gc_hist": {float(k): v for k, v in gc_hist.items()}}

    def readDupRate(self,q_cut, outfile=None,up_bound=500):
        '''Calculate read duplication rates by sequence and by mapping position.

        Duplicate-flagged reads are counted (that is what is measured);
        secondary and supplementary records are not. For paired reads the
        mapping key includes the mate position. In RNA-seq, highly expressed
        genes produce many reads at the same position without PCR, so a high
        rate does not by itself indicate PCR duplication.
        '''
        if outfile is None:
            outfile1 = self.fileName + ".seq.DupRate.xls"
            outfile2 = self.fileName + ".pos.DupRate.xls"
            outfile3 = self.fileName + ".DupRate_plot.r"
        else:
            outfile1 = outfile + ".seq.DupRate.xls"
            outfile2 = outfile + ".pos.DupRate.xls"
            outfile3 = outfile +".DupRate_plot.r"
        SEQ=open(outfile1,'w')
        POS=open(outfile2,'w')
        RS=open(outfile3,'w')
        
        seqDup=collections.defaultdict(int)
        posDup=collections.defaultdict(int)
        
        seqDup_count=collections.defaultdict(int)
        posDup_count=collections.defaultdict(int)

        if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
        else:print("Load SAM file ... ", end=' ', file=sys.stderr)

        try:
            while(1):
                exon_boundary=""
                aligned_read = next(self.samfile)
                if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                if aligned_read.is_unmapped:continue    #skip unmapped read
                if aligned_read.is_qcfail:continue      #skip low quality
                if aligned_read.mapq < q_cut: continue
                RNA_read = aligned_read.seq.upper()     
                seqDup[RNA_read] +=1                    #key is read sequence

                chrom = self.samfile.getrname(aligned_read.tid)
                hit_st = aligned_read.pos
                exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)               
                for ex in exon_blocks:
                    exon_boundary += str(ex[1]) + '-' + str(ex[2]) + ":"
                key = chrom + ":" + str(hit_st) + ":" + exon_boundary
                if aligned_read.is_paired and not aligned_read.mate_is_unmapped:
                    # A pair is a positional duplicate only if both ends match.
                    key += "mate=" + str(aligned_read.next_reference_id) + ":" + str(aligned_read.next_reference_start)
                posDup[key] +=1

        except StopIteration:
            print("Done", file=sys.stderr)

        print("report duplicte rate based on sequence ...", file=sys.stderr)
        print("Occurrence\tUniqReadNumber", file=SEQ)
        for i in list(seqDup.values()):         #key is occurence, value is uniq reads number (based on seq)
            seqDup_count[i] +=1
        for k in sorted(seqDup_count.keys()):   
            print(str(k) +'\t'+ str(seqDup_count[k]), file=SEQ)
        SEQ.close()
        
        print("report duplicte rate based on mapping  ...", file=sys.stderr)
        print("Occurrence\tUniqReadNumber", file=POS)
        for i in list(posDup.values()):         #key is occurence, value is uniq reads number (based on coord)
            posDup_count[i] +=1
        for k in sorted(posDup_count.keys()):   
            print(str(k) +'\t'+ str(posDup_count[k]), file=POS)
        POS.close()
        
        
        print("generate R script ...", file=sys.stderr)
        print("pdf(%s)" % r_string(outfile +".DupRate_plot.pdf"), file=RS)
        print("par(mar=c(5,4,4,5),las=0)", file=RS)
        print("seq_occ=c(" + ','.join([str(i) for i in sorted(seqDup_count.keys()) ]) + ')', file=RS)
        print("seq_uniqRead=c(" + ','.join([str(seqDup_count[i]) for i in sorted(seqDup_count.keys()) ]) + ')', file=RS)
        print("pos_occ=c(" + ','.join([str(i) for i in sorted(posDup_count.keys()) ]) + ')', file=RS)
        print("pos_uniqRead=c(" + ','.join([str(posDup_count[i]) for i in sorted(posDup_count.keys()) ]) + ')', file=RS)
        print("plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Occurrence of read',pch=4,cex=0.8,col='blue',xlim=c(1,%d),yaxt='n')" % up_bound, file=RS)
        print("points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')", file=RS)
        print('ym=floor(max(log10(pos_uniqRead)))', file=RS)
        print("legend(%d,ym,legend=c('Sequence-based','Mapping-based'),col=c('blue','red'),pch=c(4,20))" % max(up_bound-200,1), file=RS)
        print('axis(side=2,at=0:ym,labels=0:ym)', file=RS)
        print('axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[2]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[3]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[4]*100/sum(pos_uniqRead*pos_occ))))', file=RS)
        print('mtext(4, text = "Reads %", line = 2)', file=RS)
        print('dev.off()', file=RS)
        RS.close()
        return {"sequence": dict(seqDup_count), "position": dict(posDup_count), "up_bound": up_bound}

    def clipping_profile(self, outfile, q_cut, PE, type="S"):
        """Calculate the positional profile of soft clipping.

        For variable-length reads, the denominator at each position is the
        number of reads that actually contain that position. For paired-end
        data, read 1 and read 2 lengths are tracked independently.
        """
        out_file1 = outfile + ".clipping_profile.xls"
        out_file2 = outfile + ".clipping_profile.r"
        OUT = open(out_file1, "w")
        ROUT = open(out_file2, "w")

        print("Position\tClipped_nt\tNon_clipped_nt", file=OUT)

        if self.bam_format:
            print("Load BAM file ... ", end=" ", file=sys.stderr)
        else:
            print("Load SAM file ... ", end=" ", file=sys.stderr)

        # single-end sequencing
        if PE is False:
            total_read = 0
            max_len = 0
            soft_clip_profile = collections.defaultdict(int)
            position_coverage = collections.defaultdict(int)

            try:
                while True:
                    aligned_read = next(self.samfile)
                    if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                    if aligned_read.mapq < q_cut:
                        continue
                    if aligned_read.is_unmapped:
                        continue
                    if aligned_read.is_qcfail:
                        continue

                    total_read += 1
                    cigar_str = bam_cigar.list2longstr(aligned_read.cigar)
                    read_len = len(cigar_str)
                    max_len = max(max_len, read_len)

                    # Each accepted read contributes to the denominator only
                    # at positions that exist in that read.
                    for indx in range(read_len):
                        position_coverage[indx] += 1

                    if aligned_read.is_reverse:
                        cigar_str = cigar_str[::-1]

                    if type not in cigar_str:
                        continue

                    for indx, symbol in enumerate(cigar_str):
                        if symbol == type:
                            soft_clip_profile[indx] += 1
            except StopIteration:
                print("Done", file=sys.stderr)

            print("Total reads used: %d" % total_read, file=sys.stderr)
            read_pos = list(range(max_len))
            clip_count = []
            nonclip_count = []

            for i in read_pos:
                clipped = soft_clip_profile[i]
                nonclipped = position_coverage[i] - clipped
                print(
                    str(i) + "\t" + str(clipped) + "\t" + str(nonclipped),
                    file=OUT,
                )
                clip_count.append(clipped)
                nonclip_count.append(nonclipped)

            print(
                "pdf(%s)" % r_string(outfile + ".clipping_profile.pdf"),
                file=ROUT,
            )
            print(
                "read_pos=c(%s)"
                % ",".join(str(i) for i in read_pos),
                file=ROUT,
            )
            print(
                "clip_count=c(%s)"
                % ",".join(str(i) for i in clip_count),
                file=ROUT,
            )
            print(
                "nonclip_count=c(%s)"
                % ",".join(str(i) for i in nonclip_count),
                file=ROUT,
            )
            print(
                'plot(read_pos, nonclip_count*100/(clip_count+nonclip_count),'
                'col="blue",main="clipping profile",xlab="Position of read",'
                'ylab="Non-clipped %",type="b")',
                file=ROUT,
            )
            print("dev.off()", file=ROUT)

        # paired-end sequencing
        if PE is True:
            total_read1 = 0
            total_read2 = 0
            max_len_r1 = 0
            max_len_r2 = 0
            r1_soft_clip_profile = collections.defaultdict(int)
            r2_soft_clip_profile = collections.defaultdict(int)
            r1_position_coverage = collections.defaultdict(int)
            r2_position_coverage = collections.defaultdict(int)

            try:
                while True:
                    aligned_read = next(self.samfile)
                    if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                    if aligned_read.mapq < q_cut:
                        continue
                    if aligned_read.is_unmapped:
                        continue
                    if aligned_read.is_qcfail:
                        continue
                    if not aligned_read.is_paired:
                        continue

                    cigar_str = bam_cigar.list2longstr(aligned_read.cigar)
                    read_len = len(cigar_str)

                    # Orient read positions consistently from the 5' end
                    # before counting clipped positions.
                    if aligned_read.is_reverse:
                        cigar_str = cigar_str[::-1]

                    if aligned_read.is_read1:
                        total_read1 += 1
                        max_len_r1 = max(max_len_r1, read_len)
                        for indx in range(read_len):
                            r1_position_coverage[indx] += 1
                        if type in cigar_str:
                            for indx, symbol in enumerate(cigar_str):
                                if symbol == type:
                                    r1_soft_clip_profile[indx] += 1

                    elif aligned_read.is_read2:
                        total_read2 += 1
                        max_len_r2 = max(max_len_r2, read_len)
                        for indx in range(read_len):
                            r2_position_coverage[indx] += 1
                        if type in cigar_str:
                            for indx, symbol in enumerate(cigar_str):
                                if symbol == type:
                                    r2_soft_clip_profile[indx] += 1

            except StopIteration:
                print("Done", file=sys.stderr)

            read_pos_r1 = list(range(max_len_r1))
            read_pos_r2 = list(range(max_len_r2))
            r1_clip_count = []
            r1_nonclip_count = []
            r2_clip_count = []
            r2_nonclip_count = []

            print("Total read-1 used: %d" % total_read1, file=sys.stderr)
            print("Total read-2 used: %d" % total_read2, file=sys.stderr)

            print("Read-1:", file=OUT)
            for i in read_pos_r1:
                clipped = r1_soft_clip_profile[i]
                nonclipped = r1_position_coverage[i] - clipped
                print(
                    str(i) + "\t" + str(clipped) + "\t" + str(nonclipped),
                    file=OUT,
                )
                r1_clip_count.append(clipped)
                r1_nonclip_count.append(nonclipped)

            print("Read-2:", file=OUT)
            for i in read_pos_r2:
                clipped = r2_soft_clip_profile[i]
                nonclipped = r2_position_coverage[i] - clipped
                print(
                    str(i) + "\t" + str(clipped) + "\t" + str(nonclipped),
                    file=OUT,
                )
                r2_clip_count.append(clipped)
                r2_nonclip_count.append(nonclipped)

            print(
                "pdf(%s)" % r_string(outfile + ".clipping_profile.R1.pdf"),
                file=ROUT,
            )
            print(
                "read_pos_r1=c(%s)"
                % ",".join(str(i) for i in read_pos_r1),
                file=ROUT,
            )
            print(
                "r1_clip_count=c(%s)"
                % ",".join(str(i) for i in r1_clip_count),
                file=ROUT,
            )
            print(
                "r1_nonclip_count=c(%s)"
                % ",".join(str(i) for i in r1_nonclip_count),
                file=ROUT,
            )
            print(
                'plot(read_pos_r1, '
                'r1_nonclip_count*100/(r1_clip_count+r1_nonclip_count),'
                'col="blue",main="clipping profile",'
                'xlab="Position of read (read-1)",'
                'ylab="Non-clipped %",type="b")',
                file=ROUT,
            )
            print("dev.off()", file=ROUT)

            print(
                "pdf(%s)" % r_string(outfile + ".clipping_profile.R2.pdf"),
                file=ROUT,
            )
            print(
                "read_pos_r2=c(%s)"
                % ",".join(str(i) for i in read_pos_r2),
                file=ROUT,
            )
            print(
                "r2_clip_count=c(%s)"
                % ",".join(str(i) for i in r2_clip_count),
                file=ROUT,
            )
            print(
                "r2_nonclip_count=c(%s)"
                % ",".join(str(i) for i in r2_nonclip_count),
                file=ROUT,
            )
            print(
                'plot(read_pos_r2, '
                'r2_nonclip_count*100/(r2_clip_count+r2_nonclip_count),'
                'col="blue",main="clipping profile",'
                'xlab="Position of read (read-2)",'
                'ylab="Non-clipped %",type="b")',
                file=ROUT,
            )
            print("dev.off()", file=ROUT)

        OUT.close()
        ROUT.close()
        if PE:
            return {
                "R1": (read_pos_r1, r1_clip_count, r1_nonclip_count),
                "R2": (read_pos_r2, r2_clip_count, r2_nonclip_count),
            }
        return {"": (read_pos, clip_count, nonclip_count)}

    def insertion_profile(self,outfile, q_cut, PE, type="I"):
        '''calculate profile of insertion'''
        
        out_file1 = outfile + ".insertion_profile.xls"
        out_file2 = outfile + ".insertion_profile.r"
        OUT = open(out_file1,'w')
        ROUT = open(out_file2,'w')
        
        print("Position\tInsert_nt\tNon_insert_nt", file=OUT)
        
        if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
        else:print("Load SAM file ... ", end=' ', file=sys.stderr)
        
        cigar_str = ""
        max_len = 0
        #single end sequencing
        if PE is False:
            total_read = 0.0        
            soft_clip_profile = collections.defaultdict(int)
            try:
                while(1):
                    aligned_read = next(self.samfile)
                    if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                    if aligned_read.mapq < q_cut: continue
                    if aligned_read.is_unmapped:continue    #skip unmapped read
                    if aligned_read.is_qcfail:continue      #skip low quality
                
                    total_read +=1
                    cigar_str = bam_cigar.list2longstr(aligned_read.cigar)  # ([(0, 9), (4, 1)] ==> MMMMMMMMMS
                    max_len = max(max_len, len(cigar_str))
                    if type not in cigar_str:   # no insertion
                        continue
                    if aligned_read.is_reverse:
                        cigar_str = cigar_str[::-1]     

                    for indx,symbl in enumerate(cigar_str):
                        if symbl == type:
                            soft_clip_profile[indx] += 1.0          
            except StopIteration:
                print("Done", file=sys.stderr)
            
            print("Totoal reads used: %d" % int(total_read), file=sys.stderr)
            read_pos = list(range(0, max_len))
            clip_count = []
            for i in read_pos:
                print(str(i) + '\t' + str(soft_clip_profile[i]) + '\t' + str(total_read - soft_clip_profile[i]), file=OUT)
                clip_count.append(soft_clip_profile[i])
            
            
            print("pdf(%s)" % r_string(outfile + '.insertion_profile.pdf'), file=ROUT)
            print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
            print("insert_count=c(%s)" % ','.join([str(i) for i in clip_count]), file=ROUT)
            print("noninsert_count= %d - insert_count" % (total_read), file=ROUT)
            print('plot(read_pos, insert_count*100/(insert_count+noninsert_count),col="blue",main="Insertion profile",xlab="Position of read",ylab="Insertion %",type="b")', file=ROUT)
            print("dev.off()", file=ROUT)

        
        if PE is True:
            total_read1 = 0.0
            total_read2 = 0.0
            r1_soft_clip_profile = collections.defaultdict(int)
            r2_soft_clip_profile = collections.defaultdict(int)
            try:
                while(1):
                    aligned_read = next(self.samfile)
                    if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                    if aligned_read.mapq < q_cut: continue
                    if aligned_read.is_unmapped:continue    #skip unmapped read
                    if aligned_read.is_qcfail:continue      #skip low quality
                    if not aligned_read.is_paired: continue
                    if aligned_read.is_read1:
                        total_read1 += 1
                    if aligned_read.is_read2:
                        total_read2 += 1                                            
                    cigar_str = bam_cigar.list2longstr(aligned_read.cigar)  # ([(0, 9), (4, 1)] ==> MMMMMMMMMS
                    max_len = max(max_len, len(cigar_str))
                    if aligned_read.is_reverse:
                        cigar_str = cigar_str[::-1]                 
                    
                    if type not in cigar_str:   # no clipping
                        continue
                    
                    if aligned_read.is_read1:
                        for indx,symbl in enumerate(cigar_str):
                            if symbl == type:
                                r1_soft_clip_profile[indx] += 1.0
                    if aligned_read.is_read2:
                        for indx,symbl in  enumerate(cigar_str):
                            if symbl == type:
                                r2_soft_clip_profile[indx] += 1.0
            except StopIteration:
                print("Done", file=sys.stderr)
            
            read_pos = list(range(0, max_len))
            r1_clip_count = []
            r2_clip_count = []
            
            print("Totoal read-1 used: %d" % int(total_read1), file=sys.stderr)
            print("Totoal read-2 used: %d" % int(total_read2), file=sys.stderr)
            print("Read-1:", file=OUT)
            for i in read_pos:
                print(str(i) + '\t' + str(r1_soft_clip_profile[i]) + '\t' + str(total_read1 - r1_soft_clip_profile[i]), file=OUT)       
                r1_clip_count.append(r1_soft_clip_profile[i])
                

            print("Read-2:", file=OUT)
            for i in read_pos:
                print(str(i) + '\t' + str(r2_soft_clip_profile[i]) + '\t' + str(total_read2 - r2_soft_clip_profile[i]), file=OUT)
                r2_clip_count.append(r2_soft_clip_profile[i])
        
            
            print("pdf(%s)" % r_string(outfile + '.insertion_profile.R1.pdf'), file=ROUT)
            print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
            print("r1_insert_count=c(%s)" % ','.join([str(i) for i in r1_clip_count]), file=ROUT)
            print("r1_noninsert_count = %d - r1_insert_count" % (total_read1), file=ROUT)
            print('plot(read_pos, r1_insert_count*100/(r1_insert_count + r1_noninsert_count),col="blue",main="Insertion profile",xlab="Position of read (read-1)",ylab="Insertion %",type="b")', file=ROUT)
            print("dev.off()", file=ROUT)

            print("pdf(%s)" % r_string(outfile + '.insertion_profile.R2.pdf'), file=ROUT)
            print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
            print("r2_insert_count=c(%s)" % ','.join([str(i) for i in r2_clip_count]), file=ROUT)
            print("r2_noninsert_count = %d - r2_insert_count" % (total_read2), file=ROUT)
            print('plot(read_pos, r2_insert_count*100/(r2_insert_count + r2_noninsert_count),col="blue",main="Insertion profile",xlab="Position of read (read-2)",ylab="Insertion %",type="b")', file=ROUT)
            print("dev.off()", file=ROUT)
        OUT.close()
        ROUT.close()
        if PE:
            return {
                "R1": (read_pos, r1_clip_count, [total_read1 - c for c in r1_clip_count]),
                "R2": (read_pos, r2_clip_count, [total_read2 - c for c in r2_clip_count]),
            }
        return {"": (read_pos, clip_count, [total_read - c for c in clip_count])}

    def mRNA_inner_distance(self,outfile,refbed,low_bound=0,up_bound=1000,step=10,sample_size=1000000, q_cut=30, seed=0):
        """Estimate the inner distance (insert size minus both read lengths)
        of paired-end fragments.

        For each sampled pair the genomic gap between the end of the leftmost
        read and the start of its mate is measured. If both reads lie in the
        same transcript, only exonic bases of that transcript's exons (all its
        isoforms) are counted, giving the distance on the mRNA. Overlapping
        mates give a negative distance. Pairs are sampled across the genome
        (see :mod:`qcmodule.sampling`).

        The mean, median and SD printed by the R script use only pairs inside
        ``[low_bound, up_bound)``; the numbers of pairs outside that window
        are reported on stderr.
        """
        out_file1 = outfile + ".inner_distance.txt"
        out_file2 = outfile + ".inner_distance_freq.txt"
        out_file3 = outfile + ".inner_distance_plot.r"

        window_left_bound = list(range(low_bound, up_bound, step))
        distances = []
        sizes = []
        counts = []

        print("Get exon regions from " + refbed + " ...", file=sys.stderr)
        transcripts = load_transcript_exons(refbed)
        transcript_ranges = {}
        regions = []
        for index, (chrom, tx_start, tx_end, _exons) in enumerate(transcripts):
            transcript_ranges.setdefault(chrom.upper(), Intersecter()).add_interval(
                Interval(tx_start, tx_end, value=index)
            )
            regions.append((chrom, tx_start, tx_end))

        def transcripts_at(chrom, start, end):
            tree = transcript_ranges.get(chrom)
            return {iv.value for iv in tree.find(start, end)} if tree is not None else set()

        def usable(read):
            if not readfilter.passes(read, min_mapq=q_cut):
                return False
            if not read.is_paired or read.mate_is_unmapped:
                return False
            # Each pair is measured once, from its leftmost read.
            if read.next_reference_id == read.reference_id:
                if read.next_reference_start < read.reference_start:
                    return False
                if read.next_reference_start == read.reference_start and read.is_read1:
                    return False
            return True

        if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
        else:print("Load SAM file ... ", end=' ', file=sys.stderr)

        pair_num = 0
        with open(out_file1, 'w') as FO:
            for aligned_read in sampling.sample_alignments(self.samfile, regions, sample_size, usable, seed=seed):
                pair_num += 1
                if aligned_read.next_reference_id != aligned_read.reference_id:
                    FO.write(aligned_read.query_name + '\tNA\tsameChrom=No\n')
                    continue

                chrom = aligned_read.reference_name.upper()
                read1_end = aligned_read.reference_end
                read2_start = aligned_read.next_reference_start

                if read2_start >= read1_end:
                    inner_distance = read2_start - read1_end
                else:
                    # Mates overlap: minus the read-1 aligned bases at or after the mate start.
                    inner_distance = -sum(
                        max(0, end - max(start, read2_start))
                        for start, end in aligned_read.get_blocks()
                    )

                shared = transcripts_at(chrom, read1_end - 1, read1_end) & transcripts_at(chrom, read2_start, read2_start + 1)
                if not shared:
                    label = 'sameTranscript=No,dist=genomic'
                    value = inner_distance
                elif inner_distance <= 0:
                    label = 'readPairOverlap'
                    value = inner_distance
                else:
                    exons = BED.unionBed3([[chrom, st, end] for i in shared for st, end in transcripts[i][3]])
                    size = sum(
                        max(0, min(end, read2_start) - max(st, read1_end))
                        for _, st, end in exons
                    )
                    if size == inner_distance:
                        label = 'sameTranscript=Yes,sameExon=Yes,dist=mRNA'
                        value = size
                    elif size > 0:
                        label = 'sameTranscript=Yes,sameExon=No,dist=mRNA'
                        value = size
                    else:
                        label = 'sameTranscript=Yes,nonExonic=Yes,dist=genomic'
                        value = inner_distance
                FO.write(aligned_read.query_name + '\t' + str(value) + '\t' + label + '\n')
                distances.append(value)
        print("Done", file=sys.stderr)

        print("Total read pairs  used " + str(pair_num), file=sys.stderr)
        if pair_num==0:
            raise NoDataFound("Cannot find paired reads")
        below = sum(1 for d in distances if d <= low_bound)
        above = sum(1 for d in distances if d > up_bound)
        print(f"Pairs with inner distance <= {low_bound}: {below}; > {up_bound}: {above}. "
              "These are not included in the histogram or in the mean/median/SD.", file=sys.stderr)

        # Bin (st, st + step]: the same bins as earlier releases.
        FQ = open(out_file2, 'w')
        RS = open(out_file3, 'w')
        for st in window_left_bound:
            sizes.append(str(st + step/2))
            count = str(sum(1 for d in distances if st < d <= st + step))
            counts.append(count)
            print(str(st) + '\t' + str(st+step) +'\t' + count, file=FQ)
        print("out_file = %s" % r_string(outfile), file=RS)
        print("pdf(%s)" % r_string(outfile + ".inner_distance_plot.pdf"), file=RS)
        #print >>RS, "par(mfrow=c(2,1),cex.main=0.8,cex.lab=0.8,cex.axis=0.8,mar=c(4,4,4,1))"
        #print >>RS, 'pie(c(%d,%d,%d),col=rainbow(3),cex=0.5,radius=1,main="Total %d fragments",labels=c("fraSize <= %d\\n(%4.2f%%)","fragSize > %d\\n(%4.2f%%)","%d < fragSize <= %d\\n(%4.2f%%)"), density=rep(80,80,80),angle=c(90,140,170))' % (ultra_low, ultra_high, pair_num -ultra_low -ultra_high, pair_num, low_bound, ultra_low*100/pair_num, up_bound, ultra_high*100/pair_num, low_bound, up_bound, 100-ultra_low*100/pair_num - ultra_high*100/pair_num)
        print('fragsize=rep(c(' + ','.join(sizes) + '),' + 'times=c(' + ','.join(counts) + '))', file=RS)
        print('frag_sd = sd(fragsize)', file=RS)
        print('frag_mean = mean(fragsize)', file=RS)
        print('frag_median = median(fragsize)', file=RS)
        print('write(x=c("Name","Mean","Median","sd"), sep="\t", file=stdout(),ncolumns=4)', file=RS)
        print('write(c(out_file,frag_mean,frag_median,frag_sd),sep="\t", file=stdout(),ncolumns=4)', file=RS)
        print('hist(fragsize,probability=T,breaks=%d,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")' % len(window_left_bound), file=RS)
        print("lines(density(fragsize,bw=%d),col='red')" % (2*step), file=RS)
        print("dev.off()", file=RS)
        FQ.close()
        RS.close()
        return {"distances": distances, "low_bound": low_bound, "up_bound": up_bound, "step": step}

    def annotate_junction(self,refgene,outfile,min_intron=50, q_cut=30):
        """Annotate splice junctions in a BAM/SAM file against a BED12 gene model.

        A read can be spliced more than once (several splicing *events*); all
        events spanning the same intron form one *junction*. Each junction is
        classified by :func:`classify_junction`:

        * ``annotated``: the exact intron (start and end) is in the gene model;
        * ``novel_combination``: both splice sites are annotated but never
          joined in the gene model (e.g. exon skipping);
        * ``partial_novel``: one splice site is annotated;
        * ``complete_novel``: neither splice site is annotated.

        Before RSeQC 5.1.0+lbzhao ``novel_combination`` junctions were reported as
        ``annotated``. Chromosome names are matched case-insensitively and
        written as they appear in the BAM header.
        """
        out_file = outfile + ".junction.xls"
        out_file2 = outfile + ".junction_plot.r"
        if refgene is None:
            raise ValueError("You must provide reference gene model in bed format.")

        print("Reading reference bed file: ", refgene, " ... ", end=' ', file=sys.stderr)
        reference = load_reference_introns(refgene)
        print("Done", file=sys.stderr)

        categories = JUNCTION_CATEGORIES
        event_counts = dict.fromkeys(categories, 0)
        total_events = 0
        filtered_events = 0
        splicing_events = collections.defaultdict(int)   # (chrom, st, end) -> reads

        if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
        else:print("Load SAM file ... ", end=' ', file=sys.stderr)

        for aligned_read in self.samfile:
            if not readfilter.passes(aligned_read, min_mapq=q_cut):
                continue
            chrom = self.samfile.get_reference_name(aligned_read.reference_id)
            for _, i_st, i_end in bam_cigar.fetch_intron(chrom, aligned_read.reference_start, aligned_read.cigartuples):
                total_events += 1
                if i_end - i_st < min_intron:
                    filtered_events += 1
                    continue
                splicing_events[(chrom, i_st, i_end)] += 1
                event_counts[classify_junction(reference, chrom, i_st, i_end)] += 1
        print("Done", file=sys.stderr)

        print("total = " + str(total_events))
        if total_events == 0:
            raise NoDataFound("No splice junction found.")

        junction_counts = dict.fromkeys(categories, 0)
        with open(out_file, 'w') as OUT:
            print("chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation", file=OUT)
            for (chrom, i_st, i_end), n_reads in splicing_events.items():
                label = classify_junction(reference, chrom, i_st, i_end)
                junction_counts[label] += 1
                print(f"{chrom}\t{i_st}\t{i_end}\t{n_reads}\t{label}", file=OUT)
        total_junctions = sum(junction_counts.values())

        names = {
            "annotated": "Known",
            "novel_combination": "Novel Combination",
            "partial_novel": "Partial Novel",
            "complete_novel": "Novel",
        }
        print("\n===================================================================", file=sys.stderr)
        print("Total splicing  Events:\t" + str(total_events), file=sys.stderr)
        for label in categories:
            print(f"{names[label]} Splicing Events:\t{event_counts[label]}", file=sys.stderr)
        print("Filtered Splicing Events:\t" + str(filtered_events), file=sys.stderr)

        if total_junctions == 0:
            raise ValueError("No splice junction passed the minimum intron size")
        print("\nTotal splicing  Junctions:\t" + str(total_junctions), file=sys.stderr)
        for label in categories:
            print(f"{names[label]} Splicing Junctions:\t{junction_counts[label]}", file=sys.stderr)
        print("\n===================================================================", file=sys.stderr)

        kept_events = sum(event_counts.values())
        with open(out_file2, 'w') as ROUT:
            for title, counts, total, pdf in (
                ("splicing events", event_counts, kept_events, ".splice_events.pdf"),
                ("splicing junctions", junction_counts, total_junctions, ".splice_junction.pdf"),
            ):
                values = [counts[label] * 100.0 / total if total else 0 for label in categories]
                labels = ",".join(
                    '"%s %d%%"' % (label, round(v)) for label, v in zip(categories, values)
                )
                print("pdf(%s)" % r_string(outfile + pdf), file=ROUT)
                print("x=c(" + ','.join(str(v) for v in values) + ')', file=ROUT)
                print('pie(x,col=c(4,6,2,3),init.angle=30,main="%s",labels=c(%s))' % (title, labels), file=ROUT)
                print("dev.off()", file=ROUT)
        return {"events": event_counts, "junctions": junction_counts}

    def saturation_junction(self,refgene,outfile=None,sample_start=5,sample_step=5,sample_end=100,min_intron=50,recur=1, q_cut=30, seed=None):
        '''check if an RNA-seq experiment is saturated in terms of detecting known splicing junction.

        Splice events are shuffled (with *seed*, for reproducible results) and
        added cumulatively at each sampling percentage. A junction is counted
        (all, known or novel) once it is supported by at least *recur* reads.
        '''
        
        out_file = outfile + ".junctionSaturation_plot.r"
        if refgene is None:
            raise ValueError("You must provide reference gene model in bed format.")
        
        OUT = open(out_file,'w')


        #reading reference gene 
        knownSpliceSites= set()
        chrom_list=set()
        print("reading reference bed file: ",refgene, " ... ", end=' ', file=sys.stderr)
        for line in open(refgene,'r'):
            if line.startswith(('#','track','browser')):continue  
            fields = line.split()
            if(len(fields)<12):
                print("Invalid bed line (skipped):",line, end=' ', file=sys.stderr)
                continue
            chrom     = fields[0].upper()
            chrom_list.add(chrom)
            tx_start = int( fields[1] )
            #tx_end   = int( fields[2] )
            if int(fields[9]) == 1:
                continue        
            
            exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
            exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
            exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
            exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
            intron_start = exon_ends[:-1]
            intron_end=exon_starts[1:]
            for st,end in zip (intron_start, intron_end):
                knownSpliceSites.add(chrom + ":" + str(st) + "-" + str(end))
        print("Done! Total "+str(len(knownSpliceSites)) + " known splicing junctions.", file=sys.stderr)


        #read SAM file
        samSpliceSites=[]
        intron_start=[]
        intron_end=[]
        uniqSpliceSites=collections.defaultdict(int)

        if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
        else:print("Load SAM file ... ", end=' ', file=sys.stderr)
        try:
            while(1):
                aligned_read = next(self.samfile)
                if aligned_read.reference_id < 0:
                    continue                                    #no reference (unmapped)
                chrom = aligned_read.reference_name.upper()
                if chrom not in chrom_list:
                    continue                
                if not readfilter.passes(aligned_read):continue   #skip unmapped, secondary, supplementary, QC-failed and duplicate records
                if aligned_read.mapq < q_cut: continue
                
                hit_st = aligned_read.pos
                intron_blocks = bam_cigar.fetch_intron(chrom, hit_st, aligned_read.cigar)           
                if len(intron_blocks)==0:
                    continue
                for intrn in intron_blocks:
                    if intrn[2] - intrn[1] < min_intron:continue
                    samSpliceSites.append(intrn[0] + ":" + str(intrn[1]) + "-" + str(intrn[2]))
        except StopIteration:
            print("Done", file=sys.stderr)
        
        print("shuffling alignments ...", end=' ', file=sys.stderr)
        random.Random(seed).shuffle(samSpliceSites)
        print("Done", file=sys.stderr)
                
        #resampling
        SR_num = len(samSpliceSites)
        sample_size=0
        all_junctionNum = 0 
        known_junc=[]
        all_junc=[]
        unknown_junc=[]
        #=========================sampling uniquely mapped reads from population
        tmp = sampling_percentages(sample_start, sample_end, sample_step)
        index_st = 0
        for pertl in tmp:
            index_end = int(SR_num * (pertl/100.0))
            sample_size = index_end

            print("sampling " + str(pertl) +"% (" + str(sample_size) + ") splicing reads.", end=' ', file=sys.stderr)

            for i in range(index_st, index_end):
                uniqSpliceSites[samSpliceSites[i]] +=1
            index_st = index_end

            supported = [sj for sj, n in uniqSpliceSites.items() if n >= recur]
            all_junctionNum = len(supported)
            known_junctionNum = sum(1 for sj in supported if sj in knownSpliceSites)
            unknown_junctionNum = all_junctionNum - known_junctionNum
            all_junc.append(str(all_junctionNum))
            known_junc.append(str(known_junctionNum))
            unknown_junc.append(str(unknown_junctionNum))
            print(str(all_junctionNum) + " splicing junctions. " + str(known_junctionNum)
                  + " known splicing junctions. " + str(unknown_junctionNum)
                  + " novel splicing junctions.", file=sys.stderr)

        #for j in uniq_SJ:
            #print >>OUT, j + "\t" + str(uniq_SJ[j])
        print("pdf(%s)" % r_string(outfile + '.junctionSaturation_plot.pdf'), file=OUT)
        print("x=c(" + ','.join([str(i) for i in tmp]) + ')', file=OUT)
        print("y=c(" + ','.join(known_junc) + ')', file=OUT)
        print("z=c(" + ','.join(all_junc) + ')', file=OUT)
        print("w=c(" + ','.join(unknown_junc) + ')', file=OUT)
        print("m=max(%d,%d,%d)" % (int(int(known_junc[-1])/1000), int(int(all_junc[-1])/1000),int(int(unknown_junc[-1])/1000)), file=OUT)
        print("n=min(%d,%d,%d)" % (int(int(known_junc[0])/1000), int(int(all_junc[0])/1000),int(int(unknown_junc[0])/1000)), file=OUT)
        print("plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))", file=OUT)
        print("points(x,y/1000,type='o',col='red')", file=OUT)
        print("points(x,w/1000,type='o',col='green')", file=OUT)
        print('legend(5,%d, legend=c("All junctions","known junctions", "novel junctions"),col=c("blue","red","green"),lwd=1,pch=1)' % int(int(all_junc[-1])/1000), file=OUT)
        print("dev.off()", file=OUT)
        OUT.close()
        return {
            "percent": tmp,
            "all": [int(i) for i in all_junc],
            "known": [int(i) for i in known_junc],
            "novel": [int(i) for i in unknown_junc],
        }

    def saturation_RPKM(self,refbed,outfile,sample_start=5,sample_step=5,sample_end=100,skip_multi=True, strand_rule=None, q_cut=30, seed=None):
        '''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
        
        if refbed is None:
            raise ValueError("You must specify a bed file representing gene model")
        rpkm_file = outfile + ".eRPKM.xls"
        raw_file = outfile + ".rawCount.xls"
        
        RPKM_OUT = open(rpkm_file,'w')
        RAW_OUT = open(raw_file ,'w')
        
        ranges={}
        #totalReads=0
        cUR_num = 0 #number of fragements
        cUR_plus = 0
        cUR_minus = 0
        block_list_plus = []    #non-spliced read AS IS, splicing reads were counted multiple times
        block_list_minus = []
        block_list = []
        strandRule = readfilter.parse_strand_rule(strand_rule)


        #read SAM or BAM
        if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
        else:print("Load SAM file ... ", end=' ', file=sys.stderr)
        try:
            while(1):
                aligned_read = next(self.samfile)
                if not readfilter.passes(aligned_read):continue   #skip unmapped, secondary, supplementary, QC-failed and duplicate records
                
                if skip_multi:
                    if aligned_read.mapq < q_cut:
                        continue            
                chrom = self.samfile.getrname(aligned_read.tid).upper()
                
                #determine read_id and read_strand
                if aligned_read.is_paired:                      #pair end
                    if aligned_read.is_read1:read_id = '1'
                    if aligned_read.is_read2:read_id = '2'
                else:
                	read_id = ''                               #single end
            
                if aligned_read.is_reverse:
                	map_strand = '-'
                else:
                	map_strand = '+'               
                strand_key = read_id + map_strand               #used to determine if a read should assign to gene(+) or gene(-)

                hit_st = aligned_read.pos
                exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)   
                cUR_num += len(exon_blocks)
                
                #strand specific
                if strand_rule is not None:
                    gene_strand = strandRule.get(strand_key)     # None: read not covered by the rule
                    if gene_strand == '+': cUR_plus += len(exon_blocks)
                    if gene_strand == '-': cUR_minus += len(exon_blocks)
                    for exn in exon_blocks:
                        if gene_strand == '+': block_list_plus.append(exn[0] + ":" + str(exn[1] + int((exn[2]-exn[1])/2) ))
                        if gene_strand == '-': block_list_minus.append(exn[0] + ":" + str(exn[1] + int((exn[2]-exn[1])/2) ))
                #Not strand specific
                else:           
                    for exn in exon_blocks:
                        block_list.append(exn[0] + ":" + str(exn[1] + int((exn[2]-exn[1])/2) ))
        except StopIteration:
            print("Done", file=sys.stderr)
        
        
        print("shuffling alignments ...", end=' ', file=sys.stderr)
        rng = random.Random(seed)
        rng.shuffle(block_list_plus)
        rng.shuffle(block_list_minus)
        rng.shuffle(block_list)
        print("Done", file=sys.stderr)
        
        
        ranges_plus={}
        ranges_minus={}
        ranges={}
        sample_size=0
        RPKM_table=collections.defaultdict(list)
        rawCount_table=collections.defaultdict(list)
        RPKM_head=['#chr','start','end','name','score','strand']

        tmp = sampling_percentages(sample_start, sample_end, sample_step)
        percent_st = 0.0
        #=========================sampling uniquely mapped reads from population
        for pertl in tmp:
            percent_end = pertl/100.0
            sample_size = cUR_num * percent_end
            RPKM_head.append(str(pertl) + '%')
            
            if strand_rule is not None:
                print("sampling " + str(pertl) +"% (" + str(int(cUR_plus * percent_end)) + ") forward strand fragments ...", file=sys.stderr)
                for i in block_list_plus[int(cUR_plus*percent_st):int(cUR_plus*percent_end)]:
                    (chr,coord) = i.split(':')
                    if chr not in ranges_plus:ranges_plus[chr] = Intersecter()
                    ranges_plus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                               
                
                print("sampling " + str(pertl) +"% (" + str(int(cUR_minus * percent_end)) + ") reverse strand fragments ...", file=sys.stderr)          
                for i in block_list_minus[int(cUR_minus*percent_st):int(cUR_minus*percent_end)]:
                    (chr,coord) = i.split(':')
                    if chr not in ranges_minus:ranges_minus[chr] = Intersecter()                
                    ranges_minus[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                      
            
            else:
                print("sampling " + str(pertl) + "% (" + str(int(sample_size)) + ") fragments ...", file=sys.stderr)
                for i in block_list[int(cUR_num*percent_st):int(cUR_num*percent_end)]:
                    (chr,coord) = i.split(':')
                    if chr not in ranges:ranges[chr] = Intersecter()                        
                    ranges[chr].add_interval( Interval( int(coord), int(coord)+1 ) )                                                        
            percent_st = percent_end

            #========================= calculating RPKM based on sub-population
            print("assign reads to transcripts in " + refbed + ' ...', file=sys.stderr)
            for line in open(refbed,'r'):
                try:
                    if line.startswith(('#','track','browser')):continue  
                    # Parse fields from gene tabls
                    fields = line.split()
                    chrom     = fields[0].upper()
                    tx_start  = int( fields[1] )
                    tx_end    = int( fields[2] )
                    geneName      = fields[3]
                    strand    = fields[5]
                    exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                    exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                    exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                    exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
                    exon_sizes = list(map(int,fields[10].rstrip(',\n').split(',')))
                    key='\t'.join((fields[0],str(tx_start),str(tx_end),geneName,'0',strand))   # chromosome name as in the BED file
                except (IndexError, ValueError):
                    print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                    continue
                mRNA_count=0    #we need to initializ it to 0 for each gene
                mRNA_len=sum(exon_sizes)
                for st,end in zip(exon_starts,exon_ends):
                    #if chrom in ranges:
                    if strand_rule is not None:
                        if (strand == '+') and (chrom in ranges_plus): mRNA_count += len(ranges_plus[chrom].find(st,end))   
                        if (strand == '-') and (chrom in ranges_minus): mRNA_count += len(ranges_minus[chrom].find(st,end))
                    else:
                        if chrom in ranges:
                            mRNA_count += len(ranges[chrom].find(st,end))
                if mRNA_len ==0:
                    raise ValueError(geneName + " has 0 nucleotides")
                if sample_size == 0:
                    raise ValueError("Too few reads to sample")
                mRNA_RPKM = (mRNA_count * 1000000000.0)/(mRNA_len * sample_size)
                RPKM_table[key].append(str(mRNA_RPKM))
                rawCount_table[key].append(str(mRNA_count))
            print("", file=sys.stderr)

        #self.f.seek(0)
        print('\t'.join(RPKM_head), file=RPKM_OUT)
        print('\t'.join(RPKM_head), file=RAW_OUT)
        for key in RPKM_table:
            print(key + '\t', end=' ', file=RPKM_OUT)
            print('\t'.join(RPKM_table[key]), file=RPKM_OUT)
            print(key + '\t', end=' ', file=RAW_OUT)
            print('\t'.join(rawCount_table[key]), file=RAW_OUT)     


    def mismatchProfile(self,read_length,read_num, outfile, q_cut=30):
        '''
        Calculate mismatch profile. Note that the "MD" tag must exist.
        '''
        
        DOUT = open(outfile + '.mismatch_profile.xls','w')
        ROUT = open(outfile + '.mismatch_profile.r','w')
        
        #reading input SAM file
        if self.bam_format:print("Process BAM file ... ", end=' ', file=sys.stderr)
        else:print("Process SAM file ... ", end=' ', file=sys.stderr)
        
        MD_pat =re.compile(r'(\d+)([A-Z]+)')
        #number_base = re.compile(r'([0-9]+)([A-Z]+)', re.I)
        
        count = 0
        data = collections.defaultdict(dict)    # data[read_coord][genotype] = geno_type_number
        try:
            while(1):
                if count >= read_num:
                    print("Total reads used: " + str(count), file=sys.stderr)
                    break
                aligned_read = next(self.samfile)
                if not readfilter.passes(aligned_read):continue   #skip unmapped, secondary, supplementary, QC-failed and duplicate records
                if aligned_read.mapq < q_cut: continue
                if aligned_read.is_reverse:
                    strand = '-'
                else:
                    strand = '+'
                
                #Skip if there is no mismatch, or there is deletion
                tags = aligned_read.tags
                skip = False
                for tag in  tags:
                    if tag[0] == "NM" and tag[1] ==0:
                        skip = True     #skip reads with no mismatches
                    if (tag[0] == "MD") and ('^' in tag[1]):
                        skip = True     # skip as there is deletion from the reference
                if skip is True: continue
                
                # skip partially mapped read
                read_seq = aligned_read.seq
                if len(read_seq) != read_length:
                    continue
                if 'N' in read_seq:
                    continue
                
                matched_portion_size = 0
                for op,value in aligned_read.cigar:
                    if op == 0: matched_portion_size += value
                if matched_portion_size != read_length:
                    continue
                
                count += 1
                if strand == '+':
                    for tag in tags:
                        if tag[0] == "MD":
                            a = MD_pat.findall(tag[1])  #tag[1] = "5G19T75"; a = [('5', 'G'), ('19', 'T')]
                            read_coord = 0
                            for (match_number, ref_base) in a:
                                read_coord += int(match_number)
                                read_base = read_seq[read_coord]
                                if read_base == ref_base: continue
                                genotype = ref_base + '2' + read_base
                                if genotype not in data[read_coord]:
                                    data[read_coord][genotype] = 1
                                else:
                                    data[read_coord][genotype] += 1
                                read_coord += 1                 
                if strand == '-':
                    for tag in tags:
                        if tag[0] == "MD":
                            a = MD_pat.findall(tag[1])  #tag[1] = "5G19T75"; a = [('5', 'G'), ('19', 'T')]
                            read_coord = 0
                            for (match_number, ref_base) in a:
                                read_coord += int(match_number)
                                read_base = read_seq[read_coord]
                                if read_base == ref_base: continue
                                genotype = ref_base + '2' + read_base
                                if genotype not in data[read_length - read_coord -1]:
                                    data[read_length - read_coord -1][genotype] = 1
                                else:
                                    data[read_length - read_coord -1][genotype] += 1
                                read_coord += 1
                                if read_base == ref_base: print(aligned_read)               
        
        except StopIteration:
            print("Total reads used: " + str(count), file=DOUT)
        print('\n')
        
        if len(data) == 0:
            raise NoDataFound("No mismatches found")
        # write data out
        all_genotypes = ['A2C','A2G','A2T','C2A','C2G','C2T','G2A','G2C','G2T','T2A','T2C','T2G']
        print("read_pos\tsum\t" + '\t'.join(all_genotypes), file=DOUT)
        for indx in sorted(data):
            tmp = [indx, sum(data[indx].values())]  #read position and sum of mismatches
            for i in all_genotypes:
                if i in data[indx]:
                    tmp.append(data[indx][i])
                else:
                    tmp.append(0)
            print('\t'.join([str(i) for i in tmp]), file=DOUT)
        
        DOUT.close()
        
        # write Rscript
        r_data = collections.defaultdict(list)
        for gt in all_genotypes:
            for indx in sorted(data):
                if gt in data[indx]:
                    r_data[gt].append(data[indx][gt])
                else:
                    r_data[gt].append(0)
        for k in sorted(r_data):
            print('%s=c(%s)' % (k, ','.join([str(i) for i in r_data[k]])), file=ROUT)
        
        #print >>ROUT, 'color_code = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")'
        print('color_code = c("green","powderblue","lightseagreen","red","violetred4","mediumorchid1","blue","royalblue","steelblue1","orange","gold","black")', file=ROUT)
        
        print('y_up_bound = max(c(%s))' % (','.join(["log10(" + str(i) + "+1)" for i in all_genotypes])), file=ROUT)
        print('y_low_bound = min(c(%s))' % (','.join(["log10(" + str(i) + "+1)" for i in all_genotypes])), file=ROUT)
        
        print("pdf(%s)" % r_string(outfile + '.mismatch_profile.pdf'), file=ROUT)
        count=1
        for gt in all_genotypes:
            if count == 1:
                print("plot(log10(%s+1),type=\"l\",col=color_code[%d],ylim=c(y_low_bound,y_up_bound),ylab=\"log10(# of mismatch)\",xlab=\"Read position (5\'->3\')\")" % (gt,count), file=ROUT)
            else:
                print("lines(log10(%s+1), col=color_code[%d])" % (gt,count), file=ROUT)
            count += 1
        print("legend(13,y_up_bound,legend=c(%s), fill=color_code, border=color_code, ncol=4)" % (','.join(['"' + i + '"' for i in all_genotypes])), file=ROUT)
        print("dev.off()", file=ROUT)
        ROUT.close()
        return {"positions": sorted(data), "counts": {gt: list(r_data[gt]) for gt in all_genotypes}}

    def deletionProfile(self,read_length,read_num, outfile, q_cut=30):
        '''
        Calculate deletion profile. 
        Deletion: Deletion from the read (relative to the reference), CIGAR operator 'D'
        '''
        
        DOUT = open(outfile + '.deletion_profile.txt','w')
        ROUT = open(outfile + '.deletion_profile.r','w')
        
        #reading input SAM file
        if self.bam_format:
            print("Process BAM file ... ", end=' ', file=sys.stderr)
        else:
            print("Process SAM file ... ", end=' ', file=sys.stderr)
        
        count = 0
        del_postns = collections.defaultdict(int)   #key: position of read. value: deletion times
        #del_sizes = collections.defaultdict(int)   #key: deletion size. value: deletion frequency of this size
        try:
            while(1):
                if count >= read_num:
                    print("Total reads used: " + str(count), file=sys.stderr)
                    break
                aligned_read = next(self.samfile)
                if not readfilter.passes(aligned_read):continue   #skip unmapped, secondary, supplementary, QC-failed and duplicate records
                if aligned_read.mapq < q_cut: continue
                if aligned_read.is_reverse:
                    strand = '-'
                else:
                    strand = '+'
                
                # skip if read doesn't have deletion
                read_cigar = aligned_read.cigar
                if 2 not in [i[0] for i in read_cigar]:     #read contains no deletion
                    continue

                # skip partially mapped read
                read_seq = aligned_read.seq
                if len(read_seq) != read_length:
                    continue
                
                matched_portion_size = 0
                for op,value in aligned_read.cigar:
                    if op == 0: matched_portion_size += value   #match
                    if op == 4: matched_portion_size += value   #soft clp
                    if op == 1: matched_portion_size += value   #insertion to read
                if matched_portion_size != read_length:
                    continue
                
                count += 1
                del_positions = bam_cigar.fetch_deletion_range(read_cigar)  #[(position, size),(position, size),...]
                for (p,s) in del_positions:
                    if strand == '-': p = read_length - p
                    del_postns[p] += 1
                    #del_sizes[s] += 1
                #print aligned_read
                #print del_positions
                        
        except StopIteration:
            print("Total reads used: " + str(count), file=sys.stderr)
        print('\n')
        
        del_count=[]
        print("read_position\tdeletion_count", file=DOUT)
        for k in range(0,read_length):
            if k in del_postns:
                print(str(k) + '\t' + str(del_postns[k]), file=DOUT)
                del_count.append(str(del_postns[k]))
            else:
                print(str(k) + '\t0', file=DOUT)
                del_count.append('0')
        DOUT.close()
        
        print("pdf(%s)" % r_string(outfile + '.deletion_profile.pdf'), file=ROUT)
        print("pos=c(%s)" % ','.join([str(i) for i in range(0,read_length)]), file=ROUT)
        print("value=c(%s)" % ','.join([i for i in del_count]), file=ROUT)      
        print("plot(pos,value,type='b', col='blue',xlab=\"Read position (5\'->3\')\", ylab='Deletion count')", file=ROUT)                   
        print("dev.off()", file=ROUT)
        ROUT.close()
        return {"positions": list(range(read_length)), "counts": [int(c) for c in del_count]}

def print_bits_as_bed( bits ):
    end = 0
    while 1:
        start = bits.next_set( end )
        if start == bits.size: break
        end = bits.next_clear( start )
        print("%d\t%d" % ( start, end ))





























