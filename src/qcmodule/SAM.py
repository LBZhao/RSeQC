#!/usr/bin/env python

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

# -*- coding: utf-8 -*-
'''manipulate BAM/SAM file.'''
#import built-in modules
import os,sys
import re
import string
import collections
import random

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from bx.binned_array import BinnedArray
from bx_extras.fpconst import isNaN
from bx.bitset_utils import *
import pysam
from qcmodule import mystat
from qcmodule import fasta
from qcmodule import bam_cigar
from qcmodule import readfilter
from qcmodule import BED
from importlib.metadata import version

__author__ = "Liguo Wang"
__version__ = version("RSeQC")


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


def sampling_percentages(sample_start, sample_end, sample_step):
    """Percentages at which saturation curves are evaluated (always ends at 100).

    ``range(sample_start, sample_end, sample_step)`` plus 100. Each step adds
    the reads between the previous percentage and this one, so every read is
    used exactly once by 100%. (Before RSeQC 5.1.0 each step added the slice
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
        try:
            self.samfile = pysam.Samfile(inputFile,'rb')
            if len(self.samfile.header) ==0:
                print("BAM/SAM file has no header section. Exit!", file=sys.stderr)
                sys.exit(1)
            self.bam_format = True
        except:
            self.samfile = pysam.Samfile(inputFile,'r')
            if len(self.samfile.header) ==0:
                print("BAM/SAM file has no header section. Exit!", file=sys.stderr)
                sys.exit(1)
            self.bam_format = False

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
        print("%-40s%d" % ("Total records:",R_total), file=sys.stdout)
        print("\n", end='', file=sys.stdout)
        print("%-40s%d" % ("QC failed:",R_qc_fail), file=sys.stdout)
        print("%-40s%d" % ("Optical/PCR duplicate:", R_duplicate), file=sys.stdout)
        print("%-40s%d" % ("Non primary hits", R_nonprimary), file=sys.stdout)
        print("%-40s%d" % ("Supplementary hits:", R_supplementary), file=sys.stdout)
        print("%-40s%d" % ("Unmapped reads:",R_unmap), file=sys.stdout)
        print("%-40s%d" % ("mapq < mapq_cut (non-unique):",R_multipleHit), file=sys.stdout)
        print("\n", end='', file=sys.stdout)
        print("%-40s%d" % ("mapq >= mapq_cut (unique):",R_uniqHit), file=sys.stdout)
        print("%-40s%d" % ("Read-1:",R_read1), file=sys.stdout)
        print("%-40s%d" % ("Read-2:",R_read2), file=sys.stdout)
        print("%-40s%d" % ("Reads map to '+':",R_forward), file=sys.stdout)
        print("%-40s%d" % ("Reads map to '-':",R_reverse), file=sys.stdout)
        print("%-40s%d" % ("Non-splice reads:",R_nonSplice), file=sys.stdout)
        print("%-40s%d" % ("Splice reads:",R_splice), file=sys.stdout)  
        print("%-40s%d" % ("Reads mapped in proper pairs:",R_properPair), file=sys.stdout)
        print("%-40s%d" % ("Proper-paired reads map to different chrom:",R_pair_diff_chrom), file=sys.stdout)
    
    def configure_experiment(self,refbed,sample_size, q_cut = 30):
        '''Given a BAM/SAM file, this function will try to guess the RNA-seq experiment:
            1) single-end or pair-end
            2) strand_specific or not
            3) if it is strand-specific, what's the strand_ness of the protocol
        '''
        
            #how many reads you want to sample
        count =0
        p_strandness=collections.defaultdict(int)
        s_strandness=collections.defaultdict(int)
        #load reference gene model
        gene_ranges={}
        print("Reading reference gene model " + refbed + ' ...', end=' ', file=sys.stderr)
        for line in open(refbed,'r'):
            try:
                if line.startswith(('#','track','browser')):continue  
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0]
                txStart  = int( fields[1] )
                txEnd    = int( fields[2] )
                #geneName      = fields[3]
                strand    = fields[5]
            except:
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, file=sys.stderr)
                continue
            if chrom not in gene_ranges:
                gene_ranges[chrom]=Intersecter()
            gene_ranges[chrom].insert(txStart,txEnd,strand)                         
        print("Done", file=sys.stderr)      
        
        #read SAM/BAM file
        #current_pos = self.samfile.tell()
        print("Loading SAM/BAM file ... ", end=' ', file=sys.stderr)
        try:
            while(1):
                if count >= sample_size:
                    break
                aligned_read = next(self.samfile)
                if not readfilter.passes(aligned_read, min_mapq=q_cut):
                    continue                        #unmapped, secondary, supplementary, QC-failed, duplicate or low MAPQ
                
                chrom = self.samfile.getrname(aligned_read.tid)
                if aligned_read.is_paired:
                    if aligned_read.is_read1:
                        read_id = '1'
                    if aligned_read.is_read2:
                        read_id = '2'
                    if aligned_read.is_reverse:
                        map_strand = '-'
                    else:
                        map_strand = '+'
                    readStart = aligned_read.pos
                    readEnd = readStart + aligned_read.qlen
                    if chrom in gene_ranges:
                        tmp = set(gene_ranges[chrom].find(readStart,readEnd))
                        if len(tmp) == 0: continue
                        strand_from_gene = ':'.join(tmp)
                        p_strandness[read_id + map_strand + strand_from_gene]+=1    
                        count += 1
                else:
                    if aligned_read.is_reverse:
                        map_strand = '-'
                    else:
                        map_strand = '+'                    
                    readStart = aligned_read.pos
                    readEnd = readStart + aligned_read.qlen
                    if chrom in gene_ranges:
                        tmp = set(gene_ranges[chrom].find(readStart,readEnd))
                        if len(tmp) == 0: continue
                        strand_from_gene = ':'.join(tmp)
                        s_strandness[map_strand + strand_from_gene]+=1
                        count += 1

        except StopIteration:
            print("Finished", file=sys.stderr)      
        #self.samfile.seek(current_pos)
        
        print("Total " + str(count) + " usable reads were sampled", file=sys.stderr)
        protocol="unknown"
        #strandness=None
        spec1=0.0
        spec2=0.0
        other=0.0
        if len(p_strandness) >0 and len(s_strandness) ==0 :
            protocol="PairEnd"
            #for k,v in p_strandness.items():
            #   print >>sys.stderr, k + '\t' + str(v)
            spec1= (p_strandness['1++'] + p_strandness['1--'] + p_strandness['2+-'] + p_strandness['2-+'])/float(sum(p_strandness.values()))
            spec2= (p_strandness['1+-'] + p_strandness['1-+'] + p_strandness['2++'] + p_strandness['2--'])/float(sum(p_strandness.values()))
            other = 1-spec1-spec2
            
        elif len(s_strandness) >0 and len(p_strandness) ==0 :
            protocol="SingleEnd"
            #for k,v in s_strandness.items():
            #   print  >>sys.stderr, k + '\t' + str(v)
            spec1 = (s_strandness['++'] + s_strandness['--'])/float(sum(s_strandness.values()))
            spec2 = (s_strandness['+-'] + s_strandness['-+'])/float(sum(s_strandness.values()))
            other = 1-spec1-spec2
        else:
            protocol="Mixture"
            spec1 = 0
            spec2 = 0
            other = 0
        return [protocol,spec1,spec2,other]

    def bamTowig(self,outfile,chrom_sizes, chrom_file,skip_multi=True,strand_rule=None,WigSumFactor=None,q_cut=30):
        """Convert BAM/SAM file to wig file. chrom_size is dict with chrom as key and chrom_size as value
        strandRule should be determined from \"infer_experiment\". such as \"1++,1--,2+-,2-+\". When
        WigSumFactor is provided, output wig file will be normalized to this number """
        
        #strand_rule={'1+':'-','1-':'+','2+':'+','2-':'-'}
        strandRule={}
        if strand_rule is None:                                                 # Not strand-specific
            pass                                                                
        elif len(strand_rule.split(',')) ==4:                                   #PairEnd, strand-specific
            for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
        elif len(strand_rule.split(',')) ==2:                                   #singeEnd, strand-specific
            for i in strand_rule.split(','):strandRule[i[0]]=i[1]
        else:
            print("Unknown value of option :'strand_rule' " + strand_rule, file=sys.stderr)
            sys.exit(1)
        if len(strandRule) == 0:
            FWO = open(outfile + '.wig','w')
        else:
            FWO = open(outfile + '.Forward.wig','w')
            RVO = open(outfile + '.Reverse.wig','w')
        
        read_id=''
        
        for chr_name, chr_size in list(chrom_sizes.items()):        #iterate each chrom
            try:
                self.samfile.fetch(chr_name,0,chr_size)
            except:
                print("No alignments for " + chr_name + '. skipped', file=sys.stderr)
                continue
            print("Processing " + chr_name + " ...", file=sys.stderr)
            if len(strandRule) == 0: FWO.write('variableStep chrom='+chr_name+'\n')
            else:
                FWO.write('variableStep chrom='+chr_name+'\n')
                RVO.write('variableStep chrom='+chr_name+'\n')
            Fwig = collections.defaultdict(int)
            Rwig = collections.defaultdict(int)
            alignedReads = self.samfile.fetch(chr_name,0,chr_size)      
            for aligned_read in alignedReads:
                if not readfilter.passes(aligned_read):continue   #skip unmapped, secondary, supplementary, QC-failed and duplicate records
                
                if skip_multi:
                    if aligned_read.mapq < q_cut:
                        continue                
                if aligned_read.is_paired:
                    if aligned_read.is_read1:read_id = '1'
                    if aligned_read.is_read2:read_id = '2'
                
                if aligned_read.is_reverse:map_strand = '-'
                else:map_strand = '+'
                
                key = read_id + map_strand
                
                hit_st = aligned_read.pos
                for block in bam_cigar.fetch_exon(chr_name, hit_st, aligned_read.cigar): 
                    for pos in range(block[1]+1,block[2]+1):    
                        if len(strandRule) == 0: Fwig[pos] +=1.0    #this is NOT strand specific. everything into Fwig
                        else:                                       #this is strand specific. separate Fwig and Rwig
                            if strandRule[key] == '+':Fwig[pos] +=1.0
                            if strandRule[key] == '-':Rwig[pos] -=1.0
            if WigSumFactor is None:    #not normalize  
                if len(strandRule) == 0:                            #this is NOT strand specific.
                    for pos in sorted (Fwig.keys()):
                        print("%d\t%.2f" % (pos,Fwig[pos]), file=FWO)
                else:
                    for pos in sorted (Fwig.keys()):
                        print("%d\t%.2f" % (pos,Fwig[pos]), file=FWO)
                    for pos in sorted (Rwig.keys()):
                        print("%d\t%.2f" % (pos,Rwig[pos]), file=RVO)
            else:                       #normalize wig signal to WigSumFactor
                if len(strandRule) == 0:                            #this is NOT strand specific.
                    for pos in sorted (Fwig.keys()):
                        print("%d\t%.2f" % (pos,Fwig[pos]*WigSumFactor), file=FWO)
                else:
                    for pos in sorted (Fwig.keys()):
                        print("%d\t%.2f" % (pos,Fwig[pos]*WigSumFactor), file=FWO)
                    for pos in sorted (Rwig.keys()):
                        print("%d\t%.2f" % (pos,Rwig[pos]*WigSumFactor), file=RVO)
        FWO.close()
        if len(strandRule) != 0:
            RVO.close()
        if len(strandRule) == 0:
            try:
                import subprocess
                print("Run " + "wigToBigWig " + outfile + '.wig ' + chrom_file + ' ' +  outfile + ".bw ")
                subprocess.call("wigToBigWig -clip " + outfile + '.wig ' + chrom_file + ' ' +  outfile + ".bw ",shell=True)
            except:
                print("Failed to call \"wigToBigWig\".", file=sys.stderr)
                pass
        else:
            try:
                import subprocess
                subprocess.call("wigToBigWig -clip " + outfile + '.Forward.wig ' + chrom_file + ' ' +  outfile + ".Forward.bw ",shell=True)
                subprocess.call("wigToBigWig -clip " + outfile + '.Reverse.wig ' + chrom_file + ' ' +  outfile + ".Reverse.bw ",shell=True)
            except:
                print("Failed to call \"wigToBigWig\".", file=sys.stderr)
                pass            

    def calWigSum(self,chrom_sizes, skip_multi=True):
        """Calculate wigsum from BAM file"""
        
        print("Calcualte wigsum ... ", file=sys.stderr)
        wigsum = 0.0
        #read_id=''
        for chr_name, chr_size in list(chrom_sizes.items()):        #iterate each chrom
            try:
                self.samfile.fetch(chr_name,0,chr_size)
            except:
                print("No alignments for " + chr_name + '. skipped', file=sys.stderr)
                continue
            print("Processing " + chr_name + " ...", file=sys.stderr)

            alignedReads = self.samfile.fetch(chr_name,0,chr_size)      
            for aligned_read in alignedReads:
                flag=0
                if not readfilter.passes(aligned_read):continue   #skip unmapped, secondary, supplementary, QC-failed and duplicate records
                
                if skip_multi:
                    if len(aligned_read.tags)>0:        #( ("NM", 1),("RG", "L1") )
                        for i in aligned_read.tags:
                            if i[0] in ParseBAM.multi_hit_tags and i[1] >1:
                                flag=1                      #multiple hit read
                                break
                    if flag==1:continue                     #skip multiple map read     
                
                #if aligned_read.is_paired:
                    #if aligned_read.is_read1:read_id = '1'
                    #if aligned_read.is_read2:read_id = '2'
                
                #if aligned_read.is_reverse:map_strand = '-'
                #else:map_strand = '+'
                
                #key = read_id + map_strand
                
                hit_st = aligned_read.pos
                for block in bam_cigar.fetch_exon(chr_name, hit_st, aligned_read.cigar): 
                    wigsum += (block[2] - block[1])
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

        try:
            while(1):
                aligned_read = next(self.samfile)
                if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                if aligned_read.mapq < q_cut: continue
                #if aligned_read.is_unmapped:continue   #skip unmapped read
                #if aligned_read.is_qcfail:continue #skip low quality

                RNA_read = aligned_read.seq.upper()     
                if aligned_read.is_reverse:
                    RNA_read = RNA_read.translate(transtab)[::-1]
                for i,j in enumerate(RNA_read):
                    key = str(i) + j
                    base_freq[key] += 1
        except StopIteration:
            print("Done", file=sys.stderr)
        
        print("generating data matrix ...", file=sys.stderr)
        print("Position\tA\tC\tG\tT\tN\tX", file=FO)
        for i in range(len(RNA_read)):
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
        print("position=c(" + ','.join([str(i) for i in range(len(RNA_read))]) + ')', file=RS)
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
            
            print('pdf(\"%s\")' % (outfile +".NVC_plot.pdf"), file=RS)
            print('plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")', file=RS)
            print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
            print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
            print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
            print('lines(position,N_count/total,type="o",pch=20,col="black")', file=RS)     
            print('lines(position,X_count/total,type="o",pch=20,col="grey")', file=RS)  
            print('legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C","N","X"),col=c("dark green","red","blue","cyan","black","grey"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan","black","grey"))', file=RS)
            print("dev.off()", file=RS)
        else:
            print("total= A_count + C_count + G_count + T_count", file=RS)
            print("ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05", file=RS)
            print("yn=min(A_count/total,C_count/total,G_count/total,T_count/total)", file=RS)
        
            print('pdf(\"%s\")' % (outfile +".NVC_plot.pdf"), file=RS)
            print('plot(position,A_count/total,type="o",pch=20,ylim=c(yn,ym),col="dark green",xlab="Position of Read",ylab="Nucleotide Frequency")', file=RS)
            print('lines(position,T_count/total,type="o",pch=20,col="red")', file=RS)
            print('lines(position,G_count/total,type="o",pch=20,col="blue")', file=RS)
            print('lines(position,C_count/total,type="o",pch=20,col="cyan")', file=RS)
            print('legend('+ str(len(RNA_read)-10) + ',ym,legend=c("A","T","G","C"),col=c("dark green","red","blue","cyan"),lwd=2,pch=20,text.col=c("dark green","red","blue","cyan"))', file=RS)
            print("dev.off()", file=RS)
        
        RS.close()
        #self.f.seek(0)

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
        try:
            while(1):
                aligned_read = next(self.samfile)
                if readfilter.is_non_primary(aligned_read): continue   #skip secondary/supplementary records (same read counted twice)
                if aligned_read.mapq < q_cut: continue
                #if aligned_read.is_unmapped:continue   #skip unmapped read
                #if aligned_read.is_qcfail:continue     #skip low quality
                
                qual_str = aligned_read.qqual
                read_len = aligned_read.rlen
                if aligned_read.is_reverse:
                    qual_str = qual_str[::-1]

                for i,j in enumerate(qual_str):
                    q=ord(j)-33
                    if q > q_max: q_max = q
                    if q < q_min: q_min = q
                    try:
                        quality[i][q] += 1
                    except:
                        quality[i][q] = 1
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
        print("pdf(\'%s\')" % (outfile + ".qual.boxplot.pdf"), file=FO)
        for i in sorted(i_box):
            print('p'+str(i) + '<-' + i_box[i], file=FO)
        print('boxplot(' + ','.join(['p'+str(i) for i in i_box]) + ',xlab=\"Position of Read(5\'->3\')\",ylab=\"Phred Quality Score\",outline=F' + ')', file=FO)
        print("dev.off()", file=FO)
        
        
        #generate R script for heatmap
        print('\n', file=FO)
        print("pdf(\'%s\')" % (outfile + ".qual.heatmap.pdf"), file=FO)
        print("qual=c(" + ','.join(q_list)  + ')', file=FO)
        print("mat=matrix(qual,ncol=%s,byrow=F)" % (read_len), file=FO)
        print('Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c(\'spline\'))', file=FO)
        print("heatmap(mat,Rowv=NA,Colv=NA,xlab=\"Position of Read\",ylab=\"Phred Quality Score\",labRow=seq(from=%s,to=%s),col = Lab.palette(256),scale=\"none\" )" % (q_min,q_max), file=FO)
        print('dev.off()', file=FO)
        
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
        print("pdf(\"%s\")" % (outfile +  ".GC_plot.pdf"), file=RS)
        print('gc=rep(c(' + ','.join([i for i in list(gc_hist.keys())]) + '),' + 'times=c(' + ','.join([str(i) for i in list(gc_hist.values())]) + '))', file=RS)
        print('hist(gc,probability=T,breaks=%d,xlab="GC content (%%)",ylab="Density of Reads",border="blue",main="")' % 100, file=RS)
        #print >>RS, "lines(density(gc),col='red')"
        print("dev.off()", file=RS)     
        #self.f.seek(0)

    def readDupRate(self,q_cut, outfile=None,up_bound=500):
        '''Calculate reads's duplicate rates'''
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
        print("pdf(\'%s\')" % (outfile +".DupRate_plot.pdf"), file=RS)
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
        #self.f.seek(0)

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
                'pdf("%s")' % (outfile + ".clipping_profile.pdf"),
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
                'pdf("%s")' % (outfile + ".clipping_profile.R1.pdf"),
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
                'pdf("%s")' % (outfile + ".clipping_profile.R2.pdf"),
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
            
            
            print("pdf(\"%s\")" % (outfile + '.insertion_profile.pdf'), file=ROUT)
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
        
            
            print("pdf(\"%s\")" % (outfile + '.insertion_profile.R1.pdf'), file=ROUT)
            print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
            print("r1_insert_count=c(%s)" % ','.join([str(i) for i in r1_clip_count]), file=ROUT)
            print("r1_noninsert_count = %d - r1_insert_count" % (total_read1), file=ROUT)
            print('plot(read_pos, r1_insert_count*100/(r1_insert_count + r1_noninsert_count),col="blue",main="Insertion profile",xlab="Position of read (read-1)",ylab="Insertion %",type="b")', file=ROUT)
            print("dev.off()", file=ROUT)

            print("pdf(\"%s\")" % (outfile + '.insertion_profile.R2.pdf'), file=ROUT)
            print("read_pos=c(%s)" % ','.join([str(i) for i in read_pos]), file=ROUT)
            print("r2_insert_count=c(%s)" % ','.join([str(i) for i in r2_clip_count]), file=ROUT)
            print("r2_noninsert_count = %d - r2_insert_count" % (total_read2), file=ROUT)
            print('plot(read_pos, r2_insert_count*100/(r2_insert_count + r2_noninsert_count),col="blue",main="Insertion profile",xlab="Position of read (read-2)",ylab="Insertion %",type="b")', file=ROUT)
            print("dev.off()", file=ROUT)
    
    def mRNA_inner_distance(self,outfile,refbed,low_bound=0,up_bound=1000,step=10,sample_size=1000000, q_cut=30):
        '''estimate the inner distance of mRNA pair end fragment. fragment size = insert_size + 2 x read_length'''
        
        out_file1 = outfile + ".inner_distance.txt" 
        out_file2 = outfile + ".inner_distance_freq.txt"
        out_file3 = outfile + ".inner_distance_plot.r"
        
        FO=open(out_file1,'w')
        FQ=open(out_file2,'w')
        RS=open(out_file3,'w')
        
        fchrom="chr100"     #this is the fake chromosome
        ranges={}
        ranges[fchrom]=Intersecter()
        
        window_left_bound = list(range(low_bound,up_bound,step))
        #frag_size=0
        
        inner_distance_bitsets=BinnedBitSet()
        tmp = BinnedBitSet()
        tmp.set_range(0,0)
        pair_num=0
        sizes=[]
        counts=[]
        count=0
        
        print("Get exon regions from " + refbed + " ...", file=sys.stderr)
        bed_obj = BED.ParseBED(refbed)
        ref_exons = []
        
        for exn in bed_obj.getExon():
            ref_exons.append([exn[0].upper(), exn[1], exn[2]])
        exon_bitsets = binned_bitsets_from_list(ref_exons)
        
        transcript_ranges = {}
        for i_chr, i_st, i_end, i_strand, i_name in bed_obj.getTranscriptRanges():
            i_chr = i_chr.upper()
            if i_chr not in transcript_ranges:
                transcript_ranges[i_chr] = Intersecter()
            else:
                transcript_ranges[i_chr].add_interval(Interval(i_st, i_end, value=i_name))
        
        if self.bam_format:print("Load BAM file ... ", end=' ', file=sys.stderr)
        else:print("Load SAM file ... ", end=' ', file=sys.stderr)

        try:
            while(1):
                if pair_num >= sample_size:
                    break
                splice_intron_size=0
                aligned_read = next(self.samfile)
                if not readfilter.passes(aligned_read):continue   #skip unmapped, secondary, supplementary, QC-failed and duplicate records
                if not aligned_read.is_paired: continue     #skip single map read
                if aligned_read.mate_is_unmapped:continue   #
                if aligned_read.mapq < q_cut:continue
                
                read1_len = aligned_read.qlen
                read1_start = aligned_read.pos
                read2_start = aligned_read.mpos     #0-based, not included
                if read2_start < read1_start:
                    continue                                #because BAM file is sorted, mate_read is already processed if its coordinate is smaller
                if  read2_start == read1_start and aligned_read.is_read1:
                    inner_distance = 0
                    continue
                
                pair_num +=1

                # check if reads were mapped to diff chromsomes
                R_read1_ref = self.samfile.getrname(aligned_read.tid)
                R_read2_ref = self.samfile.getrname(aligned_read.rnext)
                if R_read1_ref != R_read2_ref:
                    FO.write(aligned_read.qname + '\t' + 'NA' + '\tsameChrom=No\n') #reads mapped to different chromosomes
                    continue
                
                
                chrom = self.samfile.getrname(aligned_read.tid).upper()
                intron_blocks = bam_cigar.fetch_intron(chrom, read1_start, aligned_read.cigar)              
                for intron in intron_blocks:
                    splice_intron_size += intron[2] - intron[1]
                read1_end = read1_start + read1_len + splice_intron_size        
                
                if read2_start >= read1_end:
                    inner_distance = read2_start - read1_end
                else:
                    exon_positions = []
                    exon_blocks = bam_cigar.fetch_exon(chrom, read1_start,aligned_read.cigar)
                    for ex in exon_blocks:
                        for i in range(ex[1]+1,ex[2]+1):
                            exon_positions.append(i)
                    inner_distance = -len([i for i in exon_positions if i > read2_start and i <= read1_end])
                #print aligned_read.qname,read1_end, read2_start
                
                read1_gene_names = set()    #read1_end
                try:
                    for gene in transcript_ranges[chrom].find(read1_end-1, read1_end):  #gene: Interval(0, 10, value=a)
                        read1_gene_names.add(gene.value)
                except:
                    pass
                
                read2_gene_names = set()    #read2_start
                try:
                    for gene in transcript_ranges[chrom].find(read2_start, read2_start +1): #gene: Interval(0, 10, value=a)
                        read2_gene_names.add(gene.value)
                except:
                    pass
                    
                if len(read1_gene_names.intersection(read2_gene_names)) == 0:   # no common gene
                    FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tsameTranscript=No,dist=genomic\n')        #reads mapped to different gene
                    ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )     
                    continue            

                if inner_distance > 0: 
                    if chrom in exon_bitsets:
                        size =0 
                        inner_distance_bitsets.set_range(read1_end, read2_start-read1_end)
                        inner_distance_bitsets.iand(exon_bitsets[chrom])
                        end=0
                        while 1:
                            start = inner_distance_bitsets.next_set( end )
                            if start == inner_distance_bitsets.size: break
                            end = inner_distance_bitsets.next_clear( start )
                            size += (end - start)
                        inner_distance_bitsets.iand(tmp)                                            #clear BinnedBitSet
                        
                        if size == inner_distance:
                            FO.write(aligned_read.qname + '\t' + str(size) + '\tsameTranscript=Yes,sameExon=Yes,dist=mRNA\n')
                            ranges[fchrom].add_interval( Interval( size-1, size ) )
                        elif size > 0 and size < inner_distance:
                            FO.write(aligned_read.qname + '\t' + str(size) + '\tsameTranscript=Yes,sameExon=No,dist=mRNA\n')
                            ranges[fchrom].add_interval( Interval( size-1, size ) ) 
                        elif size <= 0:
                            FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tsameTranscript=Yes,nonExonic=Yes,dist=genomic\n')
                            ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )     
                    else:
                        FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\tunknownChromosome,dist=genomic')
                        ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
                else:
                    FO.write(aligned_read.qname + '\t' + str(inner_distance) + '\treadPairOverlap\n')
                    ranges[fchrom].add_interval( Interval( inner_distance-1, inner_distance ) )
                
        except StopIteration:
            print("Done", file=sys.stderr)
        
        print("Total read pairs  used " + str(pair_num), file=sys.stderr)
        if pair_num==0:
            print("Cannot find paired reads", file=sys.stderr)
            sys.exit(0)
        
        for st in window_left_bound:
            sizes.append(str(st + step/2))
            count = str(len(ranges[fchrom].find(st,st + step)))
            counts.append(count)
            print(str(st) + '\t' + str(st+step) +'\t' + count, file=FQ)     
        print("out_file = \'%s\'" % outfile, file=RS)
        print("pdf(\'%s\')" % (outfile + ".inner_distance_plot.pdf"), file=RS)
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
        FO.close()
        FQ.close()
        RS.close()
        #self.f.seek(0)

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

        Before RSeQC 5.1.0 ``novel_combination`` junctions were reported as
        ``annotated``. Chromosome names are matched case-insensitively and
        written as they appear in the BAM header.
        """
        out_file = outfile + ".junction.xls"
        out_file2 = outfile + ".junction_plot.r"
        if refgene is None:
            print("You must provide reference gene model in bed format.", file=sys.stderr)
            sys.exit(1)

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
            print("No splice junction found.", file=sys.stderr)
            sys.exit()

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
            print("No splice read found", file=sys.stderr)
            sys.exit(1)
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
                print('pdf("%s")' % (outfile + pdf), file=ROUT)
                print("x=c(" + ','.join(str(v) for v in values) + ')', file=ROUT)
                print('pie(x,col=c(4,6,2,3),init.angle=30,main="%s",labels=c(%s))' % (title, labels), file=ROUT)
                print("dev.off()", file=ROUT)


    def saturation_junction(self,refgene,outfile=None,sample_start=5,sample_step=5,sample_end=100,min_intron=50,recur=1, q_cut=30, seed=None):
        '''check if an RNA-seq experiment is saturated in terms of detecting known splicing junction.

        Splice events are shuffled (with *seed*, for reproducible results) and
        added cumulatively at each sampling percentage. A junction is counted
        (all, known or novel) once it is supported by at least *recur* reads.
        '''
        
        out_file = outfile + ".junctionSaturation_plot.r"
        if refgene is None:
            print("You must provide reference gene model in bed format.", file=sys.stderr)
            sys.exit(1)
        
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
                try:
                    chrom = self.samfile.getrname(aligned_read.tid).upper()
                except:
                    continue
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
        print("pdf(\'%s\')" % (outfile + '.junctionSaturation_plot.pdf'), file=OUT)
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

    def saturation_RPKM(self,refbed,outfile,sample_start=5,sample_step=5,sample_end=100,skip_multi=True, strand_rule=None, q_cut=30, seed=None):
        '''for each gene, check if its RPKM (epxresion level) has already been saturated or not'''
        
        if refbed is None:
            print("You must specify a bed file representing gene model\n", file=sys.stderr)
            exit(0)
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
        strandRule = {}
                
        if strand_rule is None:                                                 # Not strand-specific
            pass                                                                
        elif len(strand_rule.split(',')) ==4:                                   #PairEnd, strand-specific
            for i in strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
        elif len(strand_rule.split(',')) ==2:                                   #singeEnd, strand-specific
            for i in strand_rule.split(','):strandRule[i[0]]=i[1]
        else:
            print("Unknown value of: 'strand_rule' " +  strand_rule, file=sys.stderr)
            sys.exit(1) 


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
                    if strandRule[strand_key] == '+': cUR_plus += len(exon_blocks)
                    if strandRule[strand_key] == '-': cUR_minus += len(exon_blocks)
                    for exn in exon_blocks:
                        if strandRule[strand_key] == '+': block_list_plus.append(exn[0] + ":" + str(exn[1] + int((exn[2]-exn[1])/2) ))
                        if strandRule[strand_key] == '-': block_list_minus.append(exn[0] + ":" + str(exn[1] + int((exn[2]-exn[1])/2) ))
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
                except:
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
                    print(geneName + " has 0 nucleotides. Exit!", file=sys.stderr)
                    sys.exit(1)
                if sample_size == 0:
                    print("Too few reads to sample. Exit!", file=sys.stderr)
                    sys.exit(1)
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
            print("No mismatches found", file=sys.stderr)
            sys.exit()
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
        
        print("pdf(\"%s\")" % (outfile + '.mismatch_profile.pdf'), file=ROUT)
        count=1
        for gt in all_genotypes:
            if count == 1:
                print("plot(log10(%s+1),type=\"l\",col=color_code[%d],ylim=c(y_low_bound,y_up_bound),ylab=\"log10(# of mismatch)\",xlab=\"Read position (5\'->3\')\")" % (gt,count), file=ROUT)
            else:
                print("lines(log10(%s+1), col=color_code[%d])" % (gt,count), file=ROUT)
            count += 1
        print("legend(13,y_up_bound,legend=c(%s), fill=color_code, border=color_code, ncol=4)" % (','.join(['"' + i + '"' for i in all_genotypes])), file=ROUT)
        print("dev.off()", file=ROUT)
            
            
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
        
        print("pdf(\"%s\")" % (outfile + '.deletion_profile.pdf'), file=ROUT)
        print("pos=c(%s)" % ','.join([str(i) for i in range(0,read_length)]), file=ROUT)
        print("value=c(%s)" % ','.join([i for i in del_count]), file=ROUT)      
        print("plot(pos,value,type='b', col='blue',xlab=\"Read position (5\'->3\')\", ylab='Deletion count')", file=ROUT)                   
        print("dev.off()", file=ROUT)


def print_bits_as_bed( bits ):
    end = 0
    while 1:
        start = bits.next_set( end )
        if start == bits.size: break
        end = bits.next_clear( start )
        print("%d\t%d" % ( start, end ))





























