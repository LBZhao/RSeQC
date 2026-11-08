#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

# -*- coding: utf-8 -*-
"""
Analzye 10X genomics single cell BAM files.
"""

import pysam
import sys
import collections
import logging


def _percent(count, total):
    """Return count as a percentage of total, or 0.0 when total is zero."""
    return 0.0 if total == 0 else count * 100.0 / total
import re
import pandas as pd

_GEM_SUFFIX = re.compile(r"-\d+$")


def strip_suffix(barcode):
    """Remove the Cell Ranger GEM-well suffix (``-1``, ``-2``, ...).

    Before RSeQC 5.1.0+lbzhao only ``-1`` was removed, so in aggregated BAM files
    (``-2``, ``-3``, ...) every barcode was counted as edited.
    """
    return _GEM_SUFFIX.sub("", barcode)


def diff_str(s1, s2):
    '''
    Comparing orignal barcode to the corrected barcode
    find the index and the nucleotide that has been corrected.


    Parameters
    ----------
    s1 : str
        the original barcode
    s2 : str
        the corrrected barcode

    '''
    results = []
    if len(s1) != len(s2):
        return results
    diff_positions = [i for i in range(len(s1)) if s1[i] != s2[i]]
    for pos in diff_positions:
        results.append([pos, s1[pos], s2[pos]])
    return results

def read_match_type(cigar_str):
    '''return the matching type between read and ref'''
    match_type = ''
    if bool(re.search(r'\A\d+M\Z', cigar_str)):
        match_type ='Map_consecutively'
    elif bool(re.search(r'\A\d+M\d+N\d+M\Z', cigar_str)):
        match_type = 'Map_with_splicing'
    elif bool(re.search(r'\A\d+S\d+M\Z', cigar_str)):
        match_type = 'Map_with_clipping'
    elif bool(re.search(r'\A\d+M\d+S\Z', cigar_str)):
        match_type = 'Map_with_clipping'
    elif bool(re.search(r'\A\d+M\d+N\d+M\d+S\Z', cigar_str)):
        match_type = 'Map_with_splicing_and_clipping'
    elif bool(re.search(r'\A\d+S\d+M\d+N\d+M\Z', cigar_str)):
        match_type = 'Map_with_splicing_and_clipping'
    else:
        match_type = 'Others'
    return match_type

def list2str (lst):
    '''
    translate samtools returned cigar_list into cigar_string
    '''
    code2Char={'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}
    cigar_str=''
    for i in lst:
        cigar_str += str(i[1]) + code2Char[str(i[0])]
    return cigar_str

def barcode_edits(infile, outfile, step_size=10000, limit=2000000, CR_tag = 'CR', CB_tag = 'CB', UR_tag = 'UR', UB_tag = 'UB'):
    '''
    Analzye barcode in BAM file.

    Parameters
    ----------
    infile : str
        Input BAM file. Must be sorted and indexed.
    outfile : str
        Prefix of output files.
    step_size: int
        Output progress report when step_size alignments have been processed.
    limit : int
        Only process this number of alignments and stop.
    '''
    logging.info("Reading BAM file \"%s\" ..." % infile)
    samfile = pysam.AlignmentFile(infile, 'rb')

    CB_miss = 0 #number of reads without cell barcode
    CB_same = 0 #number of reads whose original cell barcode same as edited barcode
    CB_diff = 0 #number of reads whose cell barcode has been edited
    CB_freq = collections.defaultdict(int) # cell barcode: raw reads
    CB_corrected_bases = collections.defaultdict(collections.Counter)


    UMI_miss = 0
    UMI_same = 0
    UMI_diff = 0
    UMI_freq = collections.defaultdict(int) # UMI : raw reads
    UMI_corrected_bases = collections.defaultdict(collections.Counter)

    total_alignments = 0
    try:
        while True:
            aligned_read = next(samfile)
            if aligned_read.is_secondary or aligned_read.is_supplementary:
                continue    # count each read once
            total_alignments += 1
            tag_dict = dict(aligned_read.tags) #{'NM': 1, 'RG': 'L1'}

            original_CB = ''
            corrected_CB = ''
            if CR_tag in tag_dict and CB_tag in tag_dict:
                original_CB = strip_suffix(tag_dict[CR_tag])
                corrected_CB = strip_suffix(tag_dict[CB_tag])
                CB_freq[corrected_CB] +=1
                if original_CB != corrected_CB:
                    CB_diff += 1
                    for diff in diff_str(original_CB, corrected_CB):
                        CB_corrected_bases[diff[0]][diff[1] + ":" + diff[2]] += 1
                else:
                    CB_same +=1
            else:
                CB_miss += 1

            original_UMI = ''
            corrected_UMI = ''
            if UR_tag in tag_dict and UB_tag in tag_dict:
                original_UMI = strip_suffix(tag_dict[UR_tag])
                corrected_UMI = strip_suffix(tag_dict[UB_tag])
                UMI_freq[corrected_UMI] += 1
                if original_UMI != corrected_UMI:
                    UMI_diff += 1
                    for diff in diff_str(original_UMI, corrected_UMI):
                        UMI_corrected_bases[diff[0]][diff[1] + ":" + diff[2]] += 1
                else:
                    UMI_same += 1
            else:
                UMI_miss += 1

            if total_alignments % step_size == 0:
                print("%d alignments processed.\r" % total_alignments, end=' ', file=sys.stderr)
            if limit is not None:
                if total_alignments >= limit:
                    break

    except StopIteration:
        pass
    logging.info ("Total alignments processed: %d" % total_alignments)

    logging.info("Number of alignmenets with <cell barcode> kept AS IS: %d" % CB_same)
    logging.info("Number of alignmenets with <cell barcode> edited: %d" % CB_diff)
    logging.info("Number of alignmenets with <cell barcode> missing: %d" % CB_miss)
    logging.info("Number of alignmenets with UMI kept AS IS: %d" % UMI_same)
    logging.info("Number of alignmenets with UMI edited: %d" % UMI_diff)
    logging.info("Number of alignmenets with UMI missing: %d" % UMI_miss)

    # writing cell barcode
    logging.info ("Writing cell barcode frequencies to \"%s\"" % (outfile + '.CB_freq.tsv'))
    with open(outfile + '.CB_freq.tsv','w') as CB_OUT:
        for bc,count in sorted(CB_freq.items(), key=lambda item: item[1], reverse=True):
            CB_OUT.write(bc + '\t' + str(count) + '\n')

    # writing UMI
    logging.info ("Writing UMI frequencies to \"%s\"" % (outfile + '.UMI_freq.tsv'))
    with open(outfile + '.UMI_freq.tsv','w') as UMI_OUT:
        for bc,count in sorted(UMI_freq.items(), key=lambda item: item[1], reverse=True):
            UMI_OUT.write(bc + '\t' + str(count) + '\n')

    CB_mat_file = outfile + '.CB_edits_count.csv'
    logging.info ("Writing the nucleotide editing matrix (count) of cell barcode to \"%s\"" % CB_mat_file)
    CB_diff_mat = pd.DataFrame.from_dict(CB_corrected_bases)
    CB_diff_mat = CB_diff_mat.fillna(0)
    #CB_diff_mat = CB_diff_mat.T
    CB_diff_mat.sort_index(inplace=True)
    CB_diff_mat = CB_diff_mat.sort_index(axis=1)
    CB_diff_mat.index.name='Edits'
    CB_diff_mat.to_csv(CB_mat_file, index=True, index_label="Index")

    UMI_mat_file = outfile + '.UMI_edits_count.csv'
    logging.info ("Writing the nucleotide editing matrix of molecular barcode (UMI) to \"%s\"" % UMI_mat_file)
    UMI_diff_mat = pd.DataFrame.from_dict(UMI_corrected_bases)
    UMI_diff_mat = UMI_diff_mat.fillna(0)
    #UMI_diff_mat = UMI_diff_mat.T
    UMI_diff_mat.sort_index(inplace=True)
    UMI_diff_mat = UMI_diff_mat.sort_index(axis=1)
    UMI_diff_mat.index.name='Edits'
    UMI_diff_mat.to_csv(UMI_mat_file, index=True, index_label="Index")



def mapping_stat(infile,  step_size=50000, CB_tag = 'CB', UMI_tag = 'UB',RE_tag = 'RE', TX_tag = 'TX', AN_tag = 'AN', xf_tag = 'xf', chrM_id ='chrM', n_thread = 1):
    '''
    Reads mapping statistics

    Parameters
    ----------
    infile : str
        Input BAM file. Must be sorted and indexed.
    outfile : str
        Prefix of output files. If outfile is None, only do counting and do not generate BAM files.
    step_size: int
        Output progress report when step_size alignments have been processed.
    limit : int
        Only process this number of alignments and stop.
    '''
    logging.info("Reading BAM file \"%s\" ..." % infile)
    try:
        samfile = pysam.AlignmentFile(
            infile, mode="rb", require_index=True, threads=n_thread
        )
    except TypeError:
        # Compatibility with older pysam releases.
        samfile = pysam.AlignmentFile(
            infile, mode="rb", require_index=True, thread=n_thread
        )
    if not samfile.check_index():
        samfile.close()
        raise FileNotFoundError(f"Cannot find BAM index for {infile}")
    chrom_info = zip(samfile.references, samfile.lengths) #[('chr1', 195471971), ('chr10', 130694993),...]

    total_alignments = 0
    confi_alignments = 0

    total_reads_n = 0
    confi_reads_n = 0

    ## for confidently mapped reads
    #PCR duplicate or not
    confi_reads_dup = 0
    confi_reads_nondup = 0

    #Reverse or forward
    confi_reads_rev = 0
    confi_reads_fwd = 0

    # with error-corrected CB or UMI
    confi_CB = 0
    confi_UB = 0

    exon_reads = 0
    intron_reads = 0
    intergenic_reads = 0
    other_reads1 = 0

    #sense or antisense
    sense_reads = 0
    anti_reads = 0
    other_reads2 =0

    chrM_reads = 0
    #read match type
    read_type = collections.defaultdict(int)
    # Read names are collected in memory. (Before RSeQC 5.1.0+lbzhao they went to
    # per-chromosome temporary files in the working directory, and every
    # matching *.all_reads_id.txt file there was then read and deleted, so
    # two runs in one directory corrupted each other.)
    all_read_ids = set()
    confident_read_ids = set()

    for chr_id, chr_len in chrom_info:
        logging.info("Processing \"%s\" ..." % chr_id)
        chrom_count = 0
        chrom_total_reads = all_read_ids          # total reads in BAM file
        chrom_confi_reads = confident_read_ids    # confidently mapped (xf:i:1)

        for aligned_read in samfile.fetch(chr_id):
            total_alignments += 1
            chrom_count += 1
            read_id = aligned_read.query_name
            tag_dict = dict(aligned_read.tags) #{'NM': 1, 'RG': 'L1'}
            cigar_str = list2str( aligned_read.cigar)
            chrom_total_reads.add(read_id)

            #confident alignments
            if xf_tag in tag_dict and tag_dict[xf_tag]& 0x1 != 0:

                if chr_id == chrM_id:
                    chrM_reads += 1
                #with or without CB/UMI barcode
                if CB_tag in tag_dict:
                    confi_CB += 1
                if UMI_tag in tag_dict:
                    confi_UB += 1

                #duplicate or not
                if aligned_read.is_duplicate:
                    confi_reads_dup += 1
                else:
                    confi_reads_nondup += 1

                #forward or reverse
                if aligned_read.is_reverse:
                    confi_reads_rev += 1
                else:
                    confi_reads_fwd += 1

                #Single character indicating the region type of this alignment (E = exonic, N = intronic, I = intergenic).
                if RE_tag in tag_dict:
                    if tag_dict[RE_tag] == "E":
                        exon_reads += 1
                    elif tag_dict[RE_tag] == "N":
                        intron_reads += 1
                    elif tag_dict[RE_tag] == "I":
                        intergenic_reads += 1
                    else:
                        other_reads1 += 1
                else:
                    other_reads1 += 1

                #sense or antisense
                if TX_tag in tag_dict:
                    sense_reads += 1
                elif AN_tag in tag_dict:
                    anti_reads += 1
                else:
                    other_reads2 += 1

                #map type
                cigar_str = list2str( aligned_read.cigar)
                tmp = read_match_type(cigar_str)
                read_type[tmp] += 1

                confi_alignments += 1
                chrom_confi_reads.add(read_id)

            if chrom_count % step_size == 0:
                print("%d alignments processed.\r" % chrom_count, end=' ', file=sys.stderr)
        logging.info("Processed %d alignments from \"%s\"" % (chrom_count, chr_id))

    logging.info("Processing total %d alignments mapped to all chromosomes." % total_alignments)

    logging.info("Count total mapped reads ...")
    with open("All_reads_uniqID.txt", "w") as handle:
        for read_id in sorted(all_read_ids):
            print(read_id, file=handle)
    with open("confident_reads_uniqID.txt", "w") as handle:
        for read_id in sorted(confident_read_ids):
            print(read_id, file=handle)

    total_reads_n = len(all_read_ids)
    confi_reads_n = len(confident_read_ids)
    non_confi_reads = total_reads_n - confi_reads_n



    print ('')
    print("\nTotal_alignments: %d" % total_alignments)
    print ("└--Confident_alignments: %d" % confi_alignments)
    print ('')
    print ("Total_mapped_reads:\t%d" % total_reads_n)
    print ("|--Non_confidently_mapped_reads:\t%d\t(%.2f%%)" % (non_confi_reads, _percent(non_confi_reads, total_reads_n)))
    print ("└--Confidently_mapped_reads:\t%d\t(%.2f%%)" % (confi_reads_n, _percent(confi_reads_n, total_reads_n)))

    print ("   |--Reads_with_PCR_duplicates:\t%d\t(%.2f%%)" % (confi_reads_dup, _percent(confi_reads_dup, confi_reads_n)))
    print ("   └--Reads_no_PCR_duplicates:\t%d\t(%.2f%%)" % (confi_reads_nondup, _percent(confi_reads_nondup, confi_reads_n)))
    print ('')

    print ("   |--Reads_map_to_forward(Waston)_strand:\t%d\t(%.2f%%)" % (confi_reads_fwd, _percent(confi_reads_fwd, confi_reads_n)))
    print ("   └--Reads_map_to_Reverse(Crick)_strand:\t%d\t(%.2f%%)" % (confi_reads_rev, _percent(confi_reads_rev, confi_reads_n)))
    print ('')

    print ("   |--Reads_map_to_sense_strand:\t%d\t(%.2f%%)" % (sense_reads, _percent(sense_reads, confi_reads_n)))
    print ("   └--Reads_map_to_antisense_strand:\t%d\t(%.2f%%)" % (anti_reads, _percent(anti_reads, confi_reads_n)))
    print ("   └--Other:\t%d\t(%.2f%%)" % (other_reads2, _percent(other_reads2, confi_reads_n)))
    print ('')

    print ("   |--Reads_map_to_exons:\t%d\t(%.2f%%)" % (exon_reads, _percent(exon_reads, confi_reads_n)))
    print ("   └--Reads_map_to_introns:\t%d\t(%.2f%%)" % (intron_reads, _percent(intron_reads, confi_reads_n)))
    print ("   └--Reads_map_to_intergenic:\t%d\t(%.2f%%)" % (intergenic_reads, _percent(intergenic_reads, confi_reads_n)))
    print ("   └--Other:\t%d\t(%.2f%%)" % (other_reads1, _percent(other_reads1, confi_reads_n)))
    print ('')

    print ("   |--Reads_with_error-corrected_barcode:\t%d\t(%.2f%%)" % (confi_CB, _percent(confi_CB, confi_reads_n)))
    print ("   └--Reads_no_error-corrected_barcode:\t%d\t(%.2f%%)" % ((confi_reads_n - confi_CB), _percent(confi_reads_n - confi_CB, confi_reads_n)))
    print ('')
    print ("   |--Reads_with_error-corrected_UMI:\t%d\t(%.2f%%)" % (confi_UB, _percent(confi_UB, confi_reads_n)))
    print ("   └--Reads_no_error-corrected_UMI:\t%d\t(%.2f%%)" % ((confi_reads_n - confi_UB), _percent(confi_reads_n - confi_UB, confi_reads_n)))
    print ('')
    print ("   |--Reads_map_to_mitochonrial_genome:\t%d\t(%.2f%%)" % (chrM_reads, _percent(chrM_reads, confi_reads_n)))
    print ("   └--Reads_map_to_nuclear_genome:\t%d\t(%.2f%%)" % ((confi_reads_n - chrM_reads), _percent(confi_reads_n - chrM_reads, confi_reads_n)))

    print ('')

    for i in sorted(read_type):
        if i == 'Others':continue
        print ("   |--%s:\t%d\t(%.2f%%)" % (i, read_type[i], _percent(read_type[i], confi_reads_n)))
    print ("   └--%s:\t%d\t(%.2f%%)" % ('Others', read_type['Others'], _percent(read_type['Others'], confi_reads_n)))
    print ('')



