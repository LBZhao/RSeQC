#!/usr/bin/env python

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Convert BAM CIGAR tuples into genomic and read coordinates.

CIGAR operations (code, letter, what it consumes)::

    0  M  alignment match      reference + read
    1  I  insertion            read
    2  D  deletion             reference
    3  N  skipped (intron)     reference
    4  S  soft clip            read
    5  H  hard clip            -
    6  P  padding              -
    7  =  sequence match       reference + read
    8  X  sequence mismatch    reference + read

The start position of a BAM record (``reference_start``, 0-based) is the
first *aligned* reference base: soft-clipped bases are not on the reference
and do not move the genomic position. ``=`` and ``X`` are aligned bases just
like ``M``.

Genomic intervals are 0-based, half-open. Before RSeQC 5.1.0+lbzhao soft clips
shifted every downstream block and ``=``/``X`` operations were ignored.
"""

from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"

MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, PAD, EQUAL, DIFF = range(9)

#: Operations that are aligned bases (consume both reference and read).
ALIGNED_OPS = frozenset((MATCH, EQUAL, DIFF))
#: Operations that move along the reference.
REFERENCE_OPS = frozenset((MATCH, DEL, REF_SKIP, EQUAL, DIFF))
#: Operations that move along the read sequence (SEQ).
QUERY_OPS = frozenset((MATCH, INS, SOFT_CLIP, EQUAL, DIFF))


def map_bounds(start_pos, cigar):
    """Return ``(start, end)`` of the reference span covered by a read.

    Includes deletions and introns, excludes soft clips.
    """
    span = sum(size for code, size in cigar if code in REFERENCE_OPS)
    return (start_pos, start_pos + span)


def fetch_exon(chrom, st, cigar):
    """Return aligned blocks ``[(chrom, start, end), ...]`` on the reference.

    *st* is the 0-based alignment start. Consecutive ``M``/``=``/``X``
    operations form one block; blocks are split by insertions, deletions and
    introns. Soft-clipped and deleted bases are not part of any block.
    """
    pos = st
    blocks = []
    block_start = None
    for code, size in cigar:
        if code in ALIGNED_OPS:
            if block_start is None:
                block_start = pos
            pos += size
            continue
        if block_start is not None and code in (INS, DEL, REF_SKIP):
            blocks.append((chrom, block_start, pos))
            block_start = None
        if code in REFERENCE_OPS:
            pos += size
    if block_start is not None:
        blocks.append((chrom, block_start, pos))
    return blocks


def fetch_intron(chrom, st, cigar):
    """Return introns (``N`` operations) as ``[(chrom, start, end), ...]``."""
    pos = st
    introns = []
    for code, size in cigar:
        if code == REF_SKIP:
            introns.append((chrom, pos, pos + size))
        if code in REFERENCE_OPS:
            pos += size
    return introns


def fetch_clip(chrom, st, cigar):
    """Return soft-clipped regions as ``[(chrom, start, end), ...]``.

    Soft-clipped bases are not aligned, so the interval is where they would
    sit if extended along the reference: a leading clip ends at *st* and a
    trailing clip starts at the alignment end.
    """
    pos = st
    clips = []
    for code, size in cigar:
        if code == SOFT_CLIP:
            if pos == st:       # nothing aligned yet: leading clip
                clips.append((chrom, st - size, st))
            else:
                clips.append((chrom, pos, pos + size))
        elif code in REFERENCE_OPS:
            pos += size
    return clips


def fetch_deletion(chrom, st, cigar):
    """Return deletions (``D`` operations) as ``[(chrom, start, end), ...]``."""
    pos = st
    deletions = []
    for code, size in cigar:
        if code == DEL:
            deletions.append((chrom, pos, pos + size))
        if code in REFERENCE_OPS:
            pos += size
    return deletions


def fetch_deletion_range(cigar):
    """Return deletions as ``[(read_offset, size), ...]``.

    *read_offset* counts read bases (including soft clips) before the
    deletion.
    """
    offset = 0
    deletions = []
    for code, size in cigar:
        if code == DEL:
            deletions.append((offset, size))
        if code in QUERY_OPS:
            offset += size
    return deletions


def fetch_insertion_range(cigar):
    """Return insertions as ``[(read_offset, size), ...]``.

    *read_offset* counts read bases (including soft clips) before the
    insertion.
    """
    offset = 0
    insertions = []
    for code, size in cigar:
        if code == INS:
            insertions.append((offset, size))
        if code in QUERY_OPS:
            offset += size
    return insertions


def fetch_insertion(chrom, st, cigar):
    """Return insertions as ``[(chrom, position, size), ...]``.

    Inserted bases have no reference coordinate; ``(chrom, 100, 2)`` means
    two bases inserted between reference positions 99 and 100 (0-based),
    i.e. before reference base 100.
    """
    pos = st
    insertions = []
    for code, size in cigar:
        if code == INS:
            insertions.append((chrom, pos, size))
        if code in REFERENCE_OPS:
            pos += size
    return insertions


def list2str (lst):
    '''translate samtools returned cigar_list into cigar_string
    
    [(4, 1), (0, 9)] ==> '1S9M'
    '''
    code2Char={'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}
        
    cigar_str=''
    for i in lst:
        cigar_str += str(i[1]) + code2Char[str(i[0])]
    return cigar_str

def list2longstr (lst):
    '''translate samtools returned cigar_list into LONG cigar_string
    Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
    '''
    
    code2Char={'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}
        
    cigar_str=''
    for i in lst:
        (code,length) = (i[0],i[1])
        if code in (0,1,4,7,8):
            cigar_str += code2Char[str(code)] * int(length)
        else:
            continue
    return cigar_str

    
