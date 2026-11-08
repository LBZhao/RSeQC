"""Tests for qcmodule.bam_cigar (CIGAR -> coordinates)."""

from __future__ import annotations

import random

import pysam
import pytest

from qcmodule import bam_cigar

CHROM = "chr1"
START = 1000

HEADER = pysam.AlignmentHeader.from_dict(
    {"SQ": [{"SN": CHROM, "LN": 1_000_000}]}
)


def make_read(cigar, start=START):
    read = pysam.AlignedSegment(HEADER)
    read.reference_id = 0
    read.reference_start = start
    read.cigarstring = cigar
    return read


def spans(blocks):
    return [(b[1], b[2]) for b in blocks]


def covered(blocks):
    return [pos for _, start, end in blocks for pos in range(start, end)]


# --------------------------------------------------------------------------
# Hand-picked cases
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cigar, exons, introns",
    [
        ("100M", [(1000, 1100)], []),
        # Soft clips are not on the reference and must not move blocks.
        ("10S90M", [(1000, 1090)], []),
        ("90M10S", [(1000, 1090)], []),
        ("5S45M1000N50M", [(1000, 1045), (2045, 2095)], [(1045, 2045)]),
        ("3H5S45M1000N50M2S4H", [(1000, 1045), (2045, 2095)], [(1045, 2045)]),
        # =/X are aligned bases, exactly like M.
        ("50=1X49=", [(1000, 1100)], []),
        ("30=500N20X", [(1000, 1030), (1530, 1550)], [(1030, 1530)]),
        ("10M5=500N10X10M", [(1000, 1015), (1515, 1535)], [(1015, 1515)]),
        # Deletions and insertions split blocks (as in earlier releases);
        # deleted bases are not covered.
        ("20M10D20M", [(1000, 1020), (1030, 1050)], []),
        ("20M2I20M", [(1000, 1020), (1020, 1040)], []),
        ("20M10D5S", [(1000, 1020)], []),
    ],
)
def test_exons_and_introns(cigar, exons, introns):
    read = make_read(cigar)
    assert spans(bam_cigar.fetch_exon(CHROM, START, read.cigartuples)) == exons
    assert spans(bam_cigar.fetch_intron(CHROM, START, read.cigartuples)) == introns


def test_fetch_exon_keeps_chrom_and_tuple_shape():
    blocks = bam_cigar.fetch_exon("chrX", 5, [(0, 10)])
    assert blocks == [("chrX", 5, 15)]


def test_map_bounds_ignores_clips_and_counts_gaps():
    read = make_read("5S10M2D10M100N10M3S")
    assert bam_cigar.map_bounds(START, read.cigartuples) == (
        read.reference_start, read.reference_end,
    )


def test_deletion_and_insertion():
    read = make_read("5S10M3D10M2I10=1X")
    cigar = read.cigartuples
    assert spans(bam_cigar.fetch_deletion(CHROM, START, cigar)) == [(1010, 1013)]
    # Read offsets count soft-clipped bases from the start of SEQ.
    assert bam_cigar.fetch_deletion_range(cigar) == [(15, 3)]
    assert bam_cigar.fetch_insertion_range(cigar) == [(25, 2)]
    # Insertion sits before reference base 1023.
    assert bam_cigar.fetch_insertion(CHROM, START, cigar) == [(CHROM, 1023, 2)]
    # =/X advance read offsets too.
    assert bam_cigar.fetch_insertion_range(make_read("10=1X2I5M").cigartuples) \
        == [(11, 2)]


def test_fetch_clip_positions():
    cigar = make_read("5S90M7S").cigartuples
    assert spans(bam_cigar.fetch_clip(CHROM, START, cigar)) == [
        (995, 1000), (1090, 1097),
    ]


# --------------------------------------------------------------------------
# Random CIGARs checked against pysam
# --------------------------------------------------------------------------


def random_cigar(rng):
    middle_ops = "MMMM==XIDN"
    ops = []
    if rng.random() < 0.3:
        ops.append(f"{rng.randint(1, 5)}H")
    if rng.random() < 0.5:
        ops.append(f"{rng.randint(1, 20)}S")
    ops.append(f"{rng.randint(1, 30)}{rng.choice('M=X')}")
    for _ in range(rng.randint(0, 6)):
        op = rng.choice(middle_ops)
        size = rng.randint(50, 5000) if op == "N" else rng.randint(1, 30)
        ops.append(f"{size}{op}")
        if op in "IDN":
            ops.append(f"{rng.randint(1, 30)}{rng.choice('M=X')}")
    if rng.random() < 0.5:
        ops.append(f"{rng.randint(1, 20)}S")
    if rng.random() < 0.3:
        ops.append(f"{rng.randint(1, 5)}H")
    return "".join(ops)


def test_random_cigars_match_pysam():
    rng = random.Random(20261108)
    for _ in range(2000):
        cigar = random_cigar(rng)
        read = make_read(cigar)
        tuples = read.cigartuples

        blocks = bam_cigar.fetch_exon(CHROM, START, tuples)
        assert covered(blocks) == read.get_reference_positions(), cigar
        assert bam_cigar.map_bounds(START, tuples) == (
            read.reference_start, read.reference_end,
        ), cigar

        # Introns are exactly the N gaps between aligned blocks.
        expected_introns = []
        pos = START
        for code, size in tuples:
            if code == 3:
                expected_introns.append((pos, pos + size))
            if code in (0, 2, 3, 7, 8):
                pos += size
        assert spans(bam_cigar.fetch_intron(CHROM, START, tuples)) == \
            expected_introns, cigar

        deleted = {
            ref for q, ref in read.get_aligned_pairs() if q is None and ref is not None
        } - {p for s, e in expected_introns for p in range(s, e)}
        assert set(covered(bam_cigar.fetch_deletion(CHROM, START, tuples))) == \
            deleted, cigar
