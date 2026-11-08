"""Helpers for building small synthetic BAM and BED12 fixtures."""

from __future__ import annotations

import sys
from pathlib import Path

import pysam
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CHROM = "chr1"
CHROM_LENGTH = 200_000


def bed12_line(name, strand, blocks, chrom=CHROM):
    """Return a BED12 line for genomic exon *blocks* [(start, end), ...]."""
    start = blocks[0][0]
    end = blocks[-1][1]
    sizes = ",".join(str(e - s) for s, e in blocks)
    offsets = ",".join(str(s - start) for s, _ in blocks)
    return (
        f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t0\t"
        f"{len(blocks)}\t{sizes},\t{offsets},\n"
    )


def exonic_alignment(blocks, offset, length):
    """Genomic start and CIGAR for a read covering exonic bases
    ``[offset, offset + length)`` counted from the genomic left end."""
    ops = []
    ref_start = None
    remaining = length
    position = 0
    previous_end = None
    for block_start, block_end in blocks:
        block_len = block_end - block_start
        if offset >= position + block_len:
            position += block_len
            continue
        local = max(0, offset - position)
        take = min(block_len - local, remaining)
        genomic = block_start + local
        if ref_start is None:
            ref_start = genomic
        elif genomic > previous_end:
            ops.append(f"{genomic - previous_end}N")
        ops.append(f"{take}M")
        previous_end = genomic + take
        remaining -= take
        position += block_len
        if remaining == 0:
            break
    if remaining:
        raise ValueError("read runs past the transcript end")
    return ref_start, "".join(ops)


class BamBuilder:
    """Collect alignments and write a sorted, indexed BAM file."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.header = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": CHROM_LENGTH}]}
        )
        self.reads = []
        self._serial = 0

    def _name(self):
        self._serial += 1
        return f"r{self._serial}"

    def add(self, ref_start, cigar, *, flag=0, mapq=60, name=None, mate=None):
        read = pysam.AlignedSegment(self.header)
        read.query_name = name or self._name()
        read.flag = flag
        read.reference_id = 0
        read.reference_start = ref_start
        read.mapping_quality = mapq
        read.cigarstring = cigar
        length = read.query_length or sum(
            n for op, n in read.cigartuples if op in (0, 1, 4, 7, 8)
        )
        read.query_sequence = "A" * length
        read.query_qualities = pysam.qualitystring_to_array("I" * length)
        if mate is not None:
            read.next_reference_id = 0
            read.next_reference_start = mate
        self.reads.append(read)
        return read

    def add_exonic(self, blocks, offset, length, **kwargs):
        start, cigar = exonic_alignment(blocks, offset, length)
        return self.add(start, cigar, **kwargs)

    def add_pair(self, start1, cigar1, start2, cigar2, *, flag_extra=0):
        name = self._name()
        self.add(start1, cigar1, flag=0x1 | 0x2 | 0x20 | 0x40 | flag_extra,
                 name=name, mate=start2)
        self.add(start2, cigar2, flag=0x1 | 0x2 | 0x10 | 0x80 | flag_extra,
                 name=name, mate=start1)

    def tile(self, blocks, depth, read_length=50, first=0, last=None):
        """Tile exonic bases [first, last) with *depth* non-overlapping reads."""
        total = sum(e - s for s, e in blocks)
        last = total if last is None else last
        for _ in range(depth):
            for offset in range(first, last, read_length):
                self.add_exonic(blocks, offset, min(read_length, last - offset))

    def write(self) -> Path:
        unsorted = self.path.with_suffix(".unsorted.bam")
        with pysam.AlignmentFile(str(unsorted), "wb", header=self.header) as out:
            for read in self.reads:
                out.write(read)
        pysam.sort("-o", str(self.path), str(unsorted))
        pysam.index(str(self.path))
        unsorted.unlink()
        return self.path


@pytest.fixture
def bam_builder(tmp_path):
    def make(name="sample.bam"):
        return BamBuilder(tmp_path / name)
    return make
