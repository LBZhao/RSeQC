#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate fragment counts, FPM, and FPKM for BED12 transcript models.

The counting rules intentionally match the original RSeQC implementation.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Sequence
from bx.intervals.intersection import Intersecter, Interval
from qcmodule import SAM
from qcmodule import readfilter
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"

def build_parser() -> argparse.ArgumentParser:
    """Create the command-line parser without changing legacy semantics."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=Path,
        help="Alignment file in BAM format. SAM is not supported.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        type=Path,
        help="Prefix for output files.",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="refgene_bed",
        required=True,
        type=Path,
        help="Reference gene model in BED12 format.",
    )
    parser.add_argument(
        "-d",
        "--strand",
        dest="strand_rule",
        default=None,
        help=(
            "Strand rule, for example '1++,1--,2+-,2-+'. "
            "Omit for unstranded RNA-seq."
        ),
    )
    parser.add_argument(
        "-u",
        "--skip-multi-hits",
        action="store_true",
        dest="skip_multi",
        help="Skip alignments with mapping quality below --mapq.",
    )
    parser.add_argument(
        "-e",
        "--only-exonic",
        action="store_true",
        dest="only_exon",
        help="Use exonic fragments rather than all fragments for normalization.",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        dest="map_qual",
        type=int,
        default=30,
        metavar="INT",
        help="Minimum mapping quality. Default: %(default)s",
    )
    parser.add_argument(
        "-s",
        "--single-read",
        dest="single_read",
        type=float,
        default=1.0,
        metavar="FLOAT",
        help=(
            "Weight for pairs with only one mapped end: 0, 0.5, or 1. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def find_bam_index(bam_path: Path) -> Path | None:
    """Return an existing BAM index path."""
    candidates = (
        Path(f"{bam_path}.bai"),
        bam_path.with_suffix(".bai"),
    )
    return next((path for path in candidates if path.is_file()), None)


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate files and numeric options without altering counting behavior."""
    if not args.input_file.is_file():
        parser.error(f"input BAM file does not exist: {args.input_file}")

    if args.input_file.suffix.lower() != ".bam":
        parser.error(f"input file must be BAM: {args.input_file}")

    if find_bam_index(args.input_file) is None:
        parser.error(
            "cannot find BAM index; expected either "
            f"{args.input_file}.bai or {args.input_file.with_suffix('.bai')}"
        )

    if not args.refgene_bed.is_file():
        parser.error(f"BED12 annotation does not exist: {args.refgene_bed}")

    if args.map_qual < 0:
        parser.error("--mapq must be zero or greater")

    if args.single_read < 0:
        parser.error("--single-read must be zero or greater")

    if not args.out_prefix.parent.exists():
        parser.error(f"output directory does not exist: {args.out_prefix.parent}")


def parse_strand_rule(strand_rule: str | None) -> dict[str, str]:
    """Parse the RSeQC strand-rule syntax (see qcmodule.readfilter)."""
    return readfilter.parse_strand_rule(strand_rule)


def build_exon_ranges(refgene_bed: Path) -> dict[str, Intersecter]:
    """Build the chromosome-level exon index using the original algorithm."""
    ranges: dict[str, Intersecter] = {}

    with refgene_bed.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                if line.startswith(("#", "track", "browser")):
                    continue

                fields = line.split()
                chrom = fields[0].upper()
                tx_start = int(fields[1])

                exon_starts = [
                    int(value)
                    for value in fields[11].rstrip(",\n").split(",")
                ]
                exon_starts = [
                    exon_start + tx_start
                    for exon_start in exon_starts
                ]

                exon_sizes = [
                    int(value)
                    for value in fields[10].rstrip(",\n").split(",")
                ]
                exon_ends = [
                    exon_start + exon_size
                    for exon_start, exon_size in zip(exon_starts, exon_sizes)
                ]

            except (IndexError, ValueError):
                print(
                    "[NOTE: input BED must be 12-column] skipped this line: "
                    + line,
                    end=" ",
                    file=sys.stderr,
                )
                continue

            chromosome_ranges = ranges.setdefault(chrom, Intersecter())
            for start, end in zip(exon_starts, exon_ends):
                chromosome_ranges.add_interval(Interval(start, end))

    return ranges


def alignment_passes_legacy_filters(
    aligned_read,
    *,
    skip_multi: bool,
    map_qual: int,
) -> bool:
    """Apply exactly the filters used by the original script."""
    if aligned_read.is_qcfail:
        return False

    if aligned_read.is_duplicate:
        return False

    if readfilter.is_non_primary(aligned_read):
        return False

    if skip_multi and aligned_read.mapq < map_qual:
        return False

    return True


def count_total_fragments(
    bam_path: Path,
    gene_ranges: dict[str, Intersecter],
    *,
    skip_multi: bool,
    map_qual: int,
    single_read: float,
) -> tuple[float, float]:
    """Count total and exonic fragments using the original logic exactly."""
    obj = SAM.ParseBAM(str(bam_path))
    total_frags = 0.0
    exonic_frags = 0.0

    print("Counting total fragment ... ", end=" ", file=sys.stderr)

    for aligned_read in obj.samfile:
        if not alignment_passes_legacy_filters(
            aligned_read,
            skip_multi=skip_multi,
            map_qual=map_qual,
        ):
            continue

        try:
            chrom = obj.samfile.getrname(aligned_read.tid).upper()
        except (ValueError, IndexError):
            continue

        read_st = aligned_read.reference_start
        # Reference end includes introns, so spliced reads reach their exons.
        read_end = (
            aligned_read.reference_end
            if aligned_read.reference_end is not None
            else read_st + (aligned_read.query_length or 1)
        )

        if not aligned_read.is_paired:
            total_frags += 1

            if (
                chrom in gene_ranges
                and len(gene_ranges[chrom].find(read_st, read_end)) > 0
            ):
                exonic_frags += 1

            continue

        if aligned_read.is_read2:
            continue

        mate_st = aligned_read.pnext
        mate_end = mate_st + aligned_read.rlen

        if aligned_read.is_unmapped:
            if aligned_read.mate_is_unmapped:
                continue

            total_frags += single_read

            if (
                chrom in gene_ranges
                and len(gene_ranges[chrom].find(mate_st, mate_end)) > 0
            ):
                exonic_frags += single_read

        elif aligned_read.mate_is_unmapped:
            total_frags += single_read

            if (
                chrom in gene_ranges
                and len(gene_ranges[chrom].find(read_st, read_end)) > 0
            ):
                exonic_frags += single_read

        else:
            total_frags += 1

            if (
                chrom in gene_ranges
                and len(gene_ranges[chrom].find(read_st, read_end)) > 0
                and len(gene_ranges[chrom].find(mate_st, mate_end)) > 0
            ):
                exonic_frags += 1

    print("Done", file=sys.stderr)
    return total_frags, exonic_frags


def parse_transcript_line(line: str) -> tuple[
    str,
    int,
    int,
    str,
    str,
    list[int],
    list[int],
]:
    """Parse one BED12 line using the original coordinate calculations."""
    fields = line.split()

    chrom = fields[0]
    tx_start = int(fields[1])
    tx_end = int(fields[2])
    gene_name = fields[3]
    gene_strand = fields[5].replace(" ", "_")

    exon_starts = [
        int(value)
        for value in fields[11].rstrip(",\n").split(",")
    ]
    exon_starts = [
        exon_start + tx_start
        for exon_start in exon_starts
    ]

    exon_sizes = [
        int(value)
        for value in fields[10].rstrip(",\n").split(",")
    ]
    exon_ends = [
        exon_start + exon_size
        for exon_start, exon_size in zip(exon_starts, exon_sizes)
    ]

    return (
        chrom,
        tx_start,
        tx_end,
        gene_name,
        gene_strand,
        exon_starts,
        exon_ends,
    )


def mate_reference_end(aligned_read) -> int:
    """Reference end of the mate: from its CIGAR (MC tag) when present,
    otherwise estimated from this read's length."""
    if aligned_read.has_tag("MC"):
        cigar = bam_cigar_tuples(aligned_read.get_tag("MC"))
        return aligned_read.next_reference_start + sum(
            size for op, size in cigar if op in (0, 2, 3, 7, 8)
        )
    return aligned_read.next_reference_start + (aligned_read.query_length or 1)


def bam_cigar_tuples(cigar: str) -> list[tuple[int, int]]:
    codes = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
    return [(codes[op], int(size)) for size, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def overlaps_exons(aligned_read, exon_ranges: Intersecter) -> bool:
    """True if any aligned block of the read overlaps an exon."""
    return any(
        exon_ranges.find(start, end) for start, end in aligned_read.get_blocks()
    )


def count_transcript(
    samfile,
    *,
    chrom: str,
    tx_start: int,
    tx_end: int,
    exon_ranges: Intersecter,
    strand_rule: str | None,
    strand_map: dict[str, str],
    skip_multi: bool,
    map_qual: int,
    single_read: float,
) -> tuple[float, float, float]:
    """Count fragments of one transcript: (forward, reverse, unstranded).

    Single-end reads count if an aligned block overlaps an exon. A pair
    counts if either mate starts in an exon, and is counted once: through
    read 1, or through read 2 when read 1 lies outside the transcript and
    is not fetched. A pair with one unmapped mate counts *single_read*.
    Fragments overlapping several transcripts count for each of them.
    """
    counts = {"+": 0.0, "-": 0.0, None: 0.0}

    try:
        aligned_reads = samfile.fetch(chrom, tx_start, tx_end)
    except (ValueError, OSError):
        return 0.0, 0.0, 0.0

    for aligned_read in aligned_reads:
        if not alignment_passes_legacy_filters(
            aligned_read,
            skip_multi=skip_multi,
            map_qual=map_qual,
        ):
            continue

        if strand_rule is None:
            bucket = None
        else:
            bucket = strand_map.get(readfilter.read_key(aligned_read))
            if bucket is None:
                continue

        if not aligned_read.is_paired:
            if not aligned_read.is_unmapped and overlaps_exons(aligned_read, exon_ranges):
                counts[bucket] += 1
            continue

        if aligned_read.is_unmapped and aligned_read.mate_is_unmapped:
            continue

        frag_st = aligned_read.reference_start
        mate_st = aligned_read.next_reference_start
        if (
            len(exon_ranges.find(frag_st, frag_st + 1)) < 1
            and len(exon_ranges.find(mate_st, mate_st + 1)) < 1
        ):
            continue

        if aligned_read.is_read2 and not (aligned_read.is_unmapped or aligned_read.mate_is_unmapped):
            mate_fetched = (
                aligned_read.next_reference_id == aligned_read.reference_id
                and mate_st < tx_end
                and mate_reference_end(aligned_read) > tx_start
            )
            if mate_fetched:
                continue        # counted through read 1
            # Read 1 is outside this transcript: count the pair here, with
            # the strand of read 1 (the mate of a read 2 on '+' is on '-').
            if strand_rule is not None:
                mate_key = "1" + ("-" if aligned_read.mate_is_reverse else "+")
                bucket = strand_map.get(mate_key)
                if bucket is None:
                    continue
        elif aligned_read.is_read2:
            # One mate unmapped: the mapped mate is fetched, and read 2
            # is only counted when it is that mapped mate.
            if aligned_read.is_unmapped:
                continue

        if aligned_read.is_unmapped or aligned_read.mate_is_unmapped:
            counts[bucket] += single_read
        else:
            counts[bucket] += 1

    return counts["+"], counts["-"], counts[None]


def main(argv: Sequence[str] | None = None) -> int:
    """Run fragment counting and FPKM calculation."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    try:
        strand_map = parse_strand_rule(args.strand_rule)

        print(
            f"Extract exon regions from {args.refgene_bed}...",
            file=sys.stderr,
        )
        gene_ranges = build_exon_ranges(args.refgene_bed)

        total_frags, exonic_frags = count_total_fragments(
            args.input_file,
            gene_ranges,
            skip_multi=args.skip_multi,
            map_qual=args.map_qual,
            single_read=args.single_read,
        )

        print(f"Total fragment = {total_frags:<20}", file=sys.stderr)
        print(f"Total exonic fragment = {exonic_frags:<20}", file=sys.stderr)

        if total_frags <= 0 or exonic_frags <= 0:
            raise ValueError("total and exonic fragment counts must be positive")

        denominator = exonic_frags if args.only_exon else total_frags

        output_path = Path(f"{args.out_prefix}.FPKM.xls")
        obj = SAM.ParseBAM(str(args.input_file))

        with (
            args.refgene_bed.open("r", encoding="utf-8") as bed_handle,
            output_path.open("w", encoding="utf-8") as output,
        ):
            print(
                "\t".join(
                    (
                        "#chrom",
                        "st",
                        "end",
                        "accession",
                        "mRNA_size",
                        "gene_strand",
                        "Frag_count",
                        "FPM",
                        "FPKM",
                    )
                ),
                file=output,
            )

            gene_finished = 0

            for line in bed_handle:
                if line.startswith(("#", "track", "browser")):
                    continue

                (
                    chrom,
                    tx_start,
                    tx_end,
                    gene_name,
                    gene_strand,
                    exon_starts,
                    exon_ends,
                ) = parse_transcript_line(line)

                mrna_size = 0.0
                exon_ranges = Intersecter()

                for start, end in zip(exon_starts, exon_ends):
                    mrna_size += end - start
                    exon_ranges.add_interval(Interval(start, end))

                frag_count_f, frag_count_r, frag_count_fr = count_transcript(
                    obj.samfile,
                    chrom=chrom,
                    tx_start=tx_start,
                    tx_end=tx_end,
                    exon_ranges=exon_ranges,
                    strand_rule=args.strand_rule,
                    strand_map=strand_map,
                    skip_multi=args.skip_multi,
                    map_qual=args.map_qual,
                    single_read=args.single_read,
                )

                fpm_fr = frag_count_fr * 1_000_000 / denominator
                fpm_f = frag_count_f * 1_000_000 / denominator
                fpm_r = frag_count_r * 1_000_000 / denominator

                fpkm_fr = (
                    frag_count_fr * 1_000_000_000
                    / (denominator * mrna_size)
                )
                fpkm_f = (
                    frag_count_f * 1_000_000_000
                    / (denominator * mrna_size)
                )
                fpkm_r = (
                    frag_count_r * 1_000_000_000
                    / (denominator * mrna_size)
                )

                if args.strand_rule is None:
                    values = (
                        chrom,
                        tx_start,
                        tx_end,
                        gene_name,
                        mrna_size,
                        gene_strand,
                        frag_count_fr,
                        fpm_fr,
                        fpkm_fr,
                    )
                    print("\t".join(str(value) for value in values), file=output)

                elif gene_strand == "+":
                    values = (
                        chrom,
                        tx_start,
                        tx_end,
                        gene_name,
                        mrna_size,
                        gene_strand,
                        frag_count_f,
                        fpm_f,
                        fpkm_f,
                    )
                    print("\t".join(str(value) for value in values), file=output)

                elif gene_strand == "-":
                    values = (
                        chrom,
                        tx_start,
                        tx_end,
                        gene_name,
                        mrna_size,
                        gene_strand,
                        frag_count_r,
                        fpm_r,
                        fpkm_r,
                    )
                    print("\t".join(str(value) for value in values), file=output)

                gene_finished += 1
                print(
                    f"\r{gene_finished} transcripts finished",
                    end="",
                    file=sys.stderr,
                )

        print(file=sys.stderr)
        print(f"Created {output_path}", file=sys.stderr)

    except (OSError, ValueError, IndexError, ZeroDivisionError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
