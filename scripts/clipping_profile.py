#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Estimate the clipping profile of RNA-seq reads from a BAM or SAM file.
The alignments must contain soft-clipping (``S``) CIGAR operations.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import SAM
from qcmodule import plotting
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=Path,
        help="Alignment file in BAM or SAM format.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        required=True,
        type=Path,
        help="Prefix for generated output files.",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        dest="map_qual",
        type=int,
        default=30,
        metavar="INT",
        help=(
            "Minimum mapping quality for an alignment to be considered "
            "uniquely mapped. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-s",
        "--sequencing",
        dest="layout",
        required=True,
        choices=("SE", "PE"),
        type=str.upper,
        help="Sequencing layout: SE for single-end or PE for paired-end.",
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Write the data files (and the R script) but do not draw plots.",
    )
    parser.add_argument(
        "--rscript",
        default="Rscript",
        metavar="PATH",
        help="Ignored: plots are drawn in Python since RSeQC 5.1.0 (R is not needed).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(
    parser: argparse.ArgumentParser,
    input_file: Path,
    output_prefix: Path,
    map_qual: int,
) -> None:
    """Validate command-line arguments."""
    if not input_file.is_file():
        parser.error(
            f"input alignment file does not exist or is not a file: {input_file}"
        )

    if input_file.suffix.lower() not in {".bam", ".sam", ".cram"}:
        parser.error(
            "input file must have a .bam, .sam, or .cram extension: "
            f"{input_file}"
        )

    if map_qual < 0:
        parser.error("--mapq must be zero or greater")

    output_parent = output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def main(argv: Sequence[str] | None = None) -> int:
    """Generate an RNA-seq clipping profile."""
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_args(
        parser=parser,
        input_file=args.input_file,
        output_prefix=args.out_prefix,
        map_qual=args.map_qual,
    )

    try:
        alignment = SAM.ParseBAM(str(args.input_file))
        plot_data = alignment.clipping_profile(
            outfile=str(args.out_prefix),
            q_cut=args.map_qual,
            type="S",
            PE=args.layout == "PE",
        )
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    if not args.skip_plot:
        plotting.run(parser, lambda: plotting.position_profile(plot_data, {key: f"{args.out_prefix}.clipping_profile{'.' + key if key else ''}.pdf" for key in plot_data}, title="clipping profile", ylabel="Non-clipped %", show="other"))

    return 0


if __name__ == "__main__":
    sys.exit(main())
