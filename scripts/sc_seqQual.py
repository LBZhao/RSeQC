#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Generate sequencing-quality matrices and a heatmap from a FASTQ file."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Sequence
from qcmodule import fastq
from qcmodule import heatmap
from importlib.metadata import version

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


DESCRIPTION = (
    "Generate per-cycle Phred quality-count and quality-percentage matrices "
    "from a FASTQ file, then optionally render a heatmap."
)

EPILOG = """
Examples
--------
Generate matrices and a PDF heatmap:
    sc_seqQual.py \
        -i reads.fastq.gz \
        -o sample

Generate matrices only:
    sc_seqQual.py \
        -i reads.fastq.gz \
        -o sample \
        --skip-heatmap


Notes
-----
* Input must be FASTQ.
* The percentage matrix is calculated independently for each read cycle.
* Heatmaps are drawn with matplotlib; R is not needed.
"""


SUPPORTED_FILE_TYPES = ("pdf", "png", "tiff", "bmp", "jpeg")
SUPPORTED_ANGLES = (0, 45, 90, 270, 315)


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--infile",
        dest="in_file",
        required=True,
        type=Path,
        help="Input FASTQ file; compressed input is supported by qcmodule.fastq.",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="out_file",
        required=True,
        type=Path,
        help="Prefix for generated matrix and heatmap files.",
    )
    parser.add_argument(
        "-n",
        "--nseq-limit",
        dest="max_seq",
        type=int,
        default=None,
        metavar="INT",
        help="Maximum number of sequences to process. Default: all",
    )
    parser.add_argument(
        "--cell-width",
        type=int,
        default=12,
        metavar="INT",
        help="Heatmap cell width in points. Default: %(default)s",
    )
    parser.add_argument(
        "--cell-height",
        type=int,
        default=10,
        metavar="INT",
        help="Heatmap cell height in points. Default: %(default)s",
    )
    parser.add_argument(
        "--font-size",
        type=int,
        default=6,
        metavar="INT",
        help="Heatmap font size in points. Default: %(default)s",
    )
    parser.add_argument(
        "--angle",
        dest="col_angle",
        type=int,
        choices=SUPPORTED_ANGLES,
        default=45,
        help="Column-label angle. Default: %(default)s",
    )
    parser.add_argument(
        "--text-color",
        default="black",
        metavar="COLOR",
        help="Color of numbers displayed in heatmap cells. Default: %(default)s",
    )
    parser.add_argument(
        "--file-type",
        type=str.lower,
        choices=SUPPORTED_FILE_TYPES,
        default="pdf",
        help="Heatmap output format. Default: %(default)s",
    )
    parser.add_argument(
        "--no-num",
        action="store_true",
        help="Do not print numerical values inside heatmap cells.",
    )
    parser.add_argument(
        "--skip-heatmap",
        action="store_true",
        help="Generate quality matrices without rendering a heatmap.",
    )
    parser.add_argument(
        "--rscript",
        default="Rscript",
        metavar="PATH",
        help="Ignored: heatmaps are drawn in Python since RSeQC 5.1.0+lbzhao (R is not needed).",
    )
    parser.add_argument(
        "--install-r-deps",
        action="store_true",
        help="Ignored: heatmaps are drawn in Python since RSeQC 5.1.0+lbzhao (R is not needed).",
    )
    parser.add_argument(
        "--cran-mirror",
        default="https://cloud.r-project.org",
        metavar="URL",
        help="Ignored: heatmaps are drawn in Python since RSeQC 5.1.0+lbzhao (R is not needed).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed progress logging.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def configure_logging(verbose: bool) -> None:
    """Configure command-line logging."""
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if verbose else logging.INFO,
    )


def validate_args(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> None:
    """Validate command-line arguments."""
    if not args.in_file.is_file():
        parser.error(f"input FASTQ file does not exist: {args.in_file}")

    if args.max_seq is not None and args.max_seq <= 0:
        parser.error("--nseq-limit must be greater than zero")

    if args.cell_width <= 0:
        parser.error("--cell-width must be greater than zero")

    if args.cell_height <= 0:
        parser.error("--cell-height must be greater than zero")

    if args.font_size <= 0:
        parser.error("--font-size must be greater than zero")

    if not args.text_color.strip():
        parser.error("--text-color cannot be empty")

    if not args.cran_mirror.strip():
        parser.error("--cran-mirror cannot be empty")

    output_parent = args.out_file.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def expected_output_paths(
    output_prefix: Path,
    file_type: str,
) -> tuple[Path, Path, Path]:
    """Return expected count, percentage, and heatmap output paths."""
    count_matrix = Path(f"{output_prefix}.qual_count.csv")
    percentage_matrix = Path(f"{output_prefix}.qual_percent.csv")
    heatmap_file = Path(f"{output_prefix}.qual_heatmap.{file_type}")

    return count_matrix, percentage_matrix, heatmap_file


def main(argv: Sequence[str] | None = None) -> int:
    """Generate sequencing-quality matrices and heatmap."""
    parser = build_parser()
    args = parser.parse_args(argv)

    configure_logging(args.verbose)
    validate_args(parser, args)

    count_path, percent_path, heatmap_path = expected_output_paths(
        args.out_file,
        args.file_type,
    )

    try:
        file_iter = fastq.fastq_iter(
            str(args.in_file),
            mode="qual",
        )
        quality_matrix = fastq.qual2countMat(
            file_iter,
            limit=args.max_seq,
        )

        if quality_matrix is None or quality_matrix.empty:
            raise ValueError(
                "no usable quality records were found in the FASTQ file"
            )

        quality_matrix = quality_matrix.T
        quality_matrix.sort_index(
            inplace=True,
            ascending=False,
        )

        logging.debug(
            "Sequence quality score matrix (raw read counts):\n%s",
            quality_matrix,
        )

        column_totals = quality_matrix.sum(axis=0)
        zero_total_columns = column_totals[column_totals == 0].index.tolist()

        if zero_total_columns:
            raise ValueError(
                "one or more read cycles contain zero total observations: "
                + ", ".join(map(str, zero_total_columns))
            )

        quality_percent = quality_matrix.div(
            column_totals,
            axis=1,
        )

        logging.debug(
            "Sequence quality score matrix (fraction of reads):\n%s",
            quality_percent,
        )

        quality_matrix.to_csv(
            count_path,
            index=True,
            index_label="Index",
        )
        quality_percent.to_csv(
            percent_path,
            index=True,
            index_label="Index",
        )

        if not count_path.is_file():
            raise OSError(f"quality-count matrix was not created: {count_path}")

        if not percent_path.is_file():
            raise OSError(
                f"quality-percentage matrix was not created: {percent_path}"
            )

        logging.info("Created %s", count_path)
        logging.info("Created %s", percent_path)

        if not args.skip_heatmap:

            if heatmap_path.exists():
                heatmap_path.unlink()

            heatmap.make_heatmap(
                infile=str(percent_path),
                outfile=str(Path(f"{args.out_file}.qual_heatmap")),
                filetype=args.file_type,
                cell_width=args.cell_width,
                cell_height=args.cell_height,
                col_angle=args.col_angle,
                font_size=args.font_size,
                text_color=args.text_color,
                no_numbers=args.no_num,
            )

            if not heatmap_path.is_file():
                raise RuntimeError(
                    f"heatmap generation failed; expected output was not "
                    f"created: {heatmap_path}"
                )

            logging.info("Created %s", heatmap_path)

    except (
        OSError,
        ValueError,
        RuntimeError,
        KeyError,
        IndexError,
    ) as exc:
        logging.error("%s", exc)
        return 1

    logging.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
