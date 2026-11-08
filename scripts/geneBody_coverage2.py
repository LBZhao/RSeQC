#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate RNA-seq coverage across the gene body (5' -> 3') from a BigWig.

Each transcript's exonic sequence is split into 100 equal slices and the
mean BigWig signal of every slice is measured (bases without data count as
0; introns are never read). By default the raw slice means are summed over
all transcripts, as in earlier releases. Profiles are also reported by mRNA
length, together with a 3'/5' ratio and a uniformity CV.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path
from typing import Sequence

import numpy as np

from qcmodule import genebody

__version__ = version("RSeQC")
__author__ = "Liguo Wang, Santiago Revale"


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        type=Path,
        help="Coverage signal file in BigWig format.",
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="ref_gene_model",
        required=True,
        type=Path,
        help="Reference gene model in BED12 format.",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        required=True,
        type=Path,
        help="Prefix for generated output files.",
    )
    parser.add_argument(
        "-t",
        "--graph-type",
        dest="graph_type",
        choices=("pdf", "png", "bmp", "jpeg", "tiff"),
        default="pdf",
        help="Plot file type. Default: %(default)s",
    )
    parser.add_argument(
        "-l",
        "--minimum-length",
        dest="min_mrna_length",
        type=int,
        default=100,
        metavar="INT",
        help=(
            "Minimum mRNA (exonic) length in bp. Values below 100 are not "
            "allowed. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--transcript-norm",
        choices=("none", "mean"),
        default="none",
        help=(
            "How transcripts are combined. 'none': sum raw signal across "
            "transcripts (as in earlier releases). 'mean': scale each "
            "transcript's profile to its own mean first, so every transcript "
            "has equal weight. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help=(
            "Skip transcripts whose mean exonic signal is below this. "
            "Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--length-bins",
        default="1000,2000,4000",
        metavar="EDGES",
        help=(
            "mRNA length bin edges in bp, comma separated (half-open bins, "
            "e.g. '1000,2000,4000' gives <1kb, 1-2kb, 2-4kb, >=4kb). 'auto' "
            "uses the quartiles of transcript lengths in the BED file; 'none' "
            "disables length binning. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Generate data and R code but do not execute the R script.",
    )
    parser.add_argument(
        "--rscript",
        default="Rscript",
        metavar="PATH",
        help="Rscript executable to use. Default: %(default)s",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate input files and output location."""
    if not args.input_file.is_file():
        parser.error(f"BigWig file does not exist: {args.input_file}")

    if not args.ref_gene_model.is_file():
        parser.error(
            f"reference BED12 file does not exist: {args.ref_gene_model}"
        )

    if args.min_mrna_length < 100:
        parser.error("--minimum-length cannot be smaller than 100")

    if args.min_coverage < 0:
        parser.error("--min-coverage cannot be negative")

    try:
        args.length_edges = genebody.parse_length_edges(args.length_bins)
    except ValueError as exc:
        parser.error(f"--length-bins: {exc}")

    output_parent = args.output_prefix.parent
    if not output_parent.exists():
        parser.error(f"output directory does not exist: {output_parent}")

    if not output_parent.is_dir():
        parser.error(f"output parent is not a directory: {output_parent}")


def sample_name(path: Path) -> str:
    """Sample label for output tables: the BigWig file name without suffix."""
    name = Path(path).name
    for suffix in (".bigwig", ".bigWig", ".BigWig", ".bw"):
        if name.endswith(suffix):
            return name[: -len(suffix)] or "sample"
    return Path(path).stem or "sample"


def resolve_length_edges(spec: list[int] | str, lengths: Sequence[int]) -> list[int]:
    """Turn a parsed ``--length-bins`` value into concrete bin edges."""
    if spec == "none":
        return []
    if spec == "auto":
        return genebody.auto_length_edges(lengths)
    return list(spec)


def write_coverage_table(path: Path, profile: Sequence[float]) -> None:
    """Write the 100-value profile (same layout as earlier releases)."""
    with path.open("w", encoding="utf-8") as output:
        output.write("percentile\tcount\n")
        for index, value in enumerate(profile):
            output.write(f"{index}\t{float(value)}\n")


def write_r_script(
    path: Path,
    profile: Sequence[float],
    n_transcripts: int,
    plot_file: str,
    graph_type: str,
    ylab: str = "average wigsum",
) -> None:
    """Write the R plotting script (same plot as earlier releases)."""
    y_values = ",".join(str(float(value)) for value in profile)
    with path.open("w", encoding="utf-8") as output:
        output.write(f"{graph_type}('{plot_file}')\n")
        output.write("x=1:100\n")
        output.write(f"y=c({y_values})\n")
        output.write(
            "plot("
            f"x, y/{max(n_transcripts, 1)}, "
            "xlab=\"percentile of gene body (5'->3')\", "
            f"ylab='{ylab}', "
            "type='s'"
            ")\n"
        )
        output.write("dev.off()\n")


def coverage_gene_body_bigwig(
    bigwig_file: Path,
    refbed: Path,
    output_prefix: Path,
    *,
    graph_type: str = "png",
    min_mrna_length: int = 100,
    length_bins: list[int] | str = "none",
    transcript_norm: str = "none",
    min_coverage: float = 0.0,
) -> tuple[list[Path], Path, int]:
    """Calculate gene-body coverage from a BigWig and write all outputs.

    Returns the data files written, the R script, and the number of
    transcripts that contributed to the profile.
    """
    bed = genebody.parse_bed12(refbed, min_mrna_length)
    print(
        f"Total {bed.n_records} transcripts loaded; "
        f"{len(bed.transcripts)} kept "
        f"({bed.n_duplicate} duplicates, {bed.n_short} shorter than "
        f"{min_mrna_length} bp, {bed.n_malformed} malformed lines)",
        file=sys.stderr,
    )
    if not bed.transcripts:
        raise ValueError("no usable transcripts in the BED file")

    lengths = [tx.length for tx in bed.transcripts]
    length_edges = resolve_length_edges(length_bins, lengths)
    labels = genebody.length_bin_labels(length_edges)

    print("Calculating coverage over gene body ...", file=sys.stderr)
    result = genebody.compute_bigwig_sample(
        bigwig_file,
        sample_name(bigwig_file),
        bed.transcripts,
        length_edges=length_edges,
        transcript_norm=transcript_norm,
        min_coverage=min_coverage,
        progress=sys.stderr,
    )

    print(
        f"{result.overall.n} transcripts used, "
        f"{result.n_low_coverage} below --min-coverage, "
        f"{result.n_missing_chrom} on chromosomes absent from the BigWig, "
        f"{result.n_out_of_range} extending past the chromosome end",
        file=sys.stderr,
    )
    if result.overall.n == 0:
        if result.n_missing_chrom == len(bed.transcripts):
            raise ValueError(
                "none of the BED chromosomes are in the BigWig; check "
                "chromosome naming (e.g. 'chr1' vs '1')"
            )
        raise ValueError("no transcript passed the coverage filters")

    profile = result.overall.result(transcript_norm)

    data_file = Path(f"{output_prefix}.geneBodyCoverage.txt")
    metrics_file = Path(f"{output_prefix}.geneBodyCoverage.metrics.txt")
    length_file = Path(f"{output_prefix}.geneBodyCoverage.byLength.txt")
    r_script = Path(f"{output_prefix}.geneBodyCoverage_plot.r")

    write_coverage_table(data_file, profile)
    written = [data_file]

    genebody.write_metrics(
        metrics_file, [result], labels if length_edges else [], transcript_norm
    )
    written.append(metrics_file)

    if length_edges:
        genebody.write_length_profiles(
            length_file, [result], labels, transcript_norm
        )
        written.append(length_file)

    print(
        f"\t3'/5' ratio: {genebody.fmt(genebody.ratio_3to5(profile))}\t"
        f"Uniformity CV: {genebody.fmt(genebody.uniformity_cv(profile))}",
        file=sys.stderr,
    )

    write_r_script(
        r_script,
        profile,
        result.overall.n,
        f"{output_prefix}.geneBodyCoverage.{graph_type}",
        graph_type,
        ylab=(
            "average wigsum"
            if transcript_norm == "none"
            else "signal (relative to transcript mean)"
        ),
    )

    return written, r_script, result.overall.n


def run_r_script(
    parser: argparse.ArgumentParser,
    script_path: Path,
    executable: str,
) -> None:
    """Execute the generated R script safely."""
    resolved = shutil.which(executable)

    if resolved is None:
        parser.exit(
            1,
            f'{parser.prog}: error: Rscript executable not found: '
            f'"{executable}"\n',
        )

    try:
        subprocess.run(
            [resolved, str(script_path)],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        parser.exit(
            exc.returncode or 1,
            f"{parser.prog}: error: R plotting failed for {script_path}\n",
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run BigWig-based gene-body coverage analysis."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    try:
        data_files, r_script, _n_used = coverage_gene_body_bigwig(
            bigwig_file=args.input_file,
            refbed=args.ref_gene_model,
            output_prefix=args.output_prefix,
            graph_type=args.graph_type,
            min_mrna_length=args.min_mrna_length,
            length_bins=args.length_edges,
            transcript_norm=args.transcript_norm,
            min_coverage=args.min_coverage,
        )

        if not args.skip_plot:
            run_r_script(
                parser=parser,
                script_path=r_script,
                executable=args.rscript,
            )

        for data_file in data_files:
            print(f"Created: {data_file}", file=sys.stderr)
        print(f"Created: {r_script}", file=sys.stderr)

    except (OSError, ValueError, RuntimeError, ZeroDivisionError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
