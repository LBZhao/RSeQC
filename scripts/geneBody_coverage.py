#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Calculate RNA-seq read coverage across the gene body (5' -> 3').

Each transcript's exonic sequence is split into 100 equal slices and the
mean read coverage of every slice is measured. Profiles are aggregated over
all transcripts and, separately, by mRNA length.

Notes
-----
1. Input BAM files must be coordinate-sorted and indexed.
2. Transcripts shorter than ``--minimum-length`` (at least 100 bp) are skipped.
3. Coverage counts aligned bases only (no introns, deletions, or soft clips).
   Overlapping mates of a pair are counted once.
4. Reads that are unmapped, secondary, supplementary, QC-failed, duplicates,
   or below ``--mapq`` are ignored.
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Sequence

import numpy as np

from qcmodule import genebody
from qcmodule import getBamFiles

__version__ = version("RSeQC")
__author__ = "Liguo Wang"

_LOG_FILE: Path | None = None


def valid_name(value: str) -> str:
    """Convert a string into a valid R variable name."""
    symbols = set(
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "0123456789_."
    )
    digits = set("0123456789")

    result = "_".join(value.split())

    if not result:
        return "sample"

    if result[0] in digits:
        result = "V" + result

    return "".join(character if character in symbols else "_" for character in result)


def printlog(message: str) -> None:
    """Write a timestamped message to stderr and the run's log file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted = f"@ {timestamp}: {message}"

    print(formatted, file=sys.stderr)

    if _LOG_FILE is not None:
        with _LOG_FILE.open("a", encoding="utf-8") as handle:
            print(formatted, file=handle)


def make_unique_sample_name(name: str, seen: list[str]) -> str:
    """Return a unique R-safe sample name using the original naming scheme."""
    count = seen.count(name)
    return name if count == 0 else f"{name}.{count}"


fmt = genebody.fmt

# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help=(
            "A BAM file, comma-separated BAM files, a directory containing "
            "BAM files, or a text file listing BAM paths."
        ),
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
        help="Prefix for output files.",
    )
    parser.add_argument(
        "-l",
        "--minimum-length",
        "--minimum_length",
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
        "-q",
        "--mapq",
        dest="min_mapq",
        type=int,
        default=30,
        metavar="INT",
        help="Minimum mapping quality of reads to count. Default: %(default)s",
    )
    parser.add_argument(
        "--transcript-norm",
        choices=("mean", "none"),
        default="mean",
        help=(
            "How transcripts are combined. 'mean': scale each transcript's "
            "profile to its own mean before averaging, so every transcript "
            "has equal weight. 'none': sum raw coverage across transcripts "
            "(the behavior of RSeQC <= 5.0.6; highly expressed genes "
            "dominate). Default: %(default)s"
        ),
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=None,
        metavar="FLOAT",
        help=(
            "Skip transcripts whose mean per-base coverage is below this. "
            "Default: 5 with --transcript-norm mean, 0 with none."
        ),
    )
    parser.add_argument(
        "--length-bins",
        default="1000,2000,4000",
        metavar="EDGES",
        help=(
            "mRNA length bin edges in bp, comma separated (bins are "
            "half-open, e.g. '1000,2000,4000' gives <1kb, 1-2kb, 2-4kb, "
            ">=4kb). 'auto' uses the quartiles of transcript lengths in the "
            "BED file; 'none' disables length binning. Default: %(default)s"
        ),
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        metavar="INT",
        help="Number of BAM files to process in parallel. Default: %(default)s",
    )
    parser.add_argument(
        "--per-transcript",
        action="store_true",
        help=(
            "Also write each transcript's 100-slice profile to "
            "PREFIX.SAMPLE.geneBodyCoverage.perTranscript.txt."
        ),
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="output_format",
        choices=("pdf", "png", "jpeg"),
        default="pdf",
        help="Plot output format. Default: %(default)s",
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
    """Validate command-line arguments."""
    if not args.ref_gene_model.is_file():
        parser.error(
            f"reference BED12 file does not exist: {args.ref_gene_model}"
        )

    if args.min_mrna_length < 100:
        parser.error("--minimum-length cannot be smaller than 100")

    if args.min_mapq < 0:
        parser.error("--mapq cannot be negative")

    if args.processes < 1:
        parser.error("--processes must be at least 1")

    if args.min_coverage is None:
        args.min_coverage = 5.0 if args.transcript_norm == "mean" else 0.0
    elif args.min_coverage < 0:
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


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------


def write_profiles(
    path: Path,
    results: Sequence[genebody.SampleResult],
    transcript_norm: str,
) -> None:
    """Write the overall profile table (same layout as earlier releases)."""
    with path.open("w", encoding="utf-8") as output:
        print(
            "Percentile\t" + "\t".join(str(i) for i in range(1, 101)),
            file=output,
        )
        for result in results:
            profile = result.overall.result(transcript_norm)
            print(
                result.name + "\t" + "\t".join(fmt(v) for v in profile),
                file=output,
            )


write_length_profiles = genebody.write_length_profiles
write_metrics = genebody.write_metrics


def relative_profile(profile: np.ndarray) -> np.ndarray | None:
    """Scale a profile to its own mean (1.0 = uniform coverage) for plots."""
    mean = float(np.mean(profile))
    if mean <= 0 or not math.isfinite(mean):
        return None
    return profile / mean


def _r_string(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _r_vector(values: Sequence[float]) -> str:
    return "c(" + ",".join(fmt(v) for v in values) + ")"


def write_r_code(
    results: Sequence[genebody.SampleResult],
    labels: Sequence[str],
    transcript_norm: str,
    file_prefix: Path,
    *,
    output_format: str = "pdf",
    min_transcripts: int = 100,
) -> Path:
    """Generate an R script that draws the curves, heatmap and length panels.

    All curves are scaled to each profile's own mean, so a flat line at 1.0
    means perfectly even coverage.
    """
    r_script = Path(f"{file_prefix}.r")
    ext = output_format.lower()
    device = {"pdf": "pdf", "png": "png", "jpeg": "jpeg"}[ext]

    overall: list[tuple[str, np.ndarray]] = []
    for result in results:
        rel = relative_profile(result.overall.result(transcript_norm))
        if rel is not None:
            overall.append((result.name, rel))
    overall.sort(key=lambda item: genebody.bias_rank_key(item[1]))

    lines: list[str] = []
    add = lines.append

    add("# Generated by geneBody_coverage.py (RSeQC %s)" % __version__)
    add("# Curves are scaled to each profile's mean: 1.0 = uniform coverage.")
    add("profiles <- list()")
    for name, rel in overall:
        add(f"profiles[[{_r_string(name)}]] <- {_r_vector(rel)}")

    add("by_length <- list()")
    for index, label in enumerate(labels):
        add(f"by_length[[{_r_string(label)}]] <- list()")
        add(f"by_length_n_{index} <- c()")
        for result in results:
            profile_sum = result.by_length[index]
            if not profile_sum.n:
                continue
            rel = relative_profile(profile_sum.result(transcript_norm))
            if rel is None:
                continue
            add(
                f"by_length[[{_r_string(label)}]][[{_r_string(result.name)}]]"
                f" <- {_r_vector(rel)}"
            )
            add(
                f"by_length_n_{index}[{_r_string(result.name)}] <- "
                f"{profile_sum.n}"
            )

    add(f"""
x <- 1:100
sample_names <- names(profiles)
cols <- colorRampPalette(c("#1b9e77","#d95f02","#7570b3","#e7298a",
                           "#66a61e","#e6ab02","#a6761d","#666666"))(max(1, length(sample_names)))
names(cols) <- sample_names
xlab <- "Gene body percentile (5'->3')"
ylab <- "Coverage (relative to mean)"

draw_curves <- function(pl, main) {{
  if (length(pl) == 0) {{
    plot.new(); title(main = paste(main, "(no data)")); return(invisible())
  }}
  ymax <- max(1.2, max(unlist(pl), na.rm = TRUE) * 1.05)
  plot(x, pl[[1]], type = "n", ylim = c(0, ymax), xlab = xlab, ylab = ylab,
       main = main)
  abline(h = 1, lty = 3, col = "grey50")
  for (nm in names(pl)) lines(x, pl[[nm]], col = cols[[nm]], lwd = 1.2)
}}

draw_legend_inside <- function(nms) {{
  legend("bottom", legend = nms, col = cols[nms], lwd = 2, bty = "n",
         cex = 0.8, ncol = min(3, length(nms)))
}}

## Overall curves
{device}({_r_string(f"{file_prefix}.curves.{ext}")})
if (length(sample_names) > 6) {{
  layout(matrix(c(1, 1, 1, 2), 1, 4))
  draw_curves(profiles, "Gene body coverage")
  par(mar = c(1, 0, 2, 1)); plot.new()
  legend("topleft", legend = sample_names, col = cols, lwd = 2, bty = "n",
         cex = 0.8)
}} else {{
  draw_curves(profiles, "Gene body coverage")
  if (length(sample_names) > 0) draw_legend_inside(sample_names)
}}
invisible(dev.off())

## Heatmap (3 or more samples; best coverage at the top)
if (length(sample_names) >= 3) {{
  data_matrix <- do.call(rbind, profiles)
  data_matrix <- data_matrix[nrow(data_matrix):1, , drop = FALSE]
  col_labels <- rep("", 100)
  col_labels[c(1, seq(10, 100, 10))] <- c(1, seq(10, 100, 10))
  {device}({_r_string(f"{file_prefix}.heatMap.{ext}")})
  heatmap(data_matrix, scale = "none", Rowv = NA, Colv = NA,
          labCol = col_labels, col = cm.colors(256), margins = c(6, 8),
          cexRow = 1, cexCol = 1, xlab = xlab)
  invisible(dev.off())
}}

## Curves by mRNA length
if (length(by_length) > 1) {{
  n_panels <- length(by_length)
  n_col <- if (n_panels <= 2) n_panels else 2
  n_row <- ceiling(n_panels / n_col)
  {device}({_r_string(f"{file_prefix}.byLength.{ext}")},
           width = if ({_r_string(device)} == "pdf") 5 * n_col else 480 * n_col,
           height = if ({_r_string(device)} == "pdf") 4.5 * n_row else 432 * n_row)
  par(mfrow = c(n_row, n_col))
  for (i in seq_along(by_length)) {{
    label <- names(by_length)[i]
    counts <- get(paste0("by_length_n_", i - 1))
    title_text <- label
    if (length(counts) > 0) {{
      n_text <- if (min(counts) == max(counts)) min(counts) else
        paste(range(counts), collapse = "-")
      title_text <- paste0(label, "  (n = ", n_text, " transcripts)")
      if (min(counts) < {min_transcripts}) title_text <- paste(title_text, "[low n]")
    }}
    draw_curves(by_length[[label]], title_text)
    if (i == 1 && length(sample_names) <= 6 && length(sample_names) > 0)
      draw_legend_inside(sample_names)
  }}
  invisible(dev.off())
}}
""")

    r_script.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return r_script


def run_r_script(
    parser: argparse.ArgumentParser,
    r_script: Path,
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
            [resolved, str(r_script)],
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        parser.exit(
            exc.returncode or 1,
            f"{parser.prog}: error: R plotting failed for {r_script}\n",
        )


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main(argv: Sequence[str] | None = None) -> int:
    """Run gene-body coverage analysis."""
    global _LOG_FILE

    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)

    prefix = args.output_prefix
    _LOG_FILE = Path(f"{prefix}.geneBodyCoverage.log.txt")
    coverage_file = Path(f"{prefix}.geneBodyCoverage.txt")
    length_file = Path(f"{prefix}.geneBodyCoverage.byLength.txt")
    metrics_file = Path(f"{prefix}.geneBodyCoverage.metrics.txt")
    plot_prefix = Path(f"{prefix}.geneBodyCoverage")

    try:
        printlog("Read BED file (reference gene model) ...")
        bed = genebody.parse_bed12(
            args.ref_gene_model, args.min_mrna_length
        )
        printlog(
            f"Total {bed.n_records} transcripts loaded; "
            f"{len(bed.transcripts)} kept "
            f"({bed.n_duplicate} duplicates, {bed.n_short} shorter than "
            f"{args.min_mrna_length} bp, {bed.n_malformed} malformed lines)"
        )
        if not bed.transcripts:
            raise ValueError("no usable transcripts in the BED file")

        lengths = [tx.length for tx in bed.transcripts]
        q25, q50, q75 = np.quantile(lengths, [0.25, 0.5, 0.75])
        printlog(
            "mRNA length quartiles (bp): "
            f"25%={q25:.0f}, median={q50:.0f}, 75%={q75:.0f}"
        )

        if args.length_edges == "none":
            length_edges: list[int] = []
        elif args.length_edges == "auto":
            length_edges = genebody.auto_length_edges(lengths)
        else:
            length_edges = list(args.length_edges)
        labels = genebody.length_bin_labels(length_edges)

        if length_edges:
            bin_counts = np.bincount(
                [genebody.assign_length_bin(n, length_edges) for n in lengths],
                minlength=len(labels),
            )
            printlog(
                "Transcripts per length bin: "
                + ", ".join(
                    f"{label}={count}" for label, count in zip(labels, bin_counts)
                )
            )

        clusters = genebody.build_clusters(bed.transcripts)

        printlog("Get BAM file(s) ...")
        bam_files = [
            Path(path)
            for path in getBamFiles.get_bam_files(args.input_files)
        ]
        bam_files = [path for path in bam_files if path.is_file()]
        if not bam_files:
            raise ValueError("no BAM files were found")
        for bam_file in bam_files:
            print(f"\t{bam_file}", file=sys.stderr)

        jobs = []
        seen_names: list[str] = []
        for bam_file in bam_files:
            base_name = valid_name(bam_file.stem)
            sample_name = make_unique_sample_name(base_name, seen_names)
            seen_names.append(base_name)
            jobs.append(
                dict(
                    bam_path=bam_file,
                    name=sample_name,
                    length_edges=length_edges,
                    transcript_norm=args.transcript_norm,
                    min_coverage=args.min_coverage,
                    min_mapq=args.min_mapq,
                    per_transcript=(
                        Path(
                            f"{prefix}.{sample_name}"
                            ".geneBodyCoverage.perTranscript.txt"
                        )
                        if args.per_transcript
                        else None
                    ),
                )
            )

        if args.processes > 1 and len(jobs) > 1:
            printlog(
                f"Processing {len(jobs)} BAM files with "
                f"{min(args.processes, len(jobs))} processes ..."
            )
            with ProcessPoolExecutor(
                max_workers=min(args.processes, len(jobs)),
                initializer=genebody.init_worker,
                initargs=(clusters,),
            ) as pool:
                raw_results = list(pool.map(genebody.run_worker_job, jobs))
        else:
            raw_results = []
            for job in jobs:
                printlog(f"Processing {job['bam_path'].name} ...")
                raw_results.append(
                    genebody.compute_sample(
                        clusters=clusters, progress=sys.stderr, **job
                    )
                )

        results = []
        for result in raw_results:
            if result.n_missing_chrom == len(bed.transcripts):
                print(
                    f"Warning: none of the BED chromosomes are in "
                    f"{result.bam}; check chromosome naming (e.g. 'chr1' vs "
                    f"'1'). Skipped.",
                    file=sys.stderr,
                )
                continue
            if result.overall.n == 0:
                print(
                    f"Warning: no transcript passed the coverage filters in "
                    f"{result.bam}. Skipped.",
                    file=sys.stderr,
                )
                continue
            printlog(
                f"{result.name}: {result.overall.n} transcripts used, "
                f"{result.n_low_coverage} below --min-coverage, "
                f"{result.n_missing_chrom} on chromosomes absent from the BAM"
            )
            for label, profile_sum in zip(labels, result.by_length):
                if length_edges and profile_sum.n < 100:
                    print(
                        f"Warning: {result.name} length bin {label} has only "
                        f"{profile_sum.n} transcripts; its curve is unreliable.",
                        file=sys.stderr,
                    )
            results.append(result)

        if not results:
            raise ValueError("no valid coverage profiles were generated")

        write_profiles(coverage_file, results, args.transcript_norm)
        write_metrics(metrics_file, results, labels, args.transcript_norm)
        if length_edges:
            write_length_profiles(
                length_file, results, labels, args.transcript_norm
            )

        print("\n\tSample\t3'/5' ratio\tUniformity CV", file=sys.stderr)
        for result in results:
            profile = result.overall.result(args.transcript_norm)
            print(
                f"\t{result.name}\t{fmt(genebody.ratio_3to5(profile))}\t"
                f"{fmt(genebody.uniformity_cv(profile))}",
                file=sys.stderr,
            )

        r_script = write_r_code(
            results,
            labels if length_edges else [],
            args.transcript_norm,
            plot_prefix,
            output_format=args.output_format,
        )

        if not args.skip_plot:
            printlog("Running R script ...")
            run_r_script(parser, r_script, args.rscript)

        printlog(f"Created: {coverage_file}")
        printlog(f"Created: {metrics_file}")
        if length_edges:
            printlog(f"Created: {length_file}")
        printlog(f"Created: {r_script}")

    except (OSError, ValueError) as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
