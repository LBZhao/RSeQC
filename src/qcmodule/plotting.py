# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Plots for RSeQC tools, drawn with matplotlib (no R needed).

Each function takes the data a tool computed and writes one or more image
files. The image format follows the file extension (pdf, png, jpeg, tiff,
svg; bmp is converted with Pillow). Before RSeQC 5.1.0+lbzhao the tools wrote R
scripts and ran them with Rscript; the R scripts are still written next to
the data for anyone who wants to restyle a plot in R, but they are no
longer run.
"""

from __future__ import annotations

import math
import statistics
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

_PIL_ONLY = {"bmp"}


def save(fig, path) -> Path:
    """Save *fig* to *path* (format from the extension) and close it."""
    path = Path(path)
    fmt = path.suffix.lower().lstrip(".") or "pdf"
    if fmt == "jpg":
        fmt = "jpeg"
    fig.tight_layout()
    if fmt in _PIL_ONLY:
        from io import BytesIO

        from PIL import Image

        buffer = BytesIO()
        fig.savefig(buffer, format="png", dpi=150)
        buffer.seek(0)
        Image.open(buffer).convert("RGB").save(path, format=fmt.upper())
    else:
        fig.savefig(path, format=fmt, dpi=150)
    plt.close(fig)
    return path


def _figure(width=7.0, height=7.0):
    return plt.subplots(figsize=(width, height))


# --------------------------------------------------------------------------
# Read-level QC
# --------------------------------------------------------------------------

_BASE_COLORS = {
    "A": "darkgreen", "T": "red", "G": "blue", "C": "cyan", "N": "black", "X": "grey",
}


def nucleotide_composition(data: Mapping, path) -> Path:
    """read_NVC: frequency of each base at each read position."""
    bases = "ATGCNX" if data.get("include_nx", True) else "ATGC"
    counts = {b: np.asarray(data["counts"][b], dtype=float) for b in "ACGTNX"}
    total = sum(counts[b] for b in bases)
    with np.errstate(invalid="ignore", divide="ignore"):
        freqs = {b: np.where(total > 0, counts[b] / total, np.nan) for b in bases}
    fig, ax = _figure()
    positions = np.arange(len(total))
    for base in bases:
        ax.plot(positions, freqs[base], "o-", ms=3, color=_BASE_COLORS[base], label=base)
    ax.set_xlabel("Position of Read")
    ax.set_ylabel("Nucleotide Frequency")
    ax.legend(loc="upper right", ncol=len(bases))
    return save(fig, path)


def quality_plots(data: Mapping, boxplot_path, heatmap_path) -> list[Path]:
    """read_quality: Phred score distribution per read position."""
    quality = data["quality"]
    positions = sorted(quality)
    if not positions:
        raise ValueError("no base qualities to plot")
    scores = sorted({q for pos in positions for q in quality[pos]})
    q_min, q_max = scores[0], scores[-1]

    stats = []
    for pos in positions:
        values = np.array(sorted(quality[pos]), dtype=float)
        weights = np.array([quality[pos][int(v)] for v in values], dtype=float)
        cdf = np.cumsum(weights) / weights.sum()

        def pct(p):
            return float(values[np.searchsorted(cdf, p)])

        q1, med, q3 = pct(0.25), pct(0.5), pct(0.75)
        iqr = q3 - q1
        stats.append({
            "med": med, "q1": q1, "q3": q3,
            "whislo": max(values[0], q1 - 1.5 * iqr),
            "whishi": min(values[-1], q3 + 1.5 * iqr),
            "label": str(pos + 1),
        })
    fig, ax = _figure(max(7.0, len(positions) * 0.08), 6)
    ax.bxp(stats, showfliers=False)
    step = max(1, len(positions) // 20)
    ax.set_xticks(range(1, len(positions) + 1, step))
    ax.set_xticklabels([str(p + 1) for p in positions[::step]])
    ax.set_xlabel("Position of Read (5'->3')")
    ax.set_ylabel("Phred Quality Score")
    paths = [save(fig, boxplot_path)]

    matrix = np.zeros((q_max - q_min + 1, len(positions)))
    for col, pos in enumerate(positions):
        for q, n in quality[pos].items():
            matrix[q - q_min, col] = n
    cmap = LinearSegmentedColormap.from_list("rseqc_qual", ["blue", "orange", "#cd0000", "red"])
    fig, ax = _figure(max(7.0, len(positions) * 0.08), 6)
    image = ax.imshow(matrix, aspect="auto", origin="lower", cmap=cmap,
                      extent=(0.5, len(positions) + 0.5, q_min - 0.5, q_max + 0.5))
    fig.colorbar(image, ax=ax, label="Number of bases")
    ax.set_xlabel("Position of Read")
    ax.set_ylabel("Phred Quality Score")
    paths.append(save(fig, heatmap_path))
    return paths


def gc_content(data: Mapping, path) -> Path:
    """read_GC: distribution of per-read GC content."""
    gc = data["gc_hist"]
    fig, ax = _figure()
    if gc:
        values = np.array(list(gc.keys()))
        weights = np.array(list(gc.values()), dtype=float)
        ax.hist(values, bins=100, weights=weights, density=True,
                histtype="step", color="blue")
    ax.set_xlabel("GC content (%)")
    ax.set_ylabel("Density of Reads")
    return save(fig, path)


def duplication(data: Mapping, path) -> Path:
    """read_duplication: number of reads by how often they occur."""
    fig, ax = _figure()
    for key, label, marker, color in (
        ("position", "Mapping-based", "x", "blue"),
        ("sequence", "Sequence-based", "o", "red"),
    ):
        occ = sorted(data[key])
        ax.plot(occ, [math.log10(data[key][o]) for o in occ], marker,
                ms=4, color=color, label=label, linestyle="none")
    ax.set_xlim(1, data.get("up_bound", 500))
    ax.set_xlabel("Occurrence of read")
    ax.set_ylabel("Number of Reads (log10)")
    ax.legend(loc="upper right")
    return save(fig, path)


def position_profile(
    data: Mapping[str, tuple],
    paths: Mapping[str, Path],
    *,
    title: str,
    ylabel: str,
    show: str,
) -> list[Path]:
    """clipping_profile / insertion_profile: percent of reads per position.

    *data* maps a mate label ("" for single-end, "R1", "R2") to
    ``(positions, event_counts, other_counts)``. *show* is ``"other"`` to plot
    ``other / (event + other)`` (non-clipped %) or ``"event"`` to plot
    ``event / (event + other)`` (insertion %).
    """
    written = []
    for label, (positions, events, others) in data.items():
        events = np.asarray(events, dtype=float)
        others = np.asarray(others, dtype=float)
        total = events + others
        numerator = others if show == "other" else events
        with np.errstate(invalid="ignore", divide="ignore"):
            percent = np.where(total > 0, numerator * 100 / total, np.nan)
        fig, ax = _figure()
        ax.plot(positions, percent, "o-", ms=3, color="blue")
        suffix = {"R1": " (read-1)", "R2": " (read-2)"}.get(label, "")
        ax.set_xlabel("Position of read" + suffix)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        written.append(save(fig, paths[label]))
    return written


_MISMATCH_COLORS = [
    "green", "powderblue", "lightseagreen", "red", "#8b2252", "#e066ff",
    "blue", "royalblue", "#63b8ff", "orange", "gold", "black",
]


def mismatch_profile(data: Mapping, path) -> Path:
    """mismatch_profile: mismatches by type and read position (log10)."""
    fig, ax = _figure()
    positions = data["positions"]
    for color, (genotype, counts) in zip(_MISMATCH_COLORS, data["counts"].items()):
        ax.plot(positions, np.log10(np.asarray(counts, dtype=float) + 1), color=color, label=genotype)
    ax.set_xlabel("Read position (5'->3')")
    ax.set_ylabel("log10(# of mismatch)")
    ax.legend(ncol=4, fontsize="small")
    return save(fig, path)


def deletion_profile(data: Mapping, path) -> Path:
    """deletion_profile: deletions by read position."""
    fig, ax = _figure()
    ax.plot(data["positions"], data["counts"], "o-", ms=3, color="blue")
    ax.set_xlabel("Read position (5'->3')")
    ax.set_ylabel("Deletion count")
    return save(fig, path)


# --------------------------------------------------------------------------
# Fragments and junctions
# --------------------------------------------------------------------------


def inner_distance_stats(data: Mapping) -> dict[str, float]:
    """Mean, median and SD of inner distances inside the plotting window.

    Uses the measured distances (R used the midpoints of the histogram bins).
    """
    low, up = data["low_bound"], data["up_bound"]
    values = [d for d in data["distances"] if low < d <= up]
    if not values:
        return {"mean": math.nan, "median": math.nan, "sd": math.nan, "n": 0}
    return {
        "mean": statistics.fmean(values),
        "median": statistics.median(values),
        "sd": statistics.stdev(values) if len(values) > 1 else math.nan,
        "n": len(values),
    }


def inner_distance(data: Mapping, path) -> Path:
    """inner_distance: histogram and density of inner distances."""
    low, up, step = data["low_bound"], data["up_bound"], data["step"]
    values = np.array([d for d in data["distances"] if low < d <= up], dtype=float)
    stats = inner_distance_stats(data)
    fig, ax = _figure()
    bins = np.arange(low, up + step, step)
    if values.size:
        ax.hist(values, bins=bins, density=True, histtype="step", color="blue")
        # Gaussian density with bandwidth 2 x step, as the R script used.
        # Binned at 1 bp and smoothed by convolution, so memory stays small
        # for millions of pairs.
        bandwidth = 2.0 * step
        grid = np.arange(math.floor(low), math.ceil(up) + 1)
        counts = np.bincount((values - grid[0]).astype(int), minlength=grid.size)[: grid.size]
        offsets = np.arange(-int(4 * bandwidth), int(4 * bandwidth) + 1)
        kernel = np.exp(-0.5 * (offsets / bandwidth) ** 2) / (bandwidth * math.sqrt(2 * math.pi))
        half = offsets.size // 2
        density = np.convolve(counts, kernel, mode="full")[half: half + grid.size] / values.size
        ax.plot(grid, density, color="red")
    ax.set_xlabel("mRNA insert size (bp)")
    ax.set_ylabel("Density")
    ax.set_title(f"Mean={stats['mean']:.2f}; SD={stats['sd']:.2f}")
    return save(fig, path)


_JUNCTION_COLORS = {
    "annotated": "#2166ac", "novel_combination": "#c51b7d",
    "partial_novel": "#d6604d", "complete_novel": "#4dac26",
}


def junction_pies(data: Mapping, events_path, junctions_path) -> list[Path]:
    """junction_annotation: share of each junction class."""
    written = []
    for key, title, path in (
        ("events", "splicing events", events_path),
        ("junctions", "splicing junctions", junctions_path),
    ):
        counts = data[key]
        total = sum(counts.values())
        labels = [f"{label} {round(100 * n / total) if total else 0}%" for label, n in counts.items()]
        fig, ax = _figure()
        nonzero = [(lab, n, _JUNCTION_COLORS.get(cls, "grey"))
                   for (cls, n), lab in zip(counts.items(), labels) if n > 0]
        if nonzero:
            ax.pie([n for _, n, _ in nonzero], labels=[lab for lab, _, _ in nonzero],
                   colors=[c for _, _, c in nonzero], startangle=30)
        ax.set_title(title)
        written.append(save(fig, path))
    return written


def junction_saturation(data: Mapping, path) -> Path:
    """junction_saturation: junctions detected vs. reads sampled."""
    fig, ax = _figure()
    x = data["percent"]
    for key, label, color in (("all", "All junctions", "blue"),
                              ("known", "known junctions", "red"),
                              ("novel", "novel junctions", "green")):
        ax.plot(x, np.asarray(data[key]) / 1000, "o-", color=color, label=label)
    ax.set_xlabel("percent of total reads")
    ax.set_ylabel("Number of splicing junctions (x1000)")
    ax.legend(loc="upper left")
    return save(fig, path)


def rpkm_saturation(groups: Mapping[str, Mapping[str, Sequence[float]]], path) -> Path:
    """RPKM_saturation: percent relative error by expression quartile."""
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    for ax, (quartile, by_percent) in zip(axes.flat, groups.items()):
        labels = list(by_percent)
        ax.boxplot([np.asarray(by_percent[k], dtype=float) * 100 for k in labels],
                   showfliers=False)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, rotation=90, fontsize="small")
        ax.set_title(quartile)
        ax.set_xlabel("Resampling percentage")
        ax.set_ylabel("Percent Relative Error")
    return save(fig, path)


# --------------------------------------------------------------------------
# Gene body coverage
# --------------------------------------------------------------------------


def _sample_colors(n: int):
    palette = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a",
               "#66a61e", "#e6ab02", "#a6761d", "#666666"]
    if n <= len(palette):
        return palette[:n]
    cmap = LinearSegmentedColormap.from_list("rseqc_samples", palette)
    return [cmap(i / max(1, n - 1)) for i in range(n)]


def _draw_curves(ax, profiles: Mapping[str, np.ndarray], colors, title):
    x = np.arange(1, 101)
    if not profiles:
        ax.set_title(f"{title} (no data)")
        ax.axis("off")
        return
    ymax = max(1.2, max(float(np.nanmax(p)) for p in profiles.values()) * 1.05)
    ax.axhline(1, ls=":", color="grey")
    for name, profile in profiles.items():
        ax.plot(x, profile, color=colors[name], lw=1.2, label=name)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Gene body percentile (5'->3')")
    ax.set_ylabel("Coverage (relative to mean)")
    ax.set_title(title)


def genebody_curves(
    profiles: Mapping[str, np.ndarray],
    by_length: Mapping[str, Mapping[str, np.ndarray]],
    counts: Mapping[str, Mapping[str, int]],
    prefix: str,
    ext: str = "pdf",
    *,
    min_transcripts: int = 100,
) -> list[Path]:
    """geneBody_coverage: curves, heatmap (>= 3 samples) and length panels.

    *profiles* are already scaled to their mean (1.0 = uniform).
    """
    names = list(profiles)
    colors = dict(zip(names, _sample_colors(len(names))))
    written = []

    fig, ax = _figure(9 if len(names) > 6 else 7, 7)
    _draw_curves(ax, profiles, colors, "Gene body coverage")
    if names:
        if len(names) > 6:
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize="small")
        else:
            ax.legend(loc="lower center", ncol=min(3, len(names)), fontsize="small", frameon=False)
    written.append(save(fig, f"{prefix}.curves.{ext}"))

    if len(names) >= 3:
        matrix = np.vstack([profiles[n] for n in names])
        fig, ax = _figure(8, max(3, 0.35 * len(names) + 2))
        cmap = LinearSegmentedColormap.from_list("rseqc_cm", ["cyan", "white", "magenta"])
        image = ax.imshow(matrix, aspect="auto", cmap=cmap, interpolation="nearest",
                          extent=(0.5, 100.5, len(names) - 0.5, -0.5))
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names)
        ax.set_xticks([1] + list(range(10, 101, 10)))
        ax.set_xlabel("Gene body percentile (5'->3')")
        fig.colorbar(image, ax=ax, label="Coverage (relative to mean)")
        written.append(save(fig, f"{prefix}.heatMap.{ext}"))

    if len(by_length) > 1:
        n_panels = len(by_length)
        n_col = n_panels if n_panels <= 2 else 2
        n_row = math.ceil(n_panels / n_col)
        fig, axes = plt.subplots(n_row, n_col, figsize=(5 * n_col, 4.5 * n_row), squeeze=False)
        for index, (ax, (label, panel)) in enumerate(zip(axes.flat, by_length.items())):
            n = counts.get(label, {})
            title = label
            if n:
                lo, hi = min(n.values()), max(n.values())
                title = f"{label}  (n = {lo if lo == hi else f'{lo}-{hi}'} transcripts)"
                if lo < min_transcripts:
                    title += " [low n]"
            _draw_curves(ax, panel, colors, title)
            if index == 0 and 0 < len(names) <= 6:
                ax.legend(loc="lower center", ncol=min(3, len(names)), fontsize="small", frameon=False)
        for ax in list(axes.flat)[n_panels:]:
            ax.axis("off")
        written.append(save(fig, f"{prefix}.byLength.{ext}"))
    return written


def genebody_step(profile: Sequence[float], n_transcripts: int, path, ylabel: str) -> Path:
    """geneBody_coverage2: average signal per transcript along the gene body."""
    fig, ax = _figure()
    y = np.asarray(profile, dtype=float) / max(n_transcripts, 1)
    ax.step(np.arange(1, len(y) + 1), y, where="post", color="black")
    ax.set_xlabel("percentile of gene body (5'->3')")
    ax.set_ylabel(ylabel)
    return save(fig, path)


# --------------------------------------------------------------------------
# Matrix heatmaps (sc_editMatrix, sc_seqQual)
# --------------------------------------------------------------------------


def matrix_heatmap(
    row_names: Sequence[str],
    col_names: Sequence[str],
    matrix,
    path,
    *,
    cell_width: float = 10,
    cell_height: float = 10,
    col_angle: int = 270,
    font_size: float = 10,
    text_color: str = "black",
    show_numbers: bool = True,
    log2_scale: bool = False,
) -> Path:
    """Heatmap of a labelled matrix (replaces the R pheatmap call).

    Cell sizes are in points, as in pheatmap.
    """
    values = np.asarray(matrix, dtype=float)
    shown = np.log2(values + 1) if log2_scale else values
    n_rows, n_cols = shown.shape
    width = max(3.0, n_cols * cell_width / 72 + 2.5)
    height = max(2.5, n_rows * cell_height / 72 + 2.0)
    fig, ax = plt.subplots(figsize=(width, height))
    cmap = LinearSegmentedColormap.from_list("rseqc_heat", ["#f1eef6", "#d7b5d8", "#df65b0", "#ce1256"], N=50)
    image = ax.imshow(shown, aspect="auto", cmap=cmap, interpolation="nearest")
    ax.set_xticks(range(n_cols))
    rotation = {270: 90, 315: 45}.get(col_angle, col_angle)
    ax.set_xticklabels(col_names, rotation=rotation, fontsize=font_size)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_names, fontsize=font_size)
    ax.tick_params(length=0)
    if show_numbers:
        for r in range(n_rows):
            for c in range(n_cols):
                value = shown[r, c]
                text = f"{value:.2f}" if log2_scale else f"{value:g}"
                ax.text(c, r, text, ha="center", va="center", fontsize=font_size * 0.8, color=text_color)
    fig.colorbar(image, ax=ax)
    return save(fig, path)


def run(parser, draw) -> None:
    """Call *draw* (a function that writes plots); report errors through the
    tool's argument *parser* and exit 1, as R failures did before."""
    try:
        draw()
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"{parser.prog}: error: plotting failed: {exc}\n")
