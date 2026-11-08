# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Core routines for gene-body coverage analysis.

Algorithm overview
------------------
1. Each BED12 transcript is reduced to its exonic blocks (the mature mRNA).
2. Overlapping transcripts are grouped into clusters, and every cluster's
   region is read from the BAM file once. Per-base coverage is built from
   the aligned blocks of each read (introns, deletions, soft clips, and
   insertions do not add coverage). When the two mates of a pair overlap,
   the shared bases are counted once.
3. Each transcript's exonic coverage is split into 100 equal slices along
   the mRNA and the mean coverage of every slice is taken. Minus-strand
   profiles are reversed so slice 1 is always the 5' end.
4. Profiles are aggregated across transcripts, overall and by mRNA length.

For BigWig input (geneBody_coverage2.py) step 2 is replaced by reading the
per-base signal over each exon block; steps 3 and 4 are the same.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence, TextIO

import numpy as np
import pysam

N_BINS = 100
DEFAULT_LENGTH_EDGES: tuple[int, ...] = (1000, 2000, 4000)

# Unmapped (0x4), secondary (0x100), QC fail (0x200), duplicate (0x400),
# supplementary (0x800).
EXCLUDED_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800

# Flush accumulated block coordinates into the coverage array this often so
# memory stays bounded inside very deeply covered clusters.
_FLUSH_EVERY = 1_000_000


# --------------------------------------------------------------------------
# Gene model
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Transcript:
    """A BED12 transcript reduced to its exonic blocks (0-based, half-open)."""

    tid: str
    chrom: str
    start: int
    end: int
    strand: str
    blocks: tuple[tuple[int, int], ...]
    length: int


@dataclass
class Cluster:
    """A genomic region holding one or more overlapping transcripts."""

    chrom: str
    start: int
    end: int
    transcripts: list[Transcript] = field(default_factory=list)


@dataclass
class BedSummary:
    transcripts: list[Transcript]
    n_records: int
    n_malformed: int
    n_duplicate: int
    n_short: int


def _parse_int_list(text: str) -> list[int]:
    return [int(value) for value in text.strip().rstrip(",").split(",")]


def parse_bed12(
    path: Path,
    min_length: int = 100,
    *,
    warn: TextIO | None = sys.stderr,
    max_warnings: int = 10,
) -> BedSummary:
    """Read BED12 transcripts whose exonic length is at least *min_length*.

    Transcripts with identical chrom, start, end, name and strand are kept
    once, as in earlier releases.
    """
    transcripts: list[Transcript] = []
    seen: set[str] = set()
    n_records = n_malformed = n_duplicate = n_short = 0

    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue

            fields = line.split()
            try:
                chrom = fields[0]
                tx_start = int(fields[1])
                tx_end = int(fields[2])
                name = fields[3]
                strand = fields[5]
                sizes = _parse_int_list(fields[10])
                offsets = _parse_int_list(fields[11])
                if len(sizes) != len(offsets) or not sizes:
                    raise ValueError("block size/start count mismatch")
                if any(size <= 0 for size in sizes):
                    raise ValueError("non-positive block size")
            except (IndexError, ValueError):
                n_malformed += 1
                if warn is not None and n_malformed <= max_warnings:
                    print(
                        "[NOTE: input BED must be 12-column] skipped line: "
                        + line.rstrip("\n"),
                        file=warn,
                    )
                continue

            n_records += 1
            tid = f"{chrom}_{tx_start}_{tx_end}_{name}_{strand}"
            if tid in seen:
                n_duplicate += 1
                continue
            seen.add(tid)

            blocks = tuple(
                sorted(
                    (tx_start + offset, tx_start + offset + size)
                    for offset, size in zip(offsets, sizes)
                )
            )
            length = sum(end - start for start, end in blocks)
            if length < min_length:
                n_short += 1
                continue

            transcripts.append(
                Transcript(
                    tid=tid,
                    chrom=chrom,
                    start=blocks[0][0],
                    end=blocks[-1][1],
                    strand=strand,
                    blocks=blocks,
                    length=length,
                )
            )

    if warn is not None and n_malformed > max_warnings:
        print(
            f"[NOTE] {n_malformed} malformed BED lines skipped in total",
            file=warn,
        )

    return BedSummary(transcripts, n_records, n_malformed, n_duplicate, n_short)


def build_clusters(transcripts: Iterable[Transcript]) -> list[Cluster]:
    """Group transcripts whose genomic spans overlap into clusters."""
    clusters: list[Cluster] = []
    ordered = sorted(transcripts, key=lambda tx: (tx.chrom, tx.start, tx.end))

    for tx in ordered:
        last = clusters[-1] if clusters else None
        if last is not None and last.chrom == tx.chrom and tx.start < last.end:
            last.end = max(last.end, tx.end)
            last.transcripts.append(tx)
        else:
            clusters.append(Cluster(tx.chrom, tx.start, tx.end, [tx]))

    return clusters


# --------------------------------------------------------------------------
# Length bins
# --------------------------------------------------------------------------


def format_length(bp: int) -> str:
    if bp >= 1000:
        return f"{bp / 1000:g}kb"
    return f"{bp}bp"


def length_bin_labels(edges: Sequence[int]) -> list[str]:
    """Return labels for bins defined by ascending *edges* (half-open)."""
    if not edges:
        return ["all"]
    labels = [f"<{format_length(edges[0])}"]
    labels += [
        f"{format_length(low)}-{format_length(high)}"
        for low, high in zip(edges[:-1], edges[1:])
    ]
    labels.append(f">={format_length(edges[-1])}")
    return labels


def parse_length_edges(text: str) -> list[int] | str:
    """Parse ``--length-bins``: comma-separated bp edges, 'auto' or 'none'."""
    value = text.strip().lower()
    if value in ("auto", "none"):
        return value
    try:
        edges = [int(float(item)) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise ValueError(f"invalid length bin edges: {text!r}") from exc
    if not edges:
        raise ValueError("at least one length bin edge is required")
    if any(edge <= 0 for edge in edges):
        raise ValueError("length bin edges must be positive")
    if any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError("length bin edges must be strictly increasing")
    return edges


def auto_length_edges(lengths: Sequence[int]) -> list[int]:
    """Return quartile edges of transcript lengths (rounded to 10 bp)."""
    if not lengths:
        return []
    quartiles = np.quantile(np.asarray(lengths, dtype=float), [0.25, 0.5, 0.75])
    edges = sorted({int(round(q / 10.0) * 10) for q in quartiles})
    return [edge for edge in edges if edge > 0]


def assign_length_bin(length: int, edges: Sequence[int]) -> int:
    """Index of the half-open length bin holding *length*."""
    return int(np.searchsorted(edges, length, side="right"))


# --------------------------------------------------------------------------
# Coverage
# --------------------------------------------------------------------------


def _merge_blocks(blocks: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(blocks):
        if merged and start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1][1] = end
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def region_coverage(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    *,
    min_mapq: int = 30,
) -> np.ndarray:
    """Per-base read coverage over ``[start, end)`` (0-based).

    Reads that are unmapped, secondary, supplementary, QC-failed, duplicates,
    or below *min_mapq* are ignored. Only aligned bases (CIGAR M/=/X) count.
    Bases shared by overlapping mates of the same pair are counted once.
    """
    span = end - start
    diff = np.zeros(span + 1, dtype=np.int64)
    starts: list[int] = []
    ends: list[int] = []
    pending: dict[str, list[tuple[int, int]]] = {}

    def flush() -> None:
        if starts:
            diff[:] += np.bincount(
                np.asarray(starts, dtype=np.int64), minlength=span + 1
            )
            diff[:] -= np.bincount(
                np.asarray(ends, dtype=np.int64), minlength=span + 1
            )
            starts.clear()
            ends.clear()

    def add(blocks: Iterable[tuple[int, int]]) -> None:
        for block_start, block_end in blocks:
            if block_start < start:
                block_start = start
            if block_end > end:
                block_end = end
            if block_start < block_end:
                starts.append(block_start - start)
                ends.append(block_end - start)
        if len(starts) >= _FLUSH_EVERY:
            flush()

    for read in bam.fetch(chrom, start, end):
        if read.flag & EXCLUDED_FLAGS or read.mapping_quality < min_mapq:
            continue
        blocks = read.get_blocks()
        if not blocks:
            continue

        if (
            read.is_paired
            and not read.mate_is_unmapped
            and read.next_reference_id == read.reference_id
        ):
            name = read.query_name
            mate_blocks = pending.pop(name, None)
            if mate_blocks is not None:
                add(_merge_blocks(mate_blocks + blocks))
                continue
            # The mate comes later in sorted order and may overlap this read.
            if read.reference_start <= read.next_reference_start < read.reference_end:
                pending[name] = blocks
                continue

        add(blocks)

    # Mates that never arrived (filtered out or outside the region).
    for blocks in pending.values():
        add(blocks)
    flush()

    return np.cumsum(diff[:span])


def slice_profile(
    exonic: np.ndarray,
    strand: str,
    n_bins: int = N_BINS,
    *,
    tid: str = "transcript",
) -> tuple[np.ndarray, float]:
    """Return the 5'->3' slice-mean profile and mean of per-base *exonic* values.

    *exonic* holds one value per mRNA base, in genomic (left to right) order.
    It is split into *n_bins* slices whose sizes differ by at most one base,
    and the mean of every slice is taken. Minus-strand profiles are reversed
    so slice 1 is always the 5' end.
    """
    length = exonic.size
    if length < n_bins:
        raise ValueError(f"transcript {tid} is shorter than {n_bins} bp")

    edges = (np.arange(n_bins + 1, dtype=np.int64) * length) // n_bins
    cumulative = np.concatenate(([0.0], np.cumsum(exonic, dtype=np.float64)))
    profile = (cumulative[edges[1:]] - cumulative[edges[:-1]]) / np.diff(edges)

    if strand == "-":
        profile = profile[::-1]

    return profile, float(cumulative[-1] / length)


def transcript_profile(
    coverage: np.ndarray,
    region_start: int,
    transcript: Transcript,
    n_bins: int = N_BINS,
) -> tuple[np.ndarray, float]:
    """Return the 5'->3' slice-mean profile and mean coverage of a transcript.

    *coverage* holds per-base coverage of a region beginning at
    *region_start* that contains the whole transcript.
    """
    exonic = np.concatenate(
        [
            coverage[block_start - region_start : block_end - region_start]
            for block_start, block_end in transcript.blocks
        ]
    )
    return slice_profile(exonic, transcript.strand, n_bins, tid=transcript.tid)


def bigwig_exonic_values(bigwig, transcript: Transcript) -> np.ndarray:
    """Per-base BigWig signal over a transcript's exons, in genomic order.

    *bigwig* is an open ``pyBigWig`` file. Bases without data (NaN) count
    as 0. Intronic positions are never read.
    """
    import pyBigWig

    kwargs = {"numpy": True} if pyBigWig.numpy else {}
    parts = [
        np.asarray(
            bigwig.values(transcript.chrom, block_start, block_end, **kwargs),
            dtype=np.float64,
        )
        for block_start, block_end in transcript.blocks
    ]
    return np.nan_to_num(np.concatenate(parts), nan=0.0)


def bigwig_profile(
    bigwig,
    transcript: Transcript,
    n_bins: int = N_BINS,
) -> tuple[np.ndarray, float]:
    """Return the 5'->3' slice-mean profile and mean signal from a BigWig."""
    return slice_profile(
        bigwig_exonic_values(bigwig, transcript),
        transcript.strand,
        n_bins,
        tid=transcript.tid,
    )


# --------------------------------------------------------------------------
# Per-sample aggregation
# --------------------------------------------------------------------------


@dataclass
class ProfileSum:
    values: np.ndarray = field(default_factory=lambda: np.zeros(N_BINS))
    n: int = 0

    def add(self, profile: np.ndarray) -> None:
        self.values += profile
        self.n += 1

    def result(self, transcript_norm: str) -> np.ndarray:
        """Summed profile (``none``) or mean of normalized profiles (``mean``)."""
        if transcript_norm == "mean":
            return self.values / self.n if self.n else np.zeros(N_BINS)
        return self.values.copy()


@dataclass
class SampleResult:
    name: str
    bam: str
    overall: ProfileSum
    by_length: list[ProfileSum]
    n_low_coverage: int = 0
    n_missing_chrom: int = 0
    n_out_of_range: int = 0


def compute_sample(
    bam_path: Path | str,
    name: str,
    clusters: Sequence[Cluster],
    *,
    length_edges: Sequence[int],
    transcript_norm: str = "mean",
    min_coverage: float = 5.0,
    min_mapq: int = 30,
    per_transcript: Path | None = None,
    progress: TextIO | None = None,
) -> SampleResult:
    """Compute overall and length-binned gene-body profiles for one BAM."""
    if transcript_norm not in ("mean", "none"):
        raise ValueError(f"unknown transcript normalization: {transcript_norm}")

    n_length_bins = len(length_edges) + 1
    result = SampleResult(
        name=name,
        bam=str(bam_path),
        overall=ProfileSum(),
        by_length=[ProfileSum() for _ in range(n_length_bins)],
    )
    labels = length_bin_labels(length_edges) if length_edges else ["all"]

    detail: TextIO | None = None
    if per_transcript is not None:
        detail = Path(per_transcript).open("w", encoding="utf-8")
        print(
            "Transcript\tLength\tLengthBin\tMeanCoverage\t"
            + "\t".join(str(i) for i in range(1, N_BINS + 1)),
            file=detail,
        )

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            references = set(bam.references)
            done = 0
            for cluster in clusters:
                if cluster.chrom not in references:
                    result.n_missing_chrom += len(cluster.transcripts)
                    continue

                coverage = region_coverage(
                    bam, cluster.chrom, cluster.start, cluster.end,
                    min_mapq=min_mapq,
                )

                for tx in cluster.transcripts:
                    profile, mean_cov = transcript_profile(
                        coverage, cluster.start, tx
                    )
                    if mean_cov <= 0 or mean_cov < min_coverage:
                        result.n_low_coverage += 1
                        continue

                    if transcript_norm == "mean":
                        profile = profile / profile.mean()

                    bin_index = (
                        assign_length_bin(tx.length, length_edges)
                        if length_edges
                        else 0
                    )
                    result.overall.add(profile)
                    result.by_length[bin_index].add(profile)

                    if detail is not None:
                        print(
                            f"{tx.tid}\t{tx.length}\t{labels[bin_index]}\t"
                            f"{mean_cov:.6g}\t"
                            + "\t".join(f"{value:.6g}" for value in profile),
                            file=detail,
                        )

                done += 1
                if progress is not None and done % 1000 == 0:
                    print(
                        f"\t{done} of {len(clusters)} gene clusters finished\r",
                        end="",
                        file=progress,
                    )
    finally:
        if detail is not None:
            detail.close()

    if progress is not None:
        print(file=progress)

    return result


def compute_bigwig_sample(
    bigwig_path: Path | str,
    name: str,
    transcripts: Sequence[Transcript],
    *,
    length_edges: Sequence[int],
    transcript_norm: str = "none",
    min_coverage: float = 0.0,
    progress: TextIO | None = None,
) -> SampleResult:
    """Compute overall and length-binned gene-body profiles from a BigWig.

    Same slicing as :func:`compute_sample`, but the per-base values are read
    from the BigWig signal over each transcript's exons (NaN counts as 0).
    With ``transcript_norm="none"`` (the default here, as in earlier
    releases) raw signal is summed across transcripts. Transcripts on
    chromosomes absent from the BigWig, or extending past the chromosome end,
    are skipped and counted.
    """
    import pyBigWig

    if transcript_norm not in ("mean", "none"):
        raise ValueError(f"unknown transcript normalization: {transcript_norm}")

    result = SampleResult(
        name=name,
        bam=str(bigwig_path),
        overall=ProfileSum(),
        by_length=[ProfileSum() for _ in range(len(length_edges) + 1)],
    )

    bigwig = pyBigWig.open(str(bigwig_path), "r")
    if bigwig is None:
        raise OSError(f"could not open BigWig file: {bigwig_path}")

    try:
        chrom_sizes = bigwig.chroms()
        for done, tx in enumerate(transcripts, start=1):
            size = chrom_sizes.get(tx.chrom)
            if size is None:
                result.n_missing_chrom += 1
                continue
            if tx.end > size:
                result.n_out_of_range += 1
                continue

            profile, mean_signal = bigwig_profile(bigwig, tx)
            if mean_signal <= 0 and transcript_norm == "mean":
                result.n_low_coverage += 1
                continue
            if mean_signal < min_coverage:
                result.n_low_coverage += 1
                continue

            if transcript_norm == "mean":
                profile = profile / profile.mean()

            bin_index = (
                assign_length_bin(tx.length, length_edges) if length_edges else 0
            )
            result.overall.add(profile)
            result.by_length[bin_index].add(profile)

            if progress is not None and done % 1000 == 0:
                print(
                    f"\t{done} of {len(transcripts)} transcripts finished\r",
                    end="",
                    file=progress,
                )
    finally:
        bigwig.close()

    if progress is not None:
        print(file=progress)

    return result


_WORKER_CLUSTERS: list[Cluster] = []


def init_worker(clusters: list[Cluster]) -> None:
    """Process-pool initializer: share the gene model with a worker once."""
    global _WORKER_CLUSTERS
    _WORKER_CLUSTERS = clusters


def run_worker_job(job: dict) -> SampleResult:
    """Process-pool task: compute one sample with the shared gene model."""
    return compute_sample(clusters=_WORKER_CLUSTERS, **job)


# --------------------------------------------------------------------------
# Summary metrics
# --------------------------------------------------------------------------


def _safe_ratio(numerator: float, denominator: float) -> float:
    """numerator / denominator; inf for x/0 with x > 0, NaN for 0/0."""
    if not (math.isfinite(numerator) and math.isfinite(denominator)):
        return math.nan
    if denominator == 0:
        return math.inf if numerator > 0 else math.nan
    return float(numerator / denominator)


def ratio_3to5(profile: Sequence[float]) -> float:
    """Mean coverage of slices 81-95 divided by mean of slices 6-20.

    1 means even coverage, > 1 means 3' bias and < 1 means 5' bias. The
    outermost five slices at each end are excluded because coverage there
    always falls off.
    """
    values = np.asarray(profile, dtype=float)
    return _safe_ratio(values[80:95].mean(), values[5:20].mean())


def uniformity_cv(profile: Sequence[float]) -> float:
    """Coefficient of variation over slices 10-90 (lower is more uniform)."""
    values = np.asarray(profile, dtype=float)[9:90]
    return _safe_ratio(values.std(ddof=1), values.mean())


def pearson_moment_coefficient(values: Sequence[float]) -> float:
    """Legacy RSeQC "skewness" score (deprecated).

    This describes the spread of coverage values around the 51st value and
    does not depend on position along the gene body, so it cannot tell a 5'
    bias from a 3' bias. Kept for backward compatibility only; use
    :func:`ratio_3to5` instead.
    """
    values = list(values)
    middle_value = values[int(len(values) / 2)]
    sigma = np.std(values, ddof=1)
    if sigma == 0 or not math.isfinite(sigma):
        return math.nan
    return float(
        np.mean([((value - middle_value) / sigma) ** 3 for value in values])
    )


def bias_rank_key(profile: Sequence[float]) -> float:
    """Sort key: distance of the 3'/5' ratio from 1 on a log scale."""
    ratio = ratio_3to5(profile)
    if math.isnan(ratio):
        return math.inf
    if ratio == 0 or math.isinf(ratio):
        return math.inf
    return abs(math.log2(ratio))


# --------------------------------------------------------------------------
# Shared table writers
# --------------------------------------------------------------------------


def fmt(value: float) -> str:
    """Format a number for output tables (NA / Inf for non-finite values)."""
    if math.isnan(value):
        return "NA"
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    return f"{value:.6g}"


def write_length_profiles(
    path: Path,
    results: Sequence[SampleResult],
    labels: Sequence[str],
    transcript_norm: str,
) -> None:
    with path.open("w", encoding="utf-8") as output:
        print(
            "Sample\tLengthBin\tnTranscripts\t"
            + "\t".join(str(i) for i in range(1, 101)),
            file=output,
        )
        for result in results:
            for label, profile_sum in zip(labels, result.by_length):
                profile = profile_sum.result(transcript_norm)
                values = (
                    "\t".join(fmt(v) for v in profile)
                    if profile_sum.n
                    else "\t".join(["NA"] * 100)
                )
                print(
                    f"{result.name}\t{label}\t{profile_sum.n}\t{values}",
                    file=output,
                )


def write_metrics(
    path: Path,
    results: Sequence[SampleResult],
    labels: Sequence[str],
    transcript_norm: str,
) -> None:
    with path.open("w", encoding="utf-8") as output:
        print(
            "Sample\tLengthBin\tnTranscripts\tRatio3to5\tUniformityCV\t"
            "Skewness_deprecated",
            file=output,
        )
        for result in results:
            groups = [("all", result.overall)] + list(
                zip(labels, result.by_length)
            )
            for label, profile_sum in groups:
                if profile_sum.n:
                    profile = profile_sum.result(transcript_norm)
                    ratio = ratio_3to5(profile)
                    cv = uniformity_cv(profile)
                    skew = pearson_moment_coefficient(profile)
                else:
                    ratio = cv = skew = math.nan
                print(
                    f"{result.name}\t{label}\t{profile_sum.n}\t"
                    f"{fmt(ratio)}\t{fmt(cv)}\t{fmt(skew)}",
                    file=output,
                )
