"""Tests for BigWig gene-body coverage (geneBody_coverage2.py)."""

from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pyBigWig = pytest.importorskip("pyBigWig")

from conftest import CHROM, ROOT, bed12_line  # noqa: E402
from qcmodule import genebody  # noqa: E402

SCRIPT = ROOT / "scripts" / "geneBody_coverage2.py"
BW_LENGTH = 100_000

# 150 bp exon, ~10 kb intron, 100 bp exon (250 bp mRNA). The pre-5.1.0
# percentile interpolation put percentile 60 at an intronic position here.
SPLICED_BLOCKS = [(1000, 1150), (11150, 11250)]
PLUS_BLOCKS = [(20000, 20300), (25000, 25400), (29000, 29300)]     # 1000 bp
MINUS_BLOCKS = [(40000, 40300), (45000, 45400), (49000, 49300)]    # 1000 bp


def write_bigwig(path, signal, chrom=CHROM):
    """Write a dense per-base *signal* array as a BigWig.

    NaN positions are left without data, as in a real BigWig.
    """
    signal = np.asarray(signal, dtype=float)
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(chrom, signal.size)])
    has_data = ~np.isnan(signal)
    # Runs of identical values among positions that have data.
    change = np.ones(signal.size, dtype=bool)
    change[1:] = (signal[1:] != signal[:-1]) | (has_data[1:] != has_data[:-1])
    starts = np.flatnonzero(change)
    ends = np.append(starts[1:], signal.size)
    keep = has_data[starts]
    starts, ends = starts[keep], ends[keep]
    if starts.size:
        bw.addEntries(
            [chrom] * starts.size,
            [int(s) for s in starts],
            ends=[int(e) for e in ends],
            values=[float(signal[s]) for s in starts],
        )
    bw.close()
    return Path(path)


def empty_signal():
    return np.full(BW_LENGTH, np.nan)


def make_tx(blocks, strand="+", name="tx"):
    return genebody.Transcript(
        tid=name, chrom=CHROM, start=blocks[0][0], end=blocks[-1][1],
        strand=strand, blocks=tuple(blocks),
        length=sum(e - s for s, e in blocks),
    )


def profile_of(bw_path, blocks, strand="+"):
    bw = pyBigWig.open(str(bw_path))
    try:
        return genebody.bigwig_profile(bw, make_tx(blocks, strand))
    finally:
        bw.close()


def paint(signal, blocks, value, first=0, last=None):
    """Set exonic bases [first, last) (from the genomic left end) to *value*."""
    total = sum(e - s for s, e in blocks)
    last = total if last is None else last
    position = 0
    for start, end in blocks:
        for genomic in range(start, end):
            if first <= position < last:
                signal[genomic] = value
            position += 1


# --------------------------------------------------------------------------
# Profiles
# --------------------------------------------------------------------------


def test_flat_signal_gives_flat_profile(tmp_path):
    bw = write_bigwig(tmp_path / "flat.bw", np.full(BW_LENGTH, 4.0))
    for blocks in (PLUS_BLOCKS, SPLICED_BLOCKS):
        profile, mean = profile_of(bw, blocks)
        assert profile.shape == (100,)
        np.testing.assert_allclose(profile, 4.0)
        assert mean == pytest.approx(4.0)
        assert genebody.ratio_3to5(profile) == pytest.approx(1.0)
        assert genebody.uniformity_cv(profile) == pytest.approx(0.0)


def test_intronic_signal_does_not_leak(tmp_path):
    signal = np.full(BW_LENGTH, 1000.0)          # strong signal everywhere ...
    paint(signal, SPLICED_BLOCKS, 1.0)           # ... except the exons
    bw = write_bigwig(tmp_path / "intron.bw", signal)
    profile, mean = profile_of(bw, SPLICED_BLOCKS)
    np.testing.assert_allclose(profile, 1.0)
    assert mean == pytest.approx(1.0)


def test_signal_only_in_intron_gives_zero_profile(tmp_path):
    signal = empty_signal()
    signal[1150:11150] = 50.0                    # intron of SPLICED_BLOCKS
    bw = write_bigwig(tmp_path / "intron_only.bw", signal)
    profile, mean = profile_of(bw, SPLICED_BLOCKS)
    assert mean == 0
    assert not profile.any()


def test_minus_strand_is_reversed(tmp_path):
    signal = np.zeros(BW_LENGTH)
    # Plus strand: 3' end is the genomic right end (exonic bases 700-1000).
    paint(signal, PLUS_BLOCKS, 2.0, first=700)
    # Minus strand: 3' end is the genomic left end (exonic bases 0-300).
    paint(signal, MINUS_BLOCKS, 2.0, last=300)
    bw = write_bigwig(tmp_path / "bias.bw", signal)

    plus, _ = profile_of(bw, PLUS_BLOCKS, "+")
    minus, _ = profile_of(bw, MINUS_BLOCKS, "-")
    for profile in (plus, minus):
        assert not profile[:70].any()
        np.testing.assert_allclose(profile[70:], 2.0)
        assert genebody.ratio_3to5(profile) == math.inf
    np.testing.assert_allclose(plus, minus)

    # The same minus-strand blocks read as plus strand give the mirror image.
    as_plus, _ = profile_of(bw, MINUS_BLOCKS, "+")
    np.testing.assert_allclose(as_plus, minus[::-1])


def test_missing_data_counts_as_zero(tmp_path):
    signal = empty_signal()
    paint(signal, PLUS_BLOCKS, 3.0, first=0, last=5)   # half of slice 1
    bw = write_bigwig(tmp_path / "sparse.bw", signal)
    profile, mean = profile_of(bw, PLUS_BLOCKS)
    assert profile[0] == pytest.approx(1.5)
    assert not profile[1:].any()
    assert mean == pytest.approx(15 / 1000)


def test_transcript_of_exactly_100_bp_gives_100_values(tmp_path):
    blocks = [(60000, 60040), (60500, 60560)]
    signal = np.zeros(BW_LENGTH)
    paint(signal, blocks, 1.0)
    signal[blocks[0][1]:blocks[1][0]] = 99.0
    bw = write_bigwig(tmp_path / "short.bw", signal)
    profile, _ = profile_of(bw, blocks)
    assert profile.shape == (100,)
    np.testing.assert_allclose(profile, 1.0)


# --------------------------------------------------------------------------
# Sample aggregation
# --------------------------------------------------------------------------


def test_compute_bigwig_sample_sums_and_bins(tmp_path):
    long_blocks = [(60000, 62000), (65000, 67500)]     # 4500 bp
    signal = np.zeros(BW_LENGTH)
    paint(signal, PLUS_BLOCKS, 100.0)
    paint(signal, long_blocks, 10.0, first=2250)
    bw = write_bigwig(tmp_path / "s.bw", signal)

    transcripts = [
        make_tx(PLUS_BLOCKS, name="flat"),
        make_tx(long_blocks, name="biased"),
        genebody.Transcript("other", "chrUn", 0, 200, "+", ((0, 200),), 200),
        make_tx([(BW_LENGTH - 50, BW_LENGTH + 100)], name="past_end"),
    ]
    edges = [1000, 2000, 4000]
    raw = genebody.compute_bigwig_sample(
        bw, "s", transcripts, length_edges=edges,
    )
    assert raw.overall.n == 2
    assert raw.n_missing_chrom == 1
    assert raw.n_out_of_range == 1
    summed = raw.overall.result("none")
    np.testing.assert_allclose(summed[:50], 100.0)
    np.testing.assert_allclose(summed[50:], 110.0)
    assert [s.n for s in raw.by_length] == [0, 1, 0, 1]

    mean = genebody.compute_bigwig_sample(
        bw, "s", transcripts, length_edges=edges, transcript_norm="mean",
        min_coverage=1,
    )
    profile = mean.overall.result("mean")
    np.testing.assert_allclose(profile[:50], 0.5)
    np.testing.assert_allclose(profile[50:], 1.5)

    filtered = genebody.compute_bigwig_sample(
        bw, "s", transcripts, length_edges=[], min_coverage=50,
    )
    assert filtered.overall.n == 1
    assert filtered.n_low_coverage == 1


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------


def run_cli(*args, check=True):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        capture_output=True, text=True, env=env, check=check,
    )


@pytest.fixture
def cli_dataset(tmp_path):
    genes = [
        ("spliced", "+", SPLICED_BLOCKS),
        ("plus", "+", PLUS_BLOCKS),
        ("minus", "-", MINUS_BLOCKS),
    ]
    bed = tmp_path / "genes.bed"
    bed.write_text("".join(bed12_line(*gene) for gene in genes))
    signal = np.full(BW_LENGTH, 500.0)             # intronic / intergenic noise
    for _, _, blocks in genes:
        paint(signal, blocks, 2.0)
    bw = write_bigwig(tmp_path / "sample.bw", signal)
    return tmp_path, bed, bw


def test_cli_output_layout_unchanged(cli_dataset):
    tmp_path, bed, bw = cli_dataset
    prefix = tmp_path / "out"
    run_cli("-i", bw, "-r", bed, "-o", prefix, "--skip-plot")

    lines = Path(f"{prefix}.geneBodyCoverage.txt").read_text().splitlines()
    assert lines[0] == "percentile\tcount"
    rows = [line.split("\t") for line in lines[1:]]
    assert [row[0] for row in rows] == [str(i) for i in range(100)]
    # Three flat transcripts at 2.0, summed; intronic 500s never leak in.
    assert all(float(row[1]) == pytest.approx(6.0) for row in rows)

    r_code = Path(f"{prefix}.geneBodyCoverage_plot.r").read_text()
    assert "x=1:100" in r_code and "y/3" in r_code

    metrics = [
        line.split("\t")
        for line in Path(f"{prefix}.geneBodyCoverage.metrics.txt")
        .read_text().splitlines()
    ]
    assert metrics[0][:5] == [
        "Sample", "LengthBin", "nTranscripts", "Ratio3to5", "UniformityCV",
    ]
    overall = metrics[1]
    assert overall[:3] == ["sample", "all", "3"]
    assert float(overall[3]) == pytest.approx(1.0)

    by_length = [
        line.split("\t")
        for line in Path(f"{prefix}.geneBodyCoverage.byLength.txt")
        .read_text().splitlines()
    ]
    assert {row[1]: row[2] for row in by_length[1:]} == {
        "<1kb": "1", "1kb-2kb": "2", "2kb-4kb": "0", ">=4kb": "0",
    }


def test_cli_no_length_bins(cli_dataset):
    tmp_path, bed, bw = cli_dataset
    prefix = tmp_path / "nobins"
    run_cli("-i", bw, "-r", bed, "-o", prefix, "--skip-plot",
            "--length-bins", "none")
    assert Path(f"{prefix}.geneBodyCoverage.txt").is_file()
    assert Path(f"{prefix}.geneBodyCoverage.metrics.txt").is_file()
    assert not Path(f"{prefix}.geneBodyCoverage.byLength.txt").exists()


def test_cli_reports_chromosome_mismatch(tmp_path):
    bed = tmp_path / "genes.bed"
    bed.write_text(bed12_line("x", "+", PLUS_BLOCKS, chrom="1"))
    bw = write_bigwig(tmp_path / "s.bw", np.ones(BW_LENGTH))
    result = run_cli("-i", bw, "-r", bed, "-o", tmp_path / "x",
                     "--skip-plot", check=False)
    assert result.returncode != 0
    assert "chromosome naming" in result.stderr


def test_cli_rejects_bad_length_bins(cli_dataset):
    tmp_path, bed, bw = cli_dataset
    result = run_cli("-i", bw, "-r", bed, "-o", tmp_path / "x",
                     "--length-bins", "2000,1000", check=False)
    assert result.returncode != 0
    assert "--length-bins" in result.stderr
