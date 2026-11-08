"""Tests for qcmodule.genebody (gene-body coverage core)."""

from __future__ import annotations

import math

import numpy as np
import pysam
import pytest

from conftest import CHROM, bed12_line
from qcmodule import genebody

# Three-exon transcript, 1000 bp of mRNA, with two large introns.
PLUS_BLOCKS = [(1000, 1300), (5000, 5400), (9000, 9300)]
# Same structure on the minus strand.
MINUS_BLOCKS = [(21000, 21300), (25000, 25400), (29000, 29300)]


def write_bed(tmp_path, records):
    path = tmp_path / "genes.bed"
    path.write_text("".join(bed12_line(*record) for record in records))
    return path


def single_profile(bam_path, blocks, strand="+", **kwargs):
    tx = genebody.Transcript(
        tid="tx", chrom=CHROM, start=blocks[0][0], end=blocks[-1][1],
        strand=strand, blocks=tuple(blocks),
        length=sum(e - s for s, e in blocks),
    )
    with pysam.AlignmentFile(str(bam_path)) as bam:
        coverage = genebody.region_coverage(
            bam, CHROM, tx.start, tx.end, **kwargs
        )
    return genebody.transcript_profile(coverage, tx.start, tx)


# --------------------------------------------------------------------------
# BED parsing, clustering and length bins
# --------------------------------------------------------------------------


def test_parse_bed12_lengths_duplicates_and_short(tmp_path):
    bed = write_bed(tmp_path, [
        ("a", "+", PLUS_BLOCKS),
        ("a", "+", PLUS_BLOCKS),            # exact duplicate
        ("short", "+", [(50000, 50099)]),    # 99 bp
        ("ok100", "+", [(60000, 60100)]),    # exactly 100 bp
    ])
    with bed.open("a") as handle:
        handle.write("chr1\t1\t2\tbroken\n")
    summary = genebody.parse_bed12(bed, 100, warn=None)
    lengths = {tx.tid.split("_")[3]: tx.length for tx in summary.transcripts}
    assert lengths == {"a": 1000, "ok100": 100}
    assert summary.n_duplicate == 1
    assert summary.n_short == 1
    assert summary.n_malformed == 1


def test_build_clusters_merges_overlaps_only():
    def tx(name, start, end):
        return genebody.Transcript(name, CHROM, start, end, "+",
                                   ((start, end),), end - start)
    clusters = genebody.build_clusters(
        [tx("c", 5000, 6000), tx("a", 0, 1000), tx("b", 500, 2000)]
    )
    assert [(c.start, c.end, len(c.transcripts)) for c in clusters] == [
        (0, 2000, 2), (5000, 6000, 1),
    ]


def test_length_bins_are_half_open_and_labelled():
    edges = [1000, 2000, 4000]
    assert genebody.length_bin_labels(edges) == [
        "<1kb", "1kb-2kb", "2kb-4kb", ">=4kb",
    ]
    assert [genebody.assign_length_bin(n, edges)
            for n in (999, 1000, 1999, 2000, 3999, 4000, 10**6)] == [
        0, 1, 1, 2, 2, 3, 3,
    ]
    assert genebody.length_bin_labels([500, 1500]) == [
        "<500bp", "500bp-1.5kb", ">=1.5kb",
    ]


@pytest.mark.parametrize("text", ["2000,1000", "0,100", "abc", ","])
def test_parse_length_edges_rejects_bad_input(text):
    with pytest.raises(ValueError):
        genebody.parse_length_edges(text)


def test_parse_length_edges_keywords_and_auto():
    assert genebody.parse_length_edges("auto") == "auto"
    assert genebody.parse_length_edges("none") == "none"
    assert genebody.parse_length_edges("500,1000") == [500, 1000]
    assert genebody.auto_length_edges(list(range(100, 1100))) == [350, 600, 850]


# --------------------------------------------------------------------------
# Coverage and profiles
# --------------------------------------------------------------------------


def test_uniform_coverage_gives_flat_profile(bam_builder):
    builder = bam_builder()
    builder.tile(PLUS_BLOCKS, depth=4)
    profile, mean = single_profile(builder.write(), PLUS_BLOCKS)
    assert profile.shape == (100,)
    np.testing.assert_allclose(profile, 4.0)
    assert mean == pytest.approx(4.0)
    assert genebody.ratio_3to5(profile) == pytest.approx(1.0)


def test_spliced_reads_add_no_intron_coverage(bam_builder):
    builder = bam_builder()
    builder.tile(PLUS_BLOCKS, depth=3)
    bam = builder.write()
    with pysam.AlignmentFile(str(bam)) as handle:
        coverage = genebody.region_coverage(handle, CHROM, 1000, 9300)
    assert coverage[1300 - 1000:5000 - 1000].sum() == 0
    assert coverage[5400 - 1000:9000 - 1000].sum() == 0
    assert (coverage[0:300] == 3).all()


def test_intronic_reads_do_not_count(bam_builder):
    builder = bam_builder()
    for start in range(1400, 4900, 50):     # reads inside intron 1 only
        builder.add(start, "50M")
    profile, mean = single_profile(builder.write(), PLUS_BLOCKS)
    assert mean == 0
    assert not profile.any()


def test_three_prime_bias_plus_and_minus_strand(bam_builder):
    builder = bam_builder()
    # Plus strand: 3' end is the genomic right end (exonic bases 700-1000).
    builder.tile(PLUS_BLOCKS, depth=2, first=700, last=1000)
    # Minus strand: 3' end is the genomic left end (exonic bases 0-300).
    builder.tile(MINUS_BLOCKS, depth=2, first=0, last=300)
    bam = builder.write()

    plus, _ = single_profile(bam, PLUS_BLOCKS, "+")
    minus, _ = single_profile(bam, MINUS_BLOCKS, "-")
    for profile in (plus, minus):
        assert not profile[:70].any()
        np.testing.assert_allclose(profile[70:], 2.0)
        assert genebody.ratio_3to5(profile) == math.inf
    np.testing.assert_allclose(plus, minus)


def test_transcript_of_exactly_100_bp_gives_100_values(bam_builder):
    blocks = [(60000, 60040), (60500, 60560)]
    builder = bam_builder()
    builder.add_exonic(blocks, 0, 100)
    profile, mean = single_profile(builder.write(), blocks)
    assert profile.shape == (100,)
    np.testing.assert_allclose(profile, 1.0)


def test_slices_average_over_every_base(bam_builder):
    builder = bam_builder()
    # Cover only the first half of slice 1 (bases 0-4 of a 1000 bp mRNA).
    builder.add_exonic(PLUS_BLOCKS, 0, 5)
    profile, _ = single_profile(builder.write(), PLUS_BLOCKS)
    assert profile[0] == pytest.approx(0.5)
    assert not profile[1:].any()


def test_read_filters(bam_builder):
    builder = bam_builder()
    builder.add(1000, "100M")                      # counted
    builder.add(1000, "100M", flag=0x400)          # duplicate
    builder.add(1000, "100M", flag=0x100)          # secondary
    builder.add(1000, "100M", flag=0x800)          # supplementary
    builder.add(1000, "100M", flag=0x200)          # QC fail
    builder.add(1000, "100M", mapq=10)             # low MAPQ
    bam = builder.write()
    with pysam.AlignmentFile(str(bam)) as handle:
        strict = genebody.region_coverage(handle, CHROM, 1000, 1100)
        lenient = genebody.region_coverage(handle, CHROM, 1000, 1100,
                                           min_mapq=0)
    assert (strict == 1).all()
    assert (lenient == 2).all()


def test_deletions_are_not_counted(bam_builder):
    builder = bam_builder()
    builder.add(1000, "20M10D20M")
    bam = builder.write()
    with pysam.AlignmentFile(str(bam)) as handle:
        coverage = genebody.region_coverage(handle, CHROM, 1000, 1050)
    assert coverage[:20].tolist() == [1] * 20
    assert coverage[20:30].tolist() == [0] * 10
    assert coverage[30:50].tolist() == [1] * 20


def test_overlapping_mates_counted_once(bam_builder):
    builder = bam_builder()
    builder.add_pair(1000, "60M", 1030, "60M")     # overlap 1030-1060
    builder.add_pair(1000, "50M", 1000, "50M")     # complete overlap
    builder.add_pair(2000, "50M", 2100, "50M")     # no overlap
    bam = builder.write()
    with pysam.AlignmentFile(str(bam)) as handle:
        coverage = genebody.region_coverage(handle, CHROM, 1000, 2200)
    assert coverage[0:50].tolist() == [2] * 50     # one per pair
    assert coverage[50:90].tolist() == [1] * 40
    assert coverage[1000:1050].tolist() == [1] * 50
    assert coverage[1100:1150].tolist() == [1] * 50
    assert coverage[1050:1100].sum() == 0


def test_mate_filtered_out_still_counts_other_mate(bam_builder):
    builder = bam_builder()
    name = "pair"
    builder.add(1000, "50M", flag=0x1 | 0x40, name=name, mate=1020)
    builder.add(1020, "50M", flag=0x1 | 0x80, name=name, mate=1000, mapq=0)
    bam = builder.write()
    with pysam.AlignmentFile(str(bam)) as handle:
        coverage = genebody.region_coverage(handle, CHROM, 1000, 1100)
    assert coverage[:50].tolist() == [1] * 50
    assert coverage[50:].sum() == 0


# --------------------------------------------------------------------------
# Sample aggregation
# --------------------------------------------------------------------------


def make_two_gene_sample(bam_builder, tmp_path):
    """A deep, flat gene and a shallow 3'-biased gene of different lengths."""
    long_blocks = [(40000, 42000), (45000, 47500)]    # 4500 bp
    bed = write_bed(tmp_path, [
        ("flat", "+", PLUS_BLOCKS),
        ("biased", "+", long_blocks),
    ])
    builder = bam_builder()
    builder.tile(PLUS_BLOCKS, depth=100)
    builder.tile(long_blocks, depth=10, first=2250, last=4500)
    bam = builder.write()
    clusters = genebody.build_clusters(
        genebody.parse_bed12(bed, 100, warn=None).transcripts
    )
    return bam, clusters


def test_mean_normalization_weights_transcripts_equally(bam_builder, tmp_path):
    bam, clusters = make_two_gene_sample(bam_builder, tmp_path)
    edges = [1000, 2000, 4000]
    mean = genebody.compute_sample(
        bam, "s", clusters, length_edges=edges, transcript_norm="mean",
        min_coverage=1,
    )
    raw = genebody.compute_sample(
        bam, "s", clusters, length_edges=edges, transcript_norm="none",
        min_coverage=0,
    )
    assert mean.overall.n == raw.overall.n == 2

    mean_profile = mean.overall.result("mean")
    # Average of flat (1.0) and biased (0 then 2.0) profiles.
    np.testing.assert_allclose(mean_profile[:50], 0.5)
    np.testing.assert_allclose(mean_profile[50:], 1.5)

    raw_profile = raw.overall.result("none")
    np.testing.assert_allclose(raw_profile[:50], 100.0)
    np.testing.assert_allclose(raw_profile[50:], 110.0)

    # Length bins: 1000 bp goes to 1-2kb, 4500 bp goes to >=4kb.
    assert [s.n for s in mean.by_length] == [0, 1, 0, 1]
    np.testing.assert_allclose(mean.by_length[1].result("mean"), 1.0)
    assert genebody.ratio_3to5(mean.by_length[3].result("mean")) == math.inf


def test_min_coverage_filter(bam_builder, tmp_path):
    bam, clusters = make_two_gene_sample(bam_builder, tmp_path)
    result = genebody.compute_sample(
        bam, "s", clusters, length_edges=[], transcript_norm="mean",
        min_coverage=50,
    )
    assert result.overall.n == 1
    assert result.n_low_coverage == 1


def test_missing_chromosome_is_reported(bam_builder, tmp_path):
    bed = write_bed(tmp_path, [("x", "+", [(100, 400)], "chrUn")])
    builder = bam_builder()
    builder.add(1000, "50M")
    clusters = genebody.build_clusters(
        genebody.parse_bed12(bed, 100, warn=None).transcripts
    )
    result = genebody.compute_sample(
        builder.write(), "s", clusters, length_edges=[], min_coverage=0,
    )
    assert result.n_missing_chrom == 1
    assert result.overall.n == 0


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------


def test_ratio_and_uniformity():
    ramp = np.linspace(0.5, 1.5, 100)
    assert genebody.ratio_3to5(ramp) > 1
    assert genebody.ratio_3to5(ramp[::-1]) < 1
    assert genebody.uniformity_cv(np.ones(100)) == 0
    assert genebody.uniformity_cv(ramp) > 0
    assert math.isnan(genebody.ratio_3to5(np.zeros(100)))


def test_legacy_skewness_cannot_see_direction():
    """Documents why the old score was replaced."""
    x = np.linspace(0, 1, 100)
    three_prime_tail = np.where(x > 0.8, 3.0, 1.0)
    five_prime_tail = three_prime_tail[::-1]
    assert genebody.pearson_moment_coefficient(three_prime_tail) == \
        pytest.approx(genebody.pearson_moment_coefficient(five_prime_tail))
    # The replacement metric does tell the two apart.
    assert genebody.ratio_3to5(three_prime_tail) > 1
    assert genebody.ratio_3to5(five_prime_tail) < 1
    assert math.isnan(genebody.pearson_moment_coefficient(np.ones(100)))
