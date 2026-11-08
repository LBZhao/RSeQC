"""Regression tests for the second review (splitting, BigWig, k-mer and
single-cell tools)."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pysam
import pytest

from conftest import CHROM, ROOT, bed12_line
from qcmodule import FrameKmer, fastq, scbam


def run_script(name, *args, check=True, cwd=None):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / name), *map(str, args)],
        capture_output=True, text=True, env=env, check=check, cwd=cwd,
    )


def names(bam_path):
    with pysam.AlignmentFile(str(bam_path)) as handle:
        return [(r.query_name, r.flag) for r in handle.fetch(until_eof=True)]


# --------------------------------------------------------------------------
# split_bam / split_paired_bam / divide_bam
# --------------------------------------------------------------------------


def test_split_bam_looks_up_mate_on_its_own_chromosome(tmp_path, bam_builder):
    bed = tmp_path / "rrna.bed"
    bed.write_text(bed12_line("rRNA", "+", [(1000, 1300)]))
    builder = bam_builder()
    builder.header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": CHROM, "LN": 200_000}, {"SN": "chr2", "LN": 200_000}],
    })
    # Read on chr1 far from the exon; its mate is on chr2 at 1100. The old
    # code looked 1100 up on chr1, i.e. inside the exon.
    read = builder.add(50_000, "50M", flag=0x1 | 0x40, name="pair", mate=1100)
    read.next_reference_id = 1
    builder.add(1100, "50M", name="in_exon")
    bam = builder.write()

    prefix = tmp_path / "out"
    run_script("split_bam.py", "-i", bam, "-r", bed, "-o", prefix)
    assert [n for n, _ in names(f"{prefix}.in.bam")] == ["in_exon"]
    assert [n for n, _ in names(f"{prefix}.ex.bam")] == ["pair"]


def test_split_paired_bam_keeps_supplementary_flag(tmp_path, bam_builder):
    builder = bam_builder()
    builder.add(1000, "50M", flag=0x1 | 0x40, name="r", mate=2000)
    builder.add(1500, "20H30M", flag=0x1 | 0x40 | 0x800, name="r", mate=2000)
    bam = builder.write()
    prefix = tmp_path / "split"
    run_script("split_paired_bam.py", "-i", bam, "-o", prefix)
    flags = sorted(flag for _, flag in names(f"{prefix}.R1.bam"))
    assert flags == [0, 0x800]


def test_divide_bam_keeps_reads_together_and_is_reproducible(tmp_path, bam_builder):
    builder = bam_builder()
    for i in range(200):
        builder.add_pair(1000 + i, "50M", 1300 + i, "50M")
    bam = builder.write()

    def divide(tag):
        prefix = tmp_path / tag
        run_script("divide_bam.py", "-i", bam, "-n", "3", "-o", prefix, "--seed", "7")
        return [sorted(names(f"{prefix}_{i}.bam")) for i in range(3)]

    first, second = divide("a"), divide("b")
    assert first == second                                  # same seed, same split
    for subset in first:
        counts = {}
        for name, _ in subset:
            counts[name] = counts.get(name, 0) + 1
        assert set(counts.values()) == {2}                  # both mates together
    assert all(30 < len(subset) / 2 < 110 for subset in first)


# --------------------------------------------------------------------------
# overlay_bigwig
# --------------------------------------------------------------------------


def write_bigwig(path, entries, length=10_000):
    import pyBigWig
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(CHROM, length)])
    bw.addEntries([CHROM] * len(entries), [s for s, _, _ in entries],
                  ends=[e for _, e, _ in entries], values=[float(v) for _, _, v in entries])
    bw.close()
    return path


def test_overlay_division_skips_positions_without_data(tmp_path):
    pytest.importorskip("pyBigWig")
    a = write_bigwig(tmp_path / "a.bw", [(100, 110, 4.0)])
    b = write_bigwig(tmp_path / "b.bw", [(105, 115, 1.0)])
    out = tmp_path / "ratio.wig"
    run_script("overlay_bigwig.py", "-i", a, "-j", b, "-a", "Division", "-o", out, "-c", "1000")
    values = {
        int(pos): float(v)
        for pos, v in (line.split("\t") for line in out.read_text().splitlines()
                       if not line.startswith("variableStep"))
    }
    # Only positions 101-115 (1-based) have data in either file.
    assert sorted(values) == list(range(101, 116))
    assert values[101] == pytest.approx(5.0)        # (4 + 1) / (0 + 1)
    assert values[106] == pytest.approx(2.5)        # (4 + 1) / (1 + 1)
    assert values[115] == pytest.approx(0.5)        # (0 + 1) / (1 + 1)


# --------------------------------------------------------------------------
# read_hexamer parsing
# --------------------------------------------------------------------------


def test_fastq_quality_lines_are_not_sequence(tmp_path):
    fq = tmp_path / "reads.fq"
    # Quality "CCCCGGGG" is valid Phred+33 but looks like DNA; a quality
    # line may also start with '@'. The empty read is valid FASTQ.
    fq.write_text(
        "@r1\nAAAAAAAA\n+\nCCCCGGGG\n"
        "@r2\nTTTTTTTT\n+\n@@@@@@@@\n"
        "@r3\n\n+\n\n"
        "@r4\nAAAAAAAA\n+\nIIIIIIII\n"
    )
    records = list(FrameKmer.seq_generator(str(fq)))
    assert records == [["r1", "AAAAAAAA"], ["r2", "TTTTTTTT"], ["r3", ""], ["r4", "AAAAAAAA"]]
    counts = FrameKmer.kmer_freq_file(str(fq), word_size=6)
    assert counts["AAAAAA"] == 6 and counts["TTTTTT"] == 3
    assert counts["CCCCGG"] == 0


def test_fasta_iupac_codes_do_not_create_false_kmers(tmp_path):
    fa = tmp_path / "genome.fa"
    fa.write_text(">chr\nAAAAAA\nCCRCCC\nGGGGGG\n")
    assert list(FrameKmer.seq_generator(str(fa))) == [["chr", "AAAAAACCNCCCGGGGGG"]]
    counts = FrameKmer.kmer_freq_file(str(fa), word_size=6)
    # The line with R used to be dropped, creating AAAAAAGGGGGG junctions.
    assert counts["AAAAGG"] == 0
    assert counts["AAAAAA"] == 1


# --------------------------------------------------------------------------
# FPKM-UQ
# --------------------------------------------------------------------------


def load_fpkm_uq():
    spec = importlib.util.spec_from_file_location("fpkm_uq", ROOT / "scripts" / "FPKM-UQ.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fpkm_uq_refuses_zero_upper_quartile(tmp_path):
    fpkm_uq = load_fpkm_uq()
    info = tmp_path / "info.txt"
    info.write_text("".join(
        f"g{i}\tG{i}\tchr1\t1\t100\t+\tprotein_coding\tx\tx\tx\t1000\n" for i in range(8)
    ))
    counts = tmp_path / "counts.txt"
    counts.write_text("".join(f"g{i}\t{10 if i == 0 else 0}\n" for i in range(8)))
    with pytest.raises(ValueError, match="75th percentile"):
        fpkm_uq.calculate_fpkm(counts, info, tmp_path / "out.txt")


# --------------------------------------------------------------------------
# Single-cell helpers
# --------------------------------------------------------------------------


def test_strip_suffix_handles_any_gem_well():
    assert scbam.strip_suffix("ACGTACGT-1") == "ACGTACGT"
    assert scbam.strip_suffix("ACGTACGT-12") == "ACGTACGT"
    assert scbam.strip_suffix("ACGTACGT") == "ACGTACGT"


def test_barcode_edits_aggregated_library(tmp_path, bam_builder):
    builder = bam_builder()
    read = builder.add(1000, "50M", name="a")
    read.set_tags([("CR", "ACGTACGT"), ("CB", "ACGTACGT-2"), ("UR", "AAAA"), ("UB", "AAAA")])
    read = builder.add(1100, "50M", name="b")
    read.set_tags([("CR", "ACGTACGA"), ("CB", "ACGTACGT-2"), ("UR", "AAAA"), ("UB", "AAAA")])
    read = builder.add(1200, "50M", name="b", flag=0x100)      # secondary: ignored
    read.set_tags([("CR", "ACGTACGA"), ("CB", "ACGTACGT-2"), ("UR", "AAAA"), ("UB", "AAAA")])
    bam = builder.write()
    scbam.barcode_edits(str(bam), str(tmp_path / "e"))
    freq = (tmp_path / "e.CB_freq.tsv").read_text().split()
    assert freq == ["ACGTACGT", "2"]
    edits = (tmp_path / "e.CB_edits_count.csv").read_text()
    assert "A:T" in edits                      # one real edit (original:corrected)


def test_mapping_stat_leaves_other_files_alone(tmp_path, bam_builder, monkeypatch):
    builder = bam_builder()
    read = builder.add(1000, "50M", name="a")
    read.set_tags([("xf", 1), ("CB", "ACGT-1")])
    bam = builder.write()
    workdir = tmp_path / "work"
    workdir.mkdir()
    unrelated = workdir / "chr9.all_reads_id.txt"
    unrelated.write_text("someone_else\n")
    monkeypatch.chdir(workdir)
    scbam.mapping_stat(str(bam))
    assert unrelated.read_text() == "someone_else\n"
    assert (workdir / "All_reads_uniqID.txt").read_text().split() == ["a"]


def test_fasta_iter_joins_wrapped_sequences(tmp_path):
    fa = tmp_path / "barcodes.fa"
    fa.write_text(">a\nACGT\nTTTT\n>b\nGGGG\n")
    assert list(fastq.fasta_iter(str(fa))) == ["ACGTTTTT", "GGGG"]


# --------------------------------------------------------------------------
# RNA_fragment_size
# --------------------------------------------------------------------------


def test_fragment_size_across_an_intron(tmp_path, bam_builder):
    bed = tmp_path / "g.bed"
    bed.write_text(bed12_line("g", "+", [(1000, 1300), (5000, 5400)]))
    builder = bam_builder()
    for _ in range(5):
        builder.add_pair(1200, "50M", 5050, "50M")    # 100 + 100 exonic bases
        builder.add_pair(1100, "50M", 1200, "50M")    # 150 bp in one exon
    bam = builder.write()
    out = tmp_path / "frag.txt"
    run_script("RNA_fragment_size.py", "-i", bam, "-r", bed, "-o", out, "-n", "1")
    row = out.read_text().splitlines()[1].split("\t")
    assert float(row[5]) == pytest.approx(175.0)      # mean of 200 and 150
    assert np.median([200, 150] * 5) == float(row[6])
