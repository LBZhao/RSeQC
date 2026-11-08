"""End-to-end and unit tests for the BAM-based QC tools."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import CHROM, ROOT, bed12_line
from qcmodule import SAM

# Three-exon gene: introns A = 1300-5000 and B = 5400-9000.
GENE = [(1000, 1300), (5000, 5400), (9000, 9300)]


def run_script(name, *args, check=True):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / name), *map(str, args)],
        capture_output=True, text=True, env=env, check=check,
    )


def write_bed(tmp_path, records, header=True, name="genes.bed"):
    path = tmp_path / name
    text = "track name=genes\n" if header else ""
    path.write_text(text + "".join(records))
    return path


def table(path):
    return [line.rstrip("\n").split("\t") for line in Path(path).open()]


# --------------------------------------------------------------------------
# junction_annotation
# --------------------------------------------------------------------------


def test_junction_classes(tmp_path, bam_builder):
    bed = write_bed(tmp_path, [bed12_line("g", "+", GENE)])
    builder = bam_builder()
    builder.add(1250, "50M3700N50M")      # intron A: annotated
    builder.add(1250, "50M3700N50M")
    builder.add(1250, "50M7700N50M")      # 1300-9000: exon skip, both sites known
    builder.add(1250, "50M3650N50M")      # 1300-4950: novel acceptor
    builder.add(1200, "50M3000N50M")      # 1250-4250: both sites novel
    builder.add(1250, "50M20N50M")        # 20 bp gap: filtered (< 50)
    bam = builder.write()

    prefix = tmp_path / "j"
    result = run_script("junction_annotation.py", "-i", bam, "-r", bed, "-o", prefix,
                        "--skip-plot", "--skip-bed", "--skip-interact")
    rows = {(r[1], r[2]): (r[0], r[3], r[4])
            for r in table(f"{prefix}.junction.xls")[1:]}
    assert rows == {
        ("1300", "5000"): (CHROM, "2", "annotated"),
        ("1300", "9000"): (CHROM, "1", "novel_combination"),
        ("1300", "4950"): (CHROM, "1", "partial_novel"),
        ("1250", "4250"): (CHROM, "1", "complete_novel"),
    }
    assert "Known Splicing Events:\t2" in result.stderr
    assert "Novel Combination Splicing Junctions:\t1" in result.stderr
    assert "Filtered Splicing Events:\t1" in result.stderr


def test_classify_junction_keeps_chromosome_case_insensitive(tmp_path):
    bed = write_bed(tmp_path, [bed12_line("g", "+", GENE, chrom="chrUn_KI270302v1")])
    reference = SAM.load_reference_introns(str(bed))
    assert SAM.classify_junction(reference, "chrUn_KI270302v1", 1300, 5000) == "annotated"
    assert SAM.classify_junction(reference, "CHRUN_KI270302V1", 1300, 5000) == "annotated"
    assert SAM.classify_junction(reference, "chr2", 1300, 5000) == "complete_novel"


def test_junction_output_keeps_bam_chromosome_names(tmp_path, bam_builder):
    import pysam
    name = "chrUn_KI270302v1"
    bed = write_bed(tmp_path, [bed12_line("g", "+", GENE, chrom=name)])
    builder = bam_builder()
    builder.header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": name, "LN": 200_000}]}
    )
    builder.add(1250, "50M3700N50M")
    bam = builder.write()
    prefix = tmp_path / "j"
    run_script("junction_annotation.py", "-i", bam, "-r", bed, "-o", prefix,
               "--skip-plot", "--skip-bed", "--skip-interact")
    assert table(f"{prefix}.junction.xls")[1][0] == name


# --------------------------------------------------------------------------
# Saturation tools
# --------------------------------------------------------------------------


@pytest.mark.parametrize("start, end, step, expected", [
    (5, 100, 5, list(range(5, 101, 5))),
    (10, 100, 5, list(range(10, 101, 5))),
    (5, 50, 5, [5, 10, 15, 20, 25, 30, 35, 40, 45, 100]),
    (5, 100, 30, [5, 35, 65, 95, 100]),
])
def test_sampling_percentages(start, end, step, expected):
    assert SAM.sampling_percentages(start, end, step) == expected


def make_spliced_sample(bam_builder, n_reads=40):
    builder = bam_builder()
    for _ in range(n_reads):
        builder.add(1250, "50M3700N50M")
    return builder.write()


def test_junction_saturation_uses_every_read_once(tmp_path, bam_builder):
    bed = write_bed(tmp_path, [bed12_line("g", "+", GENE)])
    bam = make_spliced_sample(bam_builder, 40)
    builder = bam_builder("novel.bam")
    for _ in range(40):
        builder.add(1250, "50M3700N50M")
    builder.add(1200, "50M3000N50M")          # one read on a novel junction
    bam = builder.write()

    prefix = tmp_path / "js"
    # Low bound 50 used to skip the first 45% of reads; upper bound 60 used to
    # skip 55-95%.
    result = run_script("junction_saturation.py", "-i", bam, "-r", bed, "-o", prefix,
                        "-l", "50", "-u", "60", "-s", "5", "--skip-plot")
    last = [line for line in result.stderr.splitlines() if line.startswith("sampling 100%")][0]
    assert "(41)" in last
    assert "2 splicing junctions" in last

    # The novel junction has one read: with -v 2 it is not counted anywhere.
    result = run_script("junction_saturation.py", "-i", bam, "-r", bed, "-o", prefix,
                        "-v", "2", "--skip-plot")
    last = [line for line in result.stderr.splitlines() if line.startswith("sampling 100%")][0]
    assert "1 splicing junctions. 1 known splicing junctions. 0 novel" in last


def test_rpkm_saturation_is_reproducible_and_complete(tmp_path, bam_builder):
    bed = write_bed(tmp_path, [bed12_line("g", "+", GENE, chrom=CHROM)])
    builder = bam_builder()
    builder.tile(GENE, depth=2)
    bam = builder.write()

    outputs = []
    for run in range(2):
        prefix = tmp_path / f"rs{run}"
        # Exits non-zero because one tiny gene fails the plot filter, but the
        # tables are written first.
        run_script("RPKM_saturation.py", "-i", bam, "-r", bed, "-o", prefix,
                   "-l", "10", "-u", "50", "--skip-plot", check=False)
        outputs.append(Path(f"{prefix}.rawCount.xls").read_text())
    assert outputs[0] == outputs[1]
    header, row = table(tmp_path / "rs0.rawCount.xls")
    assert row[0] == CHROM                       # name kept as in the BED file
    # 2 x 1000 bp of 50 bp reads = 40 blocks, all counted at 100%.
    assert header[-1] == "100%" and row[-1] == "40"


# --------------------------------------------------------------------------
# read_distribution
# --------------------------------------------------------------------------


def distribution_rows(stdout):
    rows = {}
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[1].isdigit():
            rows[parts[0]] = (int(parts[1]), int(parts[2]))
    return rows


def test_read_distribution_noncoding_mapq_and_boundaries(tmp_path, bam_builder):
    coding = bed12_line("coding", "+", [(1000, 1500), (4500, 5000)]).split("\t")
    coding[6], coding[7] = "1200", "4800"
    # Non-coding transcript, UCSC style: thickStart == thickEnd == txEnd.
    noncoding = bed12_line("lnc", "+", [(20000, 20400), (21000, 21400)]).split("\t")
    noncoding[6] = noncoding[7] = "21400"
    bed = write_bed(tmp_path, ["\t".join(coding), "\t".join(noncoding)])

    builder = bam_builder()
    builder.add(1300, "50M")              # CDS
    builder.add(1000, "2M")               # midpoint 1001: 5' UTR
    builder.add(20100, "50M")             # non-coding exon
    builder.add(21000, "1M")              # midpoint == exon start 21000
    builder.add(20100, "50M", mapq=5)     # low MAPQ: not counted
    bam = builder.write()

    out = run_script("read_distribution.py", "-i", bam, "-r", bed).stdout
    rows = distribution_rows(out)
    assert "Total Reads                   4" in out
    assert rows["CDS_Exons"][1] == 1
    assert rows["5'UTR_Exons"] == (200, 1)
    assert rows["3'UTR_Exons"] == (200, 0)
    assert rows["Noncoding_Exons"] == (800, 2)

    out = run_script("read_distribution.py", "-i", bam, "-r", bed, "-q", "0").stdout
    assert distribution_rows(out)["Noncoding_Exons"] == (800, 3)


# --------------------------------------------------------------------------
# read_NVC / read_quality
# --------------------------------------------------------------------------


def test_nvc_and_quality_cover_the_longest_read(tmp_path, bam_builder):
    builder = bam_builder()
    builder.add(1000, "80M")
    builder.add(2000, "30M")          # last read in the file is the shortest
    bam = builder.write()

    run_script("read_NVC.py", "-i", bam, "-o", tmp_path / "nvc", "--skip-plot")
    rows = table(tmp_path / "nvc.NVC.xls")
    assert len(rows) == 1 + 80
    # Positions 30-79 are covered by the long read only.
    assert [value.strip() for value in rows[-1][1:5]] == ["1", "0", "0", "0"]

    run_script("read_quality.py", "-i", bam, "-o", tmp_path / "q", "--skip-plot")
    r_code = (tmp_path / "q.qual.r").read_text()
    assert "p79<-" in r_code and "p80<-" not in r_code


# --------------------------------------------------------------------------
# infer_experiment / inner_distance
# --------------------------------------------------------------------------


def infer_fractions(stdout):
    values = {}
    for line in stdout.splitlines():
        if line.startswith("Fraction of reads explained by"):
            pattern = line.split('"')[1]
            values[pattern] = float(line.rsplit(":", 1)[1])
    return values


def test_infer_experiment_samples_across_genes(tmp_path, bam_builder):
    bed = write_bed(tmp_path, [
        bed12_line("first", "+", [(1000, 5000)]),
        bed12_line("second", "-", [(100000, 104000)]),
    ])
    builder = bam_builder()
    for i in range(3000):
        builder.add(1000 + i, "50M")       # sense reads on the first gene
        builder.add(100000 + i, "50M")     # antisense reads on the second gene
    bam = builder.write()

    # The first 2,000 reads of the sorted file all come from the first gene;
    # sampling by gene takes reads from both.
    out = run_script("infer_experiment.py", "-i", bam, "-r", bed, "-s", "2000").stdout
    assert "This is SingleEnd Data" in out
    fractions = infer_fractions(out)
    assert fractions["++,--"] == pytest.approx(0.5)
    assert fractions["+-,-+"] == pytest.approx(0.5)


def test_infer_experiment_uses_reference_end_for_spliced_reads(tmp_path, bam_builder):
    # A read that starts upstream of the gene and is spliced into its exon:
    # pos + query length ends before the gene, the reference end does not.
    bed = write_bed(tmp_path, [bed12_line("g", "-", [(5000, 6000)])])
    builder = bam_builder()
    for _ in range(10):
        builder.add(1000, "25M3990N25M")
    bam = builder.write()
    out = run_script("infer_experiment.py", "-i", bam, "-r", bed, "-s", "1000").stdout
    assert infer_fractions(out)["+-,-+"] == pytest.approx(1.0)


def test_inner_distance_on_mrna(tmp_path, bam_builder):
    bed = write_bed(tmp_path, [bed12_line("only", "+", GENE)])
    builder = bam_builder()
    builder.add_pair(1200, "50M", 5050, "50M")    # spans intron A
    builder.add_pair(1200, "50M", 1230, "50M")    # mates overlap by 20 bp
    builder.add_pair(1100, "20M5D30M", 1250, "50M")
    bam = builder.write()
    prefix = tmp_path / "id"
    result = run_script("inner_distance.py", "-i", bam, "-r", bed, "-o", prefix,
                        "--skip-plot", "-l", "-100", "-u", "100")
    rows = sorted(line.split("\t")[1:] for line in
                  Path(f"{prefix}.inner_distance.txt").read_text().splitlines())
    assert rows == sorted([
        ["100", "sameTranscript=Yes,sameExon=No,dist=mRNA"],   # 50 + 50 exonic bases
        ["-20", "readPairOverlap"],
        ["95", "sameTranscript=Yes,sameExon=Yes,dist=mRNA"],   # read end includes the deletion
    ])
    assert "> 100: 0" in result.stderr


# --------------------------------------------------------------------------
# tin
# --------------------------------------------------------------------------


def tin_scores(tmp_path, bam, bed, *args):
    out_dir = tmp_path / "tin_out"
    out_dir.mkdir(exist_ok=True)
    run_script("tin.py", "-i", bam, "-r", bed, "-o", out_dir, *args)
    rows = table(next(out_dir.glob("*.tin.xls")))
    return {row[0]: row[4] for row in rows[1:]}


def expected_tin(depths):
    import math
    total = sum(depths)
    entropy = -sum(d / total * math.log(d / total) for d in depths if d)
    return 100 * math.exp(entropy) / len(depths)


def test_tin_no_depth_cap_na_and_mapq(tmp_path, bam_builder):
    bed = write_bed(tmp_path, [
        bed12_line("deep", "+", [(1000, 1200)]),
        bed12_line("flat", "+", [(5000, 5100)]),
        bed12_line("empty", "+", [(9000, 9100)]),
    ])
    builder = bam_builder()
    # "deep": 9,000 reads on the first half, 4,500 on the second. pysam's
    # default pileup stopped at 8,000 reads per base.
    for _ in range(9000):
        builder.add(1000, "100M")
    for _ in range(4500):
        builder.add(1100, "100M")
    # "flat": even coverage.
    for _ in range(20):
        builder.add(5000, "100M")
    # Low-MAPQ reads would make "empty" pass the read check if counted.
    for start in range(9000, 9050):
        builder.add(start, "50M", mapq=0)
    bam = builder.write()

    scores = tin_scores(tmp_path, bam, bed, "-c", "0", "-n", "1000")
    assert float(scores["deep"]) == pytest.approx(expected_tin([9000] * 100 + [4500] * 100), rel=1e-6)
    # Every base of a 100 bp transcript sampled once: perfectly even = 100.
    assert float(scores["flat"]) == pytest.approx(100.0)
    assert scores["empty"] == "NA"


# --------------------------------------------------------------------------
# bam2wig
# --------------------------------------------------------------------------


def read_wig(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("variableStep"):
            continue
        pos, value = line.split("\t")
        values[int(pos)] = float(value)
    return values


def test_bamtowig_coverage_windows_and_normalization(tmp_path, bam_builder):
    builder = bam_builder()
    builder.add(995, "5S10M")             # crosses a window boundary at 1000
    builder.add(1003, "2M100N3M")
    builder.add(1003, "5M", mapq=0)       # skipped with skip_multi
    bam = builder.write()
    sizes = {CHROM: 200_000}
    chrom_file = tmp_path / "chrom.sizes"
    chrom_file.write_text(f"{CHROM}\t200000\n")

    alignment = SAM.ParseBAM(str(bam))
    total = alignment.calWigSum(sizes, skip_multi=True, q_cut=30)
    assert total == 15                     # 10 + 2 + 3 aligned bases

    prefix = str(tmp_path / "w")
    SAM.ParseBAM(str(bam)).bamTowig(prefix, sizes, str(chrom_file), skip_multi=True,
                                   q_cut=30, WigSumFactor=2.0, window=1000)
    wig = read_wig(prefix + ".wig")
    # 1-based wig positions; the soft clip does not shift the read.
    expected = {p: 2.0 for p in range(996, 1006)}
    expected.update({1004: 4.0, 1005: 4.0, 1106: 2.0, 1107: 2.0, 1108: 2.0})
    assert wig == expected


def test_bamtowig_strand_rule(tmp_path, bam_builder):
    builder = bam_builder()
    builder.add_pair(1000, "10M", 1100, "10M")   # read1 +, read2 -
    builder.add(2000, "10M")                     # unpaired: not covered by the rule
    bam = builder.write()
    prefix = str(tmp_path / "s")
    chrom_file = tmp_path / "chrom.sizes"
    chrom_file.write_text(f"{CHROM}\t200000\n")
    SAM.ParseBAM(str(bam)).bamTowig(prefix, {CHROM: 200_000}, str(chrom_file),
                                   strand_rule="1+-,1-+,2++,2--", skip_multi=False)
    # Rule 1+- : read 1 on '+' comes from a '-' gene; read 2 on '-' too.
    assert read_wig(prefix + ".Forward.wig") == {}
    reverse = read_wig(prefix + ".Reverse.wig")
    assert set(reverse) == set(range(1001, 1011)) | set(range(1101, 1111))
    assert set(reverse.values()) == {-1.0}


# --------------------------------------------------------------------------
# FPKM_count
# --------------------------------------------------------------------------


def fpkm_counts(tmp_path, bam, bed, *args):
    prefix = tmp_path / "fpkm"
    run_script("FPKM_count.py", "-i", bam, "-r", bed, "-o", prefix, *args)
    rows = [r for r in table(f"{prefix}.FPKM.xls") if not r[0].startswith("#")]
    return {row[3]: float(row[6]) for row in rows}


def test_fpkm_count_spliced_reads_and_pairs_at_transcript_edge(tmp_path, bam_builder):
    bed = write_bed(tmp_path, [bed12_line("g", "+", GENE)])
    builder = bam_builder()
    builder.add_pair(1100, "50M", 1200, "50M")      # inside the transcript
    # Read 1 starts upstream of the transcript and is not fetched for it;
    # read 2 starts in exon 1.
    builder.add_pair(900, "50M", 1100, "50M")
    # Single-end read spliced from upstream into exon 1: pos + read length
    # ends before the exon, the aligned blocks do not.
    builder.add(500, "10M600N40M")
    builder.add(1200, "50M", flag=0x800)            # supplementary: not counted
    bam = builder.write()
    assert fpkm_counts(tmp_path, bam, bed)["g"] == 3


def test_fpkm_count_stranded_pair_counted_through_read2(tmp_path, bam_builder):
    bed = write_bed(tmp_path, [bed12_line("g", "+", GENE)])
    builder = bam_builder()
    builder.add_pair(900, "50M", 1100, "50M")       # read 1 '+', read 2 '-'
    builder.add_pair(1100, "50M", 1200, "50M")      # fully exonic pair
    bam = builder.write()
    prefix = tmp_path / "s"
    run_script("FPKM_count.py", "-i", bam, "-r", bed, "-o", prefix, "-d", "1++,1--,2+-,2-+")
    header = [r for r in table(f"{prefix}.FPKM.xls") if r[0].startswith("#")][-1]
    row = [r for r in table(f"{prefix}.FPKM.xls") if not r[0].startswith("#")][0]
    values = dict(zip(header, row))
    assert float(values["Frag_count"]) == 2


# --------------------------------------------------------------------------
# read_duplication
# --------------------------------------------------------------------------


def test_positional_duplication_uses_mate_position(tmp_path, bam_builder):
    builder = bam_builder()
    builder.add_pair(1000, "50M", 1200, "50M")
    builder.add_pair(1000, "50M", 1300, "50M")     # same read-1 position, other mate
    builder.add_pair(2000, "50M", 2200, "50M")
    builder.add_pair(2000, "50M", 2200, "50M")     # true duplicate pair
    bam = builder.write()
    prefix = tmp_path / "dup"
    run_script("read_duplication.py", "-i", bam, "-o", prefix, "--skip-plot")
    occurrence = {int(r[0]): int(r[1]) for r in table(f"{prefix}.pos.DupRate.xls")[1:]}
    # Reads seen once: both 1000-pairs (4 reads); seen twice: the 2000 pairs.
    assert occurrence == {1: 4, 2: 2}
