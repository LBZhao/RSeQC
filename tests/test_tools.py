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
