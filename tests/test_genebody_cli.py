"""End-to-end tests for scripts/geneBody_coverage.py."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import ROOT, bed12_line

SCRIPT = ROOT / "scripts" / "geneBody_coverage.py"

GENES = [
    ("short", "+", [(1000, 1300), (5000, 5400)]),                 # 700 bp
    ("mid", "-", [(20000, 20800), (22000, 22700)]),               # 1.5 kb
    ("long", "+", [(40000, 42000), (45000, 47500)]),              # 4.5 kb
]


def run_cli(*args, check=True):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        capture_output=True, text=True, env=env, check=check,
    )


@pytest.fixture
def dataset(tmp_path, bam_builder):
    bed = tmp_path / "genes.bed"
    bed.write_text("".join(bed12_line(*gene) for gene in GENES))
    bams = []
    for index, name in enumerate(["S1.bam", "S2.bam", "S3.bam"]):
        builder = bam_builder(name)
        for _, strand, blocks in GENES:
            total = sum(e - s for s, e in blocks)
            builder.tile(blocks, depth=10)
            # Extra 3' coverage that grows with the sample index.
            if strand == "+":
                builder.tile(blocks, depth=5 * index, first=total // 2)
            else:
                builder.tile(blocks, depth=5 * index, last=total // 2)
        bams.append(builder.write())
    return tmp_path, bed, bams


def read_table(path):
    rows = [line.rstrip("\n").split("\t") for line in path.open()]
    return rows[0], rows[1:]


def test_cli_outputs(dataset):
    tmp_path, bed, bams = dataset
    prefix = tmp_path / "out"
    run_cli("-i", ",".join(map(str, bams)), "-r", bed, "-o", prefix,
            "--skip-plot", "--per-transcript", "-p", "2")

    header, rows = read_table(Path(f"{prefix}.geneBodyCoverage.txt"))
    assert header == ["Percentile"] + [str(i) for i in range(1, 101)]
    assert [row[0] for row in rows] == ["S1", "S2", "S3"]
    assert all(len(row) == 101 for row in rows)

    header, rows = read_table(Path(f"{prefix}.geneBodyCoverage.byLength.txt"))
    assert header[:3] == ["Sample", "LengthBin", "nTranscripts"]
    s1 = {row[1]: row[2] for row in rows if row[0] == "S1"}
    assert s1 == {"<1kb": "1", "1kb-2kb": "1", "2kb-4kb": "0", ">=4kb": "1"}

    header, rows = read_table(Path(f"{prefix}.geneBodyCoverage.metrics.txt"))
    ratios = {row[0]: float(row[3]) for row in rows if row[1] == "all"}
    assert ratios["S1"] == pytest.approx(1.0)
    assert ratios["S1"] < ratios["S2"] < ratios["S3"]

    per_tx = Path(f"{prefix}.S1.geneBodyCoverage.perTranscript.txt")
    assert len(per_tx.read_text().splitlines()) == 1 + len(GENES)

    assert Path(f"{prefix}.geneBodyCoverage.r").is_file()
    assert Path(f"{prefix}.geneBodyCoverage.log.txt").is_file()


def test_cli_legacy_mode_and_no_bins(dataset):
    tmp_path, bed, bams = dataset
    prefix = tmp_path / "legacy"
    run_cli("-i", bams[0], "-r", bed, "-o", prefix, "--skip-plot",
            "--transcript-norm", "none", "--length-bins", "none")
    _, rows = read_table(Path(f"{prefix}.geneBodyCoverage.txt"))
    # Summed raw coverage of three flat genes at depth 10.
    assert all(float(value) == pytest.approx(30.0) for value in rows[0][1:])
    assert not Path(f"{prefix}.geneBodyCoverage.byLength.txt").exists()


def test_cli_rejects_bad_length_bins(dataset):
    tmp_path, bed, bams = dataset
    result = run_cli("-i", bams[0], "-r", bed, "-o", tmp_path / "x",
                     "--length-bins", "2000,1000", check=False)
    assert result.returncode != 0
    assert "--length-bins" in result.stderr


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="R not installed")
def test_cli_plots_with_r(dataset):
    tmp_path, bed, bams = dataset
    prefix = tmp_path / "plot"
    run_cli("-i", ",".join(map(str, bams)), "-r", bed, "-o", prefix,
            "-f", "png")
    for suffix in ("curves", "heatMap", "byLength"):
        image = Path(f"{prefix}.geneBodyCoverage.{suffix}.png")
        assert image.is_file() and image.stat().st_size > 0
