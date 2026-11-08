"""Tests for qcmodule.readfilter and its use by the command-line tools."""

from __future__ import annotations

import os
import subprocess
import sys

import pysam
import pytest

from conftest import ROOT
from qcmodule import readfilter

HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 10_000}]})


def make_read(flag=0, mapq=60):
    read = pysam.AlignedSegment(HEADER)
    read.reference_id = 0
    read.reference_start = 100
    read.cigarstring = "50M"
    read.flag = flag
    read.mapping_quality = mapq
    return read


@pytest.mark.parametrize("flag", [0x4, 0x100, 0x200, 0x400, 0x800])
def test_default_filter_rejects_each_flag(flag):
    assert readfilter.passes(make_read())
    assert not readfilter.passes(make_read(flag=flag))


def test_mapq_and_custom_exclusion():
    assert not readfilter.passes(make_read(mapq=10), min_mapq=30)
    assert readfilter.passes(make_read(mapq=30), min_mapq=30)
    # Keep duplicates (e.g. when measuring duplication).
    dup = make_read(flag=0x400)
    assert readfilter.passes(dup, exclude=readfilter.DEFAULT_EXCLUDE & ~readfilter.DUPLICATE)


def test_is_non_primary():
    assert not readfilter.is_non_primary(make_read())
    assert readfilter.is_non_primary(make_read(flag=0x100))
    assert readfilter.is_non_primary(make_read(flag=0x800))
    assert not readfilter.is_non_primary(make_read(flag=0x400))


# --------------------------------------------------------------------------
# End to end: supplementary alignments are not counted as reads
# --------------------------------------------------------------------------


def run_script(name, *args):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / name), *map(str, args)],
        capture_output=True, text=True, env=env, check=True,
    )


@pytest.fixture
def chimeric_bam(bam_builder):
    builder = bam_builder("chimeric.bam")
    builder.add(1000, "100M", name="plain")
    builder.add(2000, "60M40S", name="chimera")
    builder.add(5000, "60H40M", flag=0x800, name="chimera")    # supplementary
    builder.add(3000, "100M", flag=0x100, mapq=0, name="multi")
    return builder.write()


def test_bam_stat_counts_supplementary_separately(chimeric_bam):
    out = run_script("bam_stat.py", "-i", chimeric_bam).stdout
    stats = {
        line[:40].strip(): line[40:].strip()
        for line in out.splitlines()
        if len(line) > 40
    }
    assert stats["Total records:"] == "4"
    assert stats["Supplementary hits:"] == "1"
    assert stats["Non primary hits"] == "1"
    assert stats["mapq >= mapq_cut (unique):"] == "2"


def test_bam2fq_writes_each_read_once(chimeric_bam, tmp_path):
    prefix = tmp_path / "reads"
    run_script("bam2fq.py", "-i", chimeric_bam, "-o", prefix, "-s")
    lines = (tmp_path / "reads.fastq").read_text().splitlines()
    names = [line[1:] for line in lines[0::4]]
    assert sorted(names) == ["chimera", "plain"]


@pytest.mark.parametrize("rule, expected", [
    (None, {}),
    ("1++,1--,2+-,2-+", {"1+": "+", "1-": "-", "2+": "-", "2-": "+"}),
    ("+-,-+", {"+": "-", "-": "+"}),
])
def test_parse_strand_rule(rule, expected):
    assert readfilter.parse_strand_rule(rule) == expected


@pytest.mark.parametrize("rule", ["1++,1--", "++,--,+-", "1+x,1--,2+-,2-+", "abc"])
def test_parse_strand_rule_rejects_bad_rules(rule):
    with pytest.raises(ValueError):
        readfilter.parse_strand_rule(rule)


def test_gene_strand_for_unpaired_read_under_paired_rule():
    rule = readfilter.parse_strand_rule("1++,1--,2+-,2-+")
    assert readfilter.gene_strand(make_read(flag=0x1 | 0x40), rule) == "+"
    assert readfilter.gene_strand(make_read(flag=0x1 | 0x80 | 0x10), rule) == "+"
    assert readfilter.gene_strand(make_read(), rule) is None
