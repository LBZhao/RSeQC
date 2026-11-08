"""Tests for qcmodule.BED.ParseBED."""

from __future__ import annotations

import pytest

from conftest import bed12_line
from qcmodule import BED


@pytest.fixture
def bed_with_headers(tmp_path):
    line = bed12_line("g1", "+", [(1000, 1500), (4500, 5000)]).split("\t")
    line[6], line[7] = "1200", "4800"          # CDS 1200-4800
    path = tmp_path / "genes.bed"
    path.write_text(
        "browser position chr1:1-10000\n"
        "track name=genes description='UCSC download'\n"
        "# comment\n\n" + "\t".join(line)
    )
    return path


@pytest.mark.parametrize("method, expected", [
    ("getCDSExon", [["chr1", 1200, 1500], ["chr1", 4500, 4800]]),
    ("getExon", [("chr1", 1000, 1500), ("chr1", 4500, 5000)]),
    ("getUTR", [["chr1", 1000, 1200, "g1", "0", "+"],
                ["chr1", 4800, 5000, "g1", "0", "+"]]),
    ("getIntron", [["chr1", 1500, 4500, "g1", "CDS"]]),
])
def test_header_lines_are_ignored(bed_with_headers, method, expected):
    bed = BED.ParseBED(str(bed_with_headers))
    assert getattr(bed, method)() == expected
    # Methods rewind, so a second call gives the same answer.
    assert getattr(bed, method)() == expected
