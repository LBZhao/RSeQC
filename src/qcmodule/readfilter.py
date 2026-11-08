# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Shared rules for which BAM records count as reads.

A sequenced read can appear in a BAM file more than once:

* **secondary** alignments (0x100) are alternative placements of a
  multi-mapped read;
* **supplementary** alignments (0x800) are the extra pieces of a chimeric
  or split read (e.g. a fusion or circular RNA junction).

Counting either of these counts the same molecule twice, so every RSeQC tool
that counts reads, fragments or coverage skips them. Tools differ only in
whether duplicates (0x400) and QC failures (0x200) are also skipped and in
their mapping-quality cut-off.

Before RSeQC 5.1.0 supplementary alignments were counted by every tool
except ``geneBody_coverage.py``.
"""

from __future__ import annotations

UNMAPPED = 0x4
SECONDARY = 0x100
QCFAIL = 0x200
DUPLICATE = 0x400
SUPPLEMENTARY = 0x800

#: Records that are not the primary alignment of a read.
NON_PRIMARY = SECONDARY | SUPPLEMENTARY

#: Default exclusion: unmapped, secondary, supplementary, QC-failed and
#: duplicate records.
DEFAULT_EXCLUDE = UNMAPPED | NON_PRIMARY | QCFAIL | DUPLICATE


def passes(read, *, min_mapq: int = 0, exclude: int = DEFAULT_EXCLUDE) -> bool:
    """Return True if *read* (a ``pysam.AlignedSegment``) should be counted.

    The record is rejected if it has any flag in *exclude* or if its mapping
    quality is below *min_mapq*.
    """
    if read.flag & exclude:
        return False
    return read.mapping_quality >= min_mapq


def is_non_primary(read) -> bool:
    """True for secondary or supplementary records."""
    return bool(read.flag & NON_PRIMARY)
