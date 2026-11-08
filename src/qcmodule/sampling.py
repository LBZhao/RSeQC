# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Sample alignments from across the genome instead of the start of a file.

Taking the first N reads of a coordinate-sorted BAM file samples only the
first few genes of the first chromosome. :func:`sample_alignments` instead
visits annotated regions in random order (fixed by *seed*) and takes at most
*per_region* reads from each, so no single locus, highly expressed gene or
chromosome dominates the sample.
"""

from __future__ import annotations

import random
import sys
from typing import Callable, Iterable, Iterator

DEFAULT_PER_REGION = 1000


def merge_regions(
    regions: Iterable[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    """Merge overlapping ``(chrom, start, end)`` regions."""
    merged: list[list] = []
    for chrom, start, end in sorted(regions):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    return [(chrom, start, end) for chrom, start, end in merged]


def sample_alignments(
    samfile,
    regions: Iterable[tuple[str, int, int]],
    limit: int,
    accept: Callable[[object], bool],
    *,
    per_region: int = DEFAULT_PER_REGION,
    seed: int | None = 0,
    warn=sys.stderr,
) -> Iterator:
    """Yield up to *limit* accepted alignments spread over *regions*.

    *accept* is called on every record; only records for which it returns
    True count towards the limits. Regions are merged, shuffled with *seed*
    and visited in turn, taking at most *per_region* records from each. A
    record overlapping two regions is yielded once.

    Files without an index (e.g. SAM) cannot be sampled by region; they are
    read from the start, as in earlier releases, with a warning.
    """
    if limit <= 0:
        return

    try:
        has_index = samfile.has_index()
    except (AttributeError, ValueError):
        has_index = False

    if not has_index:
        if warn is not None:
            print(
                "Warning: alignment file has no index; sampling the first "
                f"{limit} usable reads instead of reads across the genome.",
                file=warn,
            )
        count = 0
        for read in samfile:
            if accept(read):
                yield read
                count += 1
                if count >= limit:
                    return
        return

    references = set(samfile.references)
    regions = [r for r in merge_regions(regions) if r[0] in references]
    random.Random(seed).shuffle(regions)

    seen: set[tuple[str, int, int]] = set()
    count = 0
    for chrom, start, end in regions:
        taken = 0
        for read in samfile.fetch(chrom, start, end):
            if not accept(read):
                continue
            key = (read.query_name, read.flag, read.reference_start)
            if key in seen:
                continue
            seen.add(key)
            yield read
            count += 1
            taken += 1
            if count >= limit:
                return
            if taken >= per_region:
                break
