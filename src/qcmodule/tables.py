# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

"""Write machine-readable tab-separated tables (``--tsv``)."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Sequence


def _cell(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "Inf" if value > 0 else "-Inf"
        return f"{value:.6g}"
    return str(value)


def write_tsv(path, header: Sequence[str], rows: Iterable[Sequence]) -> Path:
    """Write *rows* under *header* as tab-separated text; None/NaN -> NA."""
    path = Path(path)
    with path.open("w", encoding="utf-8") as handle:
        print("\t".join(header), file=handle)
        for row in rows:
            print("\t".join(_cell(value) for value in row), file=handle)
    return path
