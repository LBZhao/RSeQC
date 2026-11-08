#!/usr/bin/env python3
"""Utilities for generating and analyzing DNA k-mers.

DNA sequences are expected to contain A, C, G, T, and optionally N.
"""

from __future__ import annotations

import itertools
import math
import re
from collections import Counter


_NON_ACGT = re.compile(r"[^ACGT]")


def word_generator(seq, word_size, step_size, frame=0):
    """Yield fixed-length words from a sequence.

    Parameters
    ----------
    seq : str
        Input sequence.
    word_size : int
        Length of each word.
    step_size : int
        Distance between successive word starts.
    frame : int, default=0
        Starting offset in the sequence.
    """
    if word_size <= 0:
        raise ValueError("word_size must be greater than zero")
    if step_size <= 0:
        raise ValueError("step_size must be greater than zero")
    if frame < 0:
        raise ValueError("frame must be non-negative")

    for i in range(frame, len(seq), step_size):
        word = seq[i:i + word_size]
        if len(word) == word_size:
            yield word


def seq_generator(fastafile):
    """Yield ``[name, sequence]`` records from a FASTA or FASTQ file.

    The format is taken from the first record: ``>`` means FASTA (sequences
    may span several lines), ``@`` means FASTQ (four lines per record; the
    quality line is never read as sequence). Sequences are upper-cased and
    any character other than A, C, G or T becomes N, so k-mers across
    ambiguous bases are skipped rather than created by joining lines.

    Before RSeQC 5.1.0, FASTQ quality lines made only of the characters
    A/C/G/T/N (valid quality codes) were counted as sequence, quality lines
    starting with ``@`` started a new record, and FASTA lines with any IUPAC
    code were dropped, joining their neighbours into false k-mers.
    """
    with open(fastafile, "r") as handle:
        raw = (line.strip() for line in handle)
        first = next((line for line in raw if line and not line.startswith("#")), None)
        if first is None:
            return
        if first.startswith("@"):
            # FASTQ: keep blank lines (zero-length reads are valid).
            yield from _fastq_records(first, raw, fastafile)
        else:
            lines = (line for line in raw if line and not line.startswith("#"))
            yield from _fasta_records(first, lines)


def _clean(seq):
    return _NON_ACGT.sub("N", seq.upper())


def _fasta_records(first, lines):
    name = first[1:].split()[0] if first.startswith(">") else ""
    parts = [] if first.startswith(">") else [_clean(first)]
    for line in lines:
        if line.startswith(">"):
            yield [name, "".join(parts)]
            name = line[1:].split()[0] if len(line) > 1 else ""
            parts = []
        else:
            parts.append(_clean(line))
    yield [name, "".join(parts)]


def _fastq_records(first, lines, path):
    header = first
    while header is not None:
        if not header:                      # trailing blank lines
            header = next(lines, None)
            continue
        seq = next(lines, None)
        plus = next(lines, None)
        qual = next(lines, None)
        if qual is None or not plus.startswith("+") or not header.startswith("@"):
            raise ValueError(f"{path}: malformed FASTQ record starting {header!r}")
        yield [header[1:].split()[0] if len(header) > 1 else "", _clean(seq)]
        header = next(lines, None)


def all_possible_kmer(l):
    """Yield all length-``l`` combinations of A, C, G, T, and N."""
    if l < 0:
        raise ValueError("k-mer length must be non-negative")

    for bases in itertools.product("ACGTN", repeat=l):
        yield "".join(bases)


def kmer_freq_file(fastafile, word_size, step_size=1, frame=0, min_count=0):
    """Return k-mer counts from sequences in a FASTA-like file.

    K-mers containing ``N`` are omitted from the returned dictionary.
    """
    count_table = Counter()

    for _name, seq in seq_generator(fastafile):
        count_table.update(
            word_generator(
                seq,
                word_size=word_size,
                step_size=step_size,
                frame=frame,
            )
        )

    ret_dict = {}
    for kmer in all_possible_kmer(word_size):
        count = count_table.get(kmer, 0)
        if count >= min_count and "N" not in kmer:
            ret_dict[kmer] = count

    return ret_dict


def kmer_freq_seq(seq, word_size, step_size=1, frame=0, min_count=0):
    """Print k-mer counts for a DNA sequence."""
    count_table = Counter(
        word_generator(
            seq,
            word_size=word_size,
            step_size=step_size,
            frame=frame,
        )
    )

    for kmer in all_possible_kmer(word_size):
        count = count_table.get(kmer, 0)
        if count >= min_count:
            print(f"{kmer}\t{count}")


def kmer_ratio(seq, word_size, step_size, coding, noncoding):
    """Calculate the frame-0 coding/noncoding k-mer log-ratio score.

    This preserves the active algorithm in the original implementation:
    only frame 0 is evaluated. K-mers absent from either reference mapping
    are skipped.

    Returns
    -------
    float
        Mean frame-0 score. Returns 0.0 if the sequence is shorter than
        ``word_size`` or if no usable k-mers are available.
    """
    if len(seq) < word_size:
        return 0.0

    sum_of_log_ratio_0 = 0.0
    frame0_count = 0

    for kmer in word_generator(
        seq=seq,
        word_size=word_size,
        step_size=step_size,
        frame=0,
    ):
        if kmer not in coding or kmer not in noncoding:
            continue

        coding_count = coding[kmer]
        noncoding_count = noncoding[kmer]

        if coding_count > 0 and noncoding_count > 0:
            sum_of_log_ratio_0 += math.log(coding_count / noncoding_count)
        elif coding_count > 0 and noncoding_count == 0:
            sum_of_log_ratio_0 += 1
        elif coding_count == 0 and noncoding_count == 0:
            continue
        elif coding_count == 0 and noncoding_count > 0:
            sum_of_log_ratio_0 -= 1
        else:
            continue

        frame0_count += 1

    if frame0_count == 0:
        return 0.0

    return sum_of_log_ratio_0 / frame0_count
