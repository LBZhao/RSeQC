# MultiQC support for RSeQC 5.1.0 output

RSeQC 5.1.0 adds two report categories that MultiQC (tested with 1.35)
does not know yet, so it drops them from its report:

| Tool | New output | Effect in MultiQC without the patch |
|---|---|---|
| `junction_annotation.py` | `Novel Combination Splicing Events/Junctions` lines (both splice sites annotated but never joined, e.g. exon skipping; earlier releases called these "known") | Junction categories no longer add up to 100% |
| `read_distribution.py` | `Noncoding_Exons` row (exons of non-coding transcripts; earlier releases counted them as UTR) | Those tags are missing from the read distribution plot |

All other RSeQC sections in MultiQC parse unchanged.

`0001-rseqc-novel-combination-noncoding.patch` adds both categories to
MultiQC's RSeQC module. The fields are optional, so older RSeQC output is
parsed exactly as before.

Tested against MultiQC `main` at `a286ece`:

- `multiqc --strict` on RSeQC 5.1.0 output: junction categories sum to
  100% and `Noncoding_Exons` appears in the plot.
- `multiqc --strict` on RSeQC 5.0.6 output: report unchanged.
- `ruff check` / `ruff format`: no new findings.

## Submitting

Submit after an official RSeQC release writes these lines, so the pull
request can name the release:

1. Fork https://github.com/MultiQC/MultiQC and create a branch from `main`.
2. `git am 0001-rseqc-novel-combination-noncoding.patch`
3. Open a pull request. MultiQC's tests use example files from
   https://github.com/MultiQC/test-data; adding RSeQC 5.1.0
   `junction_annotation` and `read_distribution` output there helps review.
