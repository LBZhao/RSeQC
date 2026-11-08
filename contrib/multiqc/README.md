# MultiQC support for RSeQC 5.1.0+lbzhao output

RSeQC 5.1.0+lbzhao adds two report categories that MultiQC (tested with 1.35)
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

- `multiqc --strict` on RSeQC 5.1.0+lbzhao output: junction categories sum to
  100% and `Noncoding_Exons` appears in the plot.
- `multiqc --strict` on RSeQC 5.0.6 output: report unchanged.
- `ruff check` / `ruff format`: no new findings.

## Using or submitting the patch

The official RSeQC does not write these lines; only this fork does. So:

- **To use it now**, apply the patch to your own MultiQC install:
  `git clone https://github.com/MultiQC/MultiQC && cd MultiQC &&
  git am <path>/0001-rseqc-novel-combination-noncoding.patch && pip install .`
- **To offer it to MultiQC**, fork https://github.com/MultiQC/MultiQC,
  `git am` the patch on a branch from `main`, and open a pull request that
  says the output comes from the LBZhao fork of RSeQC. MultiQC may decline
  support for a fork. Example output for https://github.com/MultiQC/test-data
  helps review.
