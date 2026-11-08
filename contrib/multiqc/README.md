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

## The LBZhao/MultiQC fork

The two patches here make up the fork at https://github.com/LBZhao/MultiQC
(official MultiQC plus these commits):

1. `0001-rseqc-novel-combination-noncoding.patch`: the RSeQC module change.
2. `0002-identify-lbzhao-fork.patch`: marks the fork as not official
   (version `1.36.dev0+lbzhao`, shown in every report; README notice).

To recreate it, fork https://github.com/MultiQC/MultiQC on GitHub, then:

```bash
git clone https://github.com/LBZhao/MultiQC && cd MultiQC
git am <path>/0001-rseqc-novel-combination-noncoding.patch <path>/0002-identify-lbzhao-fork.patch
git push origin main
```

This applies the patches on top of the fork's current `main` (they were
made against `a286ece`; if `git am` reports a conflict, the official
MultiQC changed the same lines and the patch needs updating).

Install it with `pip install git+https://github.com/LBZhao/MultiQC.git`.

To offer only the RSeQC change to the official MultiQC, open a pull request
with patch 0001 alone and say the output comes from the LBZhao fork of
RSeQC; MultiQC may decline support for a fork.
