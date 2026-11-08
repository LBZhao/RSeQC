# RSeQC (LBZhao fork)

[![License](https://img.shields.io/github/license/LBZhao/RSeQC)](LICENSE)
[![CI](https://github.com/LBZhao/RSeQC/actions/workflows/ci.yml/badge.svg)](https://github.com/LBZhao/RSeQC/actions/workflows/ci.yml)

This is an independent fork of [RSeQC](https://github.com/liguowang/RSeQC)
by Liguo Wang, branched from RSeQC 5.0.6. It is **not** the official RSeQC
and is not published on PyPI: `pip install rseqc` installs the original.

The fork's version is `5.1.0+lbzhao`. It fixes several bugs that change
results (e.g. soft-clipped reads were shifted, supplementary alignments were
counted as reads, exon-skipping junctions were called known), rewrites gene
body coverage, and no longer needs R. See the 5.1.0 release notes in
[`doc/index.rst`](doc/index.rst) for every change. Results differ from RSeQC
5.0.6; when you publish results, say which one you used.

## Requirements

- Python 3.9 or later
- `pip`

Check your Python and pip installations:

```bash
python3 --version
python3 -m pip --version
```

Example:

```text
Python 3.11.13
pip 25.3 from /path/to/lib/python3.11/site-packages/pip (python 3.11)
```

If `pip` is not installed, follow the official installation instructions:

https://pip.pypa.io/en/stable/installation/

---

## (Optional) Create a Virtual Environment

Using a virtual environment is recommended to keep RSeQC and its dependencies isolated from other Python packages.

Create a new environment:

```bash
python3 -m venv rseqc_env
```

Activate it:

**Linux/macOS**

```bash
source rseqc_env/bin/activate
```

**Windows (PowerShell)**

```powershell
rseqc_env\Scripts\Activate.ps1
```

---

## Install this fork

Install from GitHub (this fork is not on PyPI):

```bash
python3 -m pip install git+https://github.com/LBZhao/RSeQC.git
```

You can also install a specific branch, tag, or commit, for example a
commit hash for a reproducible analysis:

```bash
python3 -m pip install git+https://github.com/LBZhao/RSeQC.git@<commit>
```

The fork and the original use the same package name, so installing one
replaces the other in the same environment. To use the original RSeQC
instead, run `python3 -m pip install rseqc`.

---

## Upgrade RSeQC

```bash
python3 -m pip install --upgrade git+https://github.com/LBZhao/RSeQC.git
```

---

## Uninstall RSeQC

Remove RSeQC:

```bash
python3 -m pip uninstall rseqc
```

---

## Verify the Installation

Check that RSeQC is installed correctly:

```bash
bam_stat.py --help
```

or

```bash
python3 -c "from importlib.metadata import version; print('RSeQC', version('RSeQC'))"
```

---

## License

RSeQC is licensed under the **GNU General Public License v3.0 or later (GPL-3.0-or-later)**.

See the `LICENSE` file included with this project for the complete license text.

---

## Documentation

This fork's documentation, including what changed from 5.0.6, is in
[`doc/index.rst`](doc/index.rst).

Documentation for the original RSeQC:

1. http://rseqc.sourceforge.net/
2. https://www.lwang.org/RSeQC/

## References

If you use this fork, please cite the original RSeQC publications below and
mention that you used the LBZhao fork (with the version or commit).

1. Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments. Bioinformatics. 2012;28(16):2184-2185. doi:10.1093/bioinformatics/bts356
2. Wang L, Nie J, Sicotte H, et al. Measure transcript integrity using RNA-seq data. BMC Bioinformatics. 2016;17:58. doi:10.1186/s12859-016-0922-z

