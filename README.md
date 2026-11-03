# plasmopack

[![CI](https://github.com/MJBioInfo/plasmopack/actions/workflows/ci.yml/badge.svg)](https://github.com/MJBioInfo/plasmopack/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

**AnnData-native toolkit for *Plasmodium* and VEuPathDB organism bioinformatics.**

`plasmopack` brings the scanpy-style workflow (`pp`, `tl`, `pl`) to malaria and other
non-model eukaryotic pathogens. Query databases, build gene sets, run enrichment, analyse
transcription factors — all without ever leaving your `AnnData` object.

> **Status:** early development (v0.1.0.dev). API may change before v1.0.

```python
import plasmopack as pp

adata = pp.io.read_h5ad("my_parasite_data.h5ad")

pp.pp.qc_metrics(adata)                       # QC (apicoplast %, gene counts)
pp.tl.deg(adata, groupby="time")              # differential expression
pp.tl.enrich(adata, gene_sets="GO")           # pathway enrichment (3 modes)
pp.pl.deg_volcano(adata)                      # plot — results stay in AnnData

adata.write("annotated.h5ad")                 # everything saved, self-documenting
```

## Design principles

- **Minimum dependencies.** Core install needs only `anndata`, `numpy`, `pandas`, `scipy`.
  Everything else is an optional extra — so the package keeps working when the ecosystem changes.
- **AnnData is the only data model.** Every result is written back into standard AnnData slots
  with full provenance (which database version, which parameters, when).
- **Built to last.** Fragile data sources ship as versioned offline snapshots; tests never
  depend on live services.
- **VEuPathDB-wide.** *Plasmodium* first, but organism-parameterised from day one
  (Toxoplasma, Cryptosporidium, Trypanosoma, and more).

## Installation

**Requires Python 3.11, 3.12, or 3.13.** (Not yet 3.14+ — the scientific stack
has no binary wheels there yet.)

From PyPI (once released):

```bash
pip install plasmopack                # core
pip install "plasmopack[all]"         # with all optional features
```

From GitHub (current development version):

```bash
pip install "git+https://github.com/MJBioInfo/plasmopack.git"
```

### Recommended: a clean environment on a supported Python

Using a fresh virtual environment avoids conflicts and the HDF5 build error below:

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install "git+https://github.com/MJBioInfo/plasmopack.git"
```

Or with conda / mamba (gets a pre-built HDF5, no compiler needed):

```bash
conda create -n myapp python=3.12 anndata -c conda-forge
conda activate myapp
pip install "git+https://github.com/MJBioInfo/plasmopack.git"
```

### Troubleshooting: `Unable to load dependency HDF5` / `Failed building wheel for h5py`

This means you are on an unsupported Python (usually 3.14/3.15) where `h5py` has
no ready-made wheel and tries to compile from source. **Fix:** use Python 3.11–3.13
(see the clean-environment commands above). You do **not** need to install HDF5 yourself.

## Documentation

See [the documentation site](https://MJBioInfo.github.io/plasmopack) (coming soon).

## Citation

If you use `plasmopack`, please cite it — see [`CITATION.cff`](CITATION.cff).

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgement

Inspired by [plasmoRUtils](https://github.com/Rohit-Satyam/plasmoRUtils) (R) by Rohit Satyam.
