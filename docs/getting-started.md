# Getting started

## Installation

```bash
pip install plasmopack              # core (anndata, numpy, pandas, scipy)
pip install "plasmopack[all]"       # with all optional features
```

Optional feature groups:

| Extra | Adds | Enables |
|-------|------|---------|
| `plasmopack[deg]` | pydeseq2 | pseudo-bulk differential expression |
| `plasmopack[enrich]` | gseapy | GSEA enrichment mode |
| `plasmopack[comp]` | scCODA | Bayesian compositional analysis |
| `plasmopack[ml]` | scikit-learn | logistic-regression DEG |
| `plasmopack[sc]` | scanpy | scanpy interoperability |

## The namespace layout

```python
import plasmopack as pp

pp.io        # read/write: AnnData, GFF, GAF, GMT, Seurat bridge
pp.db        # database clients: VEuPathDB, UniProt, Ensembl, NCBI, ...
pp.pp        # preprocessing: QC, normalization, ID normalization
pp.tl        # tools: DEG, compositional, enrichment, TF, orthologs
pp.pl        # plotting
pp.datasets  # example data + organism registry
```

## Supported organisms

```python
from plasmopack.datasets import list_organisms, get_organism

list_organisms()
# ['pberghei', 'pfalciparum', 'pknowlesi', 'pvivax', 'tgondii']

get_organism("pb")     # aliases work too
```
