# plasmopack

**AnnData-native toolkit for *Plasmodium* and VEuPathDB organism bioinformatics.**

`plasmopack` brings the scanpy-style workflow (`pp`, `tl`, `pl`) to malaria and other
non-model eukaryotic pathogens. Query databases, build gene sets, run enrichment, and
analyse transcription factors — all without leaving your `AnnData` object.

!!! warning "Early development"
    This package is at `v0.1.0.dev`. The API may change before `v1.0`.

## Why plasmopack?

- **Minimum dependencies.** Core install needs only `anndata`, `numpy`, `pandas`, `scipy`.
- **AnnData-native.** Every result is written back into standard slots with full provenance.
- **Built to last.** Fragile sources ship as versioned offline snapshots; tests never
  depend on live services.
- **VEuPathDB-wide.** *Plasmodium* first, organism-parameterised from day one.

## Quick example

```python
import plasmopack as pp

adata = pp.io.read_h5ad("my_parasite_data.h5ad")
pp.pp.qc_metrics(adata)
pp.tl.enrich(adata, gene_sets="GO")
pp.pl.deg_volcano(adata)
```

See [Getting started](getting-started.md).
