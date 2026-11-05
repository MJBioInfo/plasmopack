# Usage guide

A complete walk-through of the analysis workflow. Everything operates on an
`AnnData` object and writes results back into it with provenance.

```python
import plasmopack as pp

adata = pp.io.read_h5ad("my_parasite_data.h5ad")   # PBANKA_/PF3D7_ gene IDs
```

---

## 1. Preprocessing — `pp.pp`

```python
pp.pp.qc_metrics(adata, gene_groups={"mito": mito_genes})  # -> adata.obs
pp.pp.filter_cells(adata, min_genes=200)
pp.pp.filter_genes(adata, min_cells=3)
pp.pp.normalize(adata, target_sum=1e4)             # -> adata.layers["lognorm"]
```

`normalize` is **non-destructive**: it writes the log-normalized matrix to a
layer and leaves raw counts in `adata.X` (pseudobulk DEG needs them).

---

## 2. Differential expression — `pp.tl.deg`

**The user chooses the method.** See what's available and what each expects:

```python
pp.tl.list_deg_methods()
```

| method | level | needs | expects (which data) |
|--------|-------|-------|----------------------|
| `wilcoxon` (default) | cell | core | log-normalized (`layer="lognorm"`) |
| `t-test` | cell | core | log-normalized |
| `pseudobulk` | replicate-aware | core | **raw counts** + `sample_key` |
| `deseq2` | replicate-aware | `[deg]` | **raw counts** + `sample_key` |
| `logreg` | cell | `[ml]` | log-normalized |
| `mast` | cell | R (todo) | log-normalized |

```python
# exploratory, cell-level (uses the lognorm layer automatically)
pp.tl.deg(adata, "time", method="wilcoxon")

# publication-grade, replicate-aware (needs a biological-replicate column)
pp.tl.deg(adata, "time", method="deseq2", sample_key="lane")
```

**Which matrix/layer?** You don't have to think about it — `deg` picks the right
default layer for the method (`lognorm` for cell-level, `counts` for pseudobulk)
and **warns** if you hand it the wrong kind of data. Override with `layer=...`.

Results: `adata.uns["deg"][groupby][method]` — a tidy table with `group, gene,
logFC, pval, padj, pct_grp, pct_rest`.

---

## 3. Pathway enrichment — `pp.tl.enrich`

### Zero-setup: just name the collection

```python
pp.tl.enrich(adata, gene_sets="GO", groupby="time", group="24h", direction="up")
```

`"GO"` loads the **bundled default gene sets** for the organism, which is
**auto-detected from the gene IDs** (`PBANKA_…` → *P. berghei*). Collections:

| `gene_sets=` | meaning |
|--------------|---------|
| `"GO"` | all Gene Ontology (BP + CC + MF) |
| `"BP"` / `"CC"` / `"MF"` | one GO aspect |
| a `GeneSets` object | your own sets |
| a path to a `.gmt` | your own file |

Bundled organisms: **P. berghei, P. falciparum, P. knowlesi, T. gondii**
(`pp.datasets.list_bundled_gene_sets()`). For others, or for the latest
annotations, build them live (below) or pass your own GMT.

### How genes are selected (your choice)

```python
# default: combined cutoff — significant AND strong
pp.tl.enrich(adata, "GO", groupby="time", group="24h",
             padj_cutoff=0.05, lfc_cutoff=1.0, direction="up")

# OR top-N ranked — pick the ranking
pp.tl.enrich(adata, "GO", groupby="time", group="24h",
             top_n=50, rank_by="combined")   # "padj" | "logfc" | "combined"
```

`"combined"` = `-log10(padj) × |logFC|` (significance weighted by effect size —
the recommended single ranking).

Results: `adata.uns["enrichment"][label]` — `set_id, name, set_size, n_overlap,
odds_ratio, pval, padj, overlap_genes`, sorted by `padj`. The test is Fisher's
exact (same as PlasmoDB's web tool), BH-adjusted.

### Getting / refreshing the default gene sets

```python
gs = pp.datasets.gene_sets("pberghei", "GO")          # load bundled snapshot
gs = pp.datasets.build_gene_sets("pberghei",          # fetch LATEST from UniProt
                                 save_to="pb_GO.gmt")
```

Every snapshot is version-stamped; `enrich` records which version was used.

---

## 4. Plotting — `pp.pl`

### Circos pathway ↔ gene chord

```python
fig = pp.pl.pathway_chord(
    adata, "time:24h:up",            # the enrichment result key
    groupby="time", group="24h",     # where gene fold-changes come from
    n_pathways=3, n_genes=20,         # how many to show (configurable)
    gene_rank="combined",            # which genes to keep
)
fig.savefig("chord.png", dpi=150, bbox_inches="tight")
```

Pathways are arcs on one side, their genes on the other; genes are coloured
**red (up)** / **blue (down)** by fold change.

---

## 5. Database lookups — `pp.db`

```python
pp.db.uniprot.fetch("C6KSZ8")                       # by accession
pp.db.uniprot.search("AMA1", organism="pfalciparum") # by gene name
pp.db.ensembl.lookup("PF3D7_0206800")                # gene location/symbol
```

All cached on disk (offline + reproducible after the first call).

---

## Reproducibility

Every function stamps `adata.uns["plasmopack"]["history"]` — the function,
parameters, data versions (e.g. the gene-set snapshot), timestamp, and package
version. A saved `.h5ad` is a complete record of how its results were produced.
