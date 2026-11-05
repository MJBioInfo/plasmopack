# API reference

The public API is organised into scanpy-style submodules, all operating on
`AnnData`.

## `pp.io` — input/output & gene sets

::: plasmopack.io.read_h5ad
::: plasmopack.io.write_h5ad
::: plasmopack.io.read_gaf
::: plasmopack.io.read_gmt
::: plasmopack.io.write_gmt
::: plasmopack.io.GeneSets

## `pp.pp` — preprocessing

::: plasmopack.pp.qc_metrics
::: plasmopack.pp.filter_cells
::: plasmopack.pp.filter_genes
::: plasmopack.pp.normalize

## `pp.tl` — tools

::: plasmopack.tl.deg
::: plasmopack.tl.list_deg_methods
::: plasmopack.tl.enrich

## `pp.pl` — plotting

::: plasmopack.pl.pathway_chord

## `pp.db` — database clients

::: plasmopack.db.uniprot.fetch
::: plasmopack.db.uniprot.search
::: plasmopack.db.ensembl.lookup

## `pp.datasets` — defaults & registry

::: plasmopack.datasets.gene_sets
::: plasmopack.datasets.build_gene_sets
::: plasmopack.datasets.detect_organism
::: plasmopack.datasets.list_bundled_gene_sets
::: plasmopack.datasets.get_organism
::: plasmopack.datasets.list_organisms
