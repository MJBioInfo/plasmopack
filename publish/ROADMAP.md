# plasmopack — Master Roadmap

**Last updated:** 2026-06-23
**Status:** Design complete for Phases 0-10. Ready to begin Phase 0 scaffolding.

This is the master build plan. Each phase has a clear deliverable ("Exit criterion").
Phases are built in order; later phases depend on earlier ones.

Detailed design documents live in `publish/decisions/`:
- `pathway_enrichment_strategy.md` (Phase 7)
- `tf_analysis_strategy.md` (Phase 8)
- `cci_strategy.md` (Phase 11)
- `deg_compositional_strategy.md` (Phases 5-6)

---

## Design principles (apply to every phase)

1. Minimum dependencies: required = anndata, numpy, pandas, scipy + stdlib only. All else optional extras.
2. HTTP via stdlib urllib (no requests/httpx).
3. AnnData is the only in-memory format. Results always written back with provenance.
4. Three-layer architecture: core (pure) → adapters (one per service) → public API.
5. Snapshot-first for fragile sources. Tests never hit network (VCR cassettes).
6. Confidence tiers on every result (experimental / inferred / predicted).
7. Provenance stamped in adata.uns["plasmopack"]["history"] on every write.

---

## Phase 0 — Foundation
**Build:** Empty but working package.
- Repo layout (src/plasmopack/, tests/, docs/, publish/)
- pixi.toml + pyproject.toml (hatchling)
- CI: ruff, mypy, pytest matrix (Mac/Linux/Windows × py3.11-3.13)
- Docs skeleton (mkdocs-material)
- Required deps locked: anndata, numpy, pandas, scipy
**Exit:** `pixi run test` passes on empty package.

## Phase 1 — Dummy data + AnnData conventions
**Build:** Test foundation.
- Synthetic P. berghei AnnData fixture (real PBANKA_ IDs, stage-correlated counts)
- AnnData slot convention as helper functions + provenance stamping
- Organism registry (organism → site → ID prefix → species_id)
**Exit:** dummy data loads, provenance helper works, registry returns correct config.

## Phase 2 — pp.io (I/O + GMT builder) [FLAGSHIP]
**Build:** Reading, writing, GMT feature.
- read_h5ad, write_h5ad, read_10x, read_gff, read_gaf
- build_gmt(): GAF → versioned GMT + provenance sidecar
- GeneSets in-memory object
- Seurat bridge (file-based, anndataR, documented)
**Exit:** build a GMT from PlasmoDB GAF and read it back.

## Phase 3 — pp.db (database clients)
**Build:** Durable adapter layer.
- Adapter framework (urllib + retry + disk cache + dataclass schemas + fallback snapshots)
- Stable REST clients: UniProt, Ensembl, NCBI, PlasmoDB/VEuPathDB, OrthoMCL, STRING
- VCR cassette tests each
**Exit:** pp.db.uniprot.fetch() works online, offline (cache), offline (snapshot).

## Phase 4 — pp.pp (preprocessing)
**Build:** QC and normalization for Plasmodium.
- QC metrics (apicoplast/mito %, gene/UMI counts) → adata.obs
- Filtering, normalization wrappers
- ID normalization (cross-namespace) → adata.var
**Exit:** dummy data runs full QC + normalization.

## Phase 5 — pp.tl: DEG
**Build:** Differential expression, multiple methods.
- Wilcoxon (scipy), Logistic regression (sklearn opt), Pseudo-bulk (pydeseq2 opt)
- deg_compare() multi-method concordance
- Results → adata.uns["deg"]
**Exit:** reproduces PSA Wilcoxon; pseudo-bulk runs on lane replicates.

## Phase 6 — pp.tl: Compositional analysis
**Build:** Cell-type proportion shifts.
- Simple proportions + chi-squared (scipy), CLR transform + tests (numpy)
- scCODA (optional, v0.2)
- Results → adata.uns["compositional"]
**Exit:** detects cluster proportion changes across PSA timepoints.

## Phase 7 — pp.tl: Enrichment (3 modes) [FLAGSHIP]
**Build:** Pathway pain-point solution.
- Mode A (PlasmoDB API), Mode B (local ORA/scipy), Mode C (local GSEA/gseapy opt)
- compare=True → confidence-tier output
- Results → adata.uns["enrichment"]
**Exit:** Mode A matches PSA plasmoDB_functional_analysis_pb.

## Phase 8 — pp.tl: TF analysis
**Build:** Transcription factor pipeline.
- tf_identify (ApicoTFdb snapshot + PlasmoDB GO + UniProt PF00847) → v0.1
- tf_build_regulons + tf_activity + tf_differential → v0.2
- Results → adata.var, adata.obsm["X_tf_activity"], adata.uns["regulons"]
**Exit:** marks AP2 factors in PSA; activity scoring where regulons exist.

## Phase 9 — pp.tl: Orthologs + unknown genes
**Build:** Annotation transfer.
- orthologs() (OrthoMCL), annotate_unknown() (ortholog-transfer, no BLAST)
- Results → adata.var
**Exit:** maps PBANKA_ unknown genes to human orthologs where they exist.

## Phase 10 — pp.pl (plotting)
**Build:** All visualizations.
- DEG: volcano, MA, heatmap, dotplot, violin, concordance, venn
- Compositional: stacked bar, box, bubble, coda, heatmap
- Enrichment: dotplot, agreement plot
- TF: activity heatmap, differential
- Combined: stage UMAP + composition, timepoint summary figure
**Exit:** every plot returns matplotlib Figure, smoke-tested.

## Phase 11 — pp.tl.cci (cell-cell interaction) [NOVEL]
**Build:** Parasite-host interaction (v0.2+).
- load_lr_database (IntAct + PHI-base + BioGRID) → v0.2
- score_interactions (CellPhoneDB method) → v0.3
- var_gene_activity + CCI plots → v0.4
**Exit:** LR database loads with evidence tiers; scoring validated against Poran et al.

## Phase 12 — Docs, tutorials, paper
**Build:** Publication package.
- Quickstart + 8 notebook tutorials
- Migration table (plasmoRUtils → plasmopack)
- JOSS paper draft
- API reference auto-generated
**Exit:** docs site live; paper draft complete.

## Phase 13 — Release + distribution
**Build:** Getting it to users.
- PyPI (trusted publishing), conda-forge, Bioconda
- Docker (multi-arch: Mac arm64, Linux amd64, Windows; GHCR)
- Zenodo DOI, CITATION.cff
**Exit:** pip install plasmopack, docker pull, JOSS submission.

## Phase 14 — Maintenance (ongoing)
- Monthly cron: refresh snapshots, run live-API tests, auto-PR on drift
- Quarterly dependency audit
- Yearly Python version bump
- Issue triage (Majeed + Claude)

---

## Version milestones

| Version | Includes |
|---|---|
| v0.1.0 | Phases 0-10 (core: io, db, pp, DEG, compositional, enrichment, TF-identify, orthologs, plotting) |
| v0.2.0 | TF full pipeline, scCODA, CCI LR database, combined figures |
| v0.3.0 | CCI scoring |
| v0.4.0 | CCI var gene + plots |
| v1.0.0 | Full feature set, JOSS paper published, stable API |

---

## Optional dependency groups (extras)

| Extra | Adds | Used by |
|---|---|---|
| plasmopack[deg] | pydeseq2 | Phase 5 pseudo-bulk |
| plasmopack[enrich] | gseapy | Phase 7 Mode C |
| plasmopack[comp] | sccoda | Phase 6 scCODA |
| plasmopack[sc] | scanpy | scanpy interop |
| plasmopack[ml] | scikit-learn | Phase 5 logistic regression |
| plasmopack[all] | everything above | full install |
