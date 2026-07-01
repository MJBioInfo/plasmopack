# DEG and Compositional Analysis Strategy — Design Decision

**Date:** 2026-06-23
**Status:** Decided
**Author:** Majeed Jamakhani + Claude
**Reference:** sc-best-practices.org (DEG and compositional chapters) — simplified for plasmopack

---

## DEG — Differential Gene Expression

### Core issue: pseudo-replication
The main pitfall in single-cell DEG is treating each cell as an independent replicate.
With 3 mice × 1000 cells, you have 3 replicates not 3000. Cell-level tests (Wilcoxon)
inflate significance. Pseudo-bulk (aggregate per replicate) is the recommended fix.

For PSA_minimal.h5ad: `lane` and `infection` are real replicates; `time` is condition.
Existing pb_deg_wilcoxon_time_* used cell-level Wilcoxon (exploratory, not ideal for final).

### Methods (simplest → most rigorous)

| Method | Test | Dep | When | Limitation |
|---|---|---|---|---|
| Wilcoxon | scipy.stats.ranksums | scipy (core) | exploration, no replicates | pseudo-replication |
| Logistic regression | sklearn | optional [ml] | control covariates | cell-level still |
| Pseudo-bulk | pydeseq2 NB-GLM | optional [deg] | publication, has replicates | needs replicates |
| NB-GLM (fallback) | scipy | scipy (core) | pydeseq2 unavailable | less established |

NOT included (need R): MAST, lme4 mixed models, muscat.

### DEG plots
volcano, MA, heatmap, dotplot, violin, ranked/waterfall, multi-method concordance scatter,
venn (method overlap). The concordance plot (Wilcoxon logFC vs pseudo-bulk logFC) is a
publishable QC figure — agreement = robust result.

### Storage
```
adata.uns["deg"][groupby][method] → DataFrame (gene, logFC, pval, padj, method, ...)
adata.uns["deg"][groupby]["comparison"] → DataFrame (gene, logFC_<m1>, logFC_<m2>, agreement)
```

---

## Compositional Analysis

### Core issue: compositionality
Cell-type proportions sum to 1 — if one goes up others must go down. Standard tests
assume independence, violated here. CLR transform or Bayesian models handle this.

Question for PSA: do cluster proportions shift between 24h / 48h / 62hpi?

### Methods

| Method | Test | Dep | When | Limitation |
|---|---|---|---|---|
| Simple proportions | scipy.stats.chi2_contingency | scipy (core) | quick look | ignores compositionality |
| CLR + test | numpy (log+center) + scipy | core | has replicates | needs replicates |
| scCODA | Bayesian Dirichlet-multinomial | optional [comp] | publication | needs ≥3 reps/condition, slower |

NOT included: ANCOM-BC (microbiome-focused), Dirichlet regression (redundant with scCODA).

### Compositional plots
stacked bar (always), box/violin of proportions, bubble (cell type × condition),
scCODA credible-interval plot, proportion heatmap.

### Storage
```
adata.uns["compositional"][condition_key][method] → DataFrame
   (cell_type, condition, proportion, effect, pval, padj)
```

---

## Combined DEG + Compositional figures (Plasmodium value)

| Figure | Shows | Value |
|---|---|---|
| Stage UMAP + composition pies | identity + abundance per condition | one figure, two questions |
| DEG heatmap + composition bar | what genes change + how many cells in that state | links expression to abundance |
| Timepoint summary | per timepoint: dominant clusters + top DEGs + enriched pathways | the paper's main narrative figure |

---

## Version plan

| Feature | v0.1 | v0.2 |
|---|---|---|
| DEG Wilcoxon, Logistic, Pseudo-bulk | ✓ | |
| DEG concordance/venn plots | ✓ | |
| Compositional simple + CLR | ✓ | |
| Compositional scCODA | | ✓ |
| Combined summary figures | | ✓ |

DEG and compositional both in v0.1 (only need scipy/numpy for core methods).
PSA Wilcoxon results validate the DEG module immediately.

---

## Paper talking points
1. Multi-method DEG with concordance QC — helps users see when pseudo-replication inflates results.
2. Compositional analysis aware of compositionality (CLR/scCODA) — often skipped in Plasmodium studies.
3. Combined figures tell the full story (identity + abundance + expression + pathway) in one place.
4. All simplified vs sc-best-practices — fewer methods, clearer defaults, AnnData-native, Plasmodium-ready.
