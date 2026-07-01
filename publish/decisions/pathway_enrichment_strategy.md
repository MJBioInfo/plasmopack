# Pathway Enrichment Strategy — Design Decision

**Date:** 2026-06-23  
**Status:** Decided  
**Author:** Majeed Jamakhani + Claude  

---

## The Problem We Solve

Standard Plasmodium pathway analysis requires:
1. Manually extracting DEG gene lists from AnnData
2. Uploading gene lists to PlasmoDB website by hand
3. Downloading results as tables
4. Manually pasting results back into AnnData

This workflow is:
- Not reproducible (no record of which PlasmoDB release was used)
- Not automated (manual steps cannot be scripted)
- Limited to one method (ORA only — misses subtler pathway effects)
- Not comparable across methods (no way to assess result confidence)

plasmopack eliminates all of these problems.

---

## Three Core Concepts (Must Understand Before Reading Further)

### 1. Gene-Pathway Mapping Source

Before computing enrichment, you need a table of "which genes belong to which pathway".
For Plasmodium this comes from PlasmoDB, which publishes:
- **GOA file** — Gene Ontology Annotations: maps each gene to GO terms (BP, CC, MF)
- **Pathway mapping file** — maps each gene to MetaCyc/KEGG pathways

PlasmoDB publishes these files per release (release 66, 67, 68...).

### 2. GMT File (Gene Matrix Transposed)

A GMT file is the pathway-centric view of the GOA file:
- Each row = one pathway + all its member genes
- Standard format accepted by all enrichment tools (GSEA, gseapy, decoupler)
- Built once from the GOA file, cached with version stamp

```
GO:0006412  translation  PBANKA_0314500  PBANKA_0405400  PBANKA_0807600 ...
GO:0003735  ribosome     PBANKA_0314500  PBANKA_0619100  PBANKA_0941500 ...
```

**plasmopack's GMT builder converts PlasmoDB GOA files → versioned GMT files.**
This is a standalone feature useful even outside the enrichment workflow.

### 3. Statistical Method

Two fundamentally different approaches:

**ORA (Over-Representation Analysis):**
- Input: list of significant DEGs (after logFC/padj cutoff)
- Test: Fisher's exact test per pathway
- Question: "Is pathway X over-represented in my DEG list?"
- Used by: PlasmoDB web tool, topGO (R), plasmoRUtils

**GSEA (Gene Set Enrichment Analysis / fGSEA):**
- Input: ALL genes ranked by logFC (no cutoff)
- Test: permutation-based enrichment score
- Question: "Do pathway X genes tend to be at top or bottom of my ranked list?"
- More sensitive — detects subtle pathway shifts where individual genes don't pass significance cutoff
- Used by: GSEA software, fGSEA (R), gseapy (Python)

---

## The Three Modes in plasmopack

### Mode A — PlasmoDB API (Remote ORA)
- Sends DEG list to PlasmoDB's REST API
- PlasmoDB runs Fisher's exact test on their servers with their internal annotations
- Returns results in same format as manual web submission
- **Advantage:** Identical to what community is used to; validates against existing results
- **Disadvantage:** Requires internet; database version not explicitly recorded; opaque background

### Mode B — Local ORA (GMT-based)
- Uses GMT file built from PlasmoDB GOA (cached, versioned)
- Runs Fisher's exact test locally with scipy.stats.fisher_exact (stdlib + scipy, no new dep)
- Background = all genes in the GMT file (same as PlasmoDB's background)
- **Advantage:** Fully offline after first download; version-stamped; reproducible; fast
- **Disadvantage:** Small differences from Mode A possible (rounding, exact background definition)

### Mode C — Local GSEA (GMT-based)
- Uses same GMT file as Mode B
- Runs fGSEA on all genes ranked by logFC (requires gseapy as optional dep)
- No DEG cutoff needed — uses all genes
- **Advantage:** More sensitive; catches pathways missed by ORA; complementary to Modes A and B
- **Disadvantage:** Needs ranked list of all genes; requires optional gseapy install

---

## Multi-Mode Comparison — The Key Scientific Advantage

Running all three modes and comparing their output gives a **confidence tier** for each pathway:

| Found in | Confidence | Interpretation |
|---|---|---|
| Mode A + B + C | Highest | Robust finding, consistent across methods |
| Mode A + B only | High | ORA-specific; strong effect in DEG list |
| Mode C only | Medium | Subtle shift; real but below DEG cutoff |
| Mode A only | Low | Check: may be background/version difference |
| Mode B only | Low | Check: may be background/version difference |

This confidence tiering is novel — no existing tool for Plasmodium provides it. It is a genuine scientific contribution and a clear talking point for the paper.

---

## Validation Against Existing Data

The PSA_minimal.h5ad dataset (P. berghei, liver stage, 3081 cells × 5254 genes) contains:
- `adata.uns['plasmoDB_functional_analysis_pb']` — results from manual PlasmoDB submission
- Structure: timepoint (24/48/62hpi) → category (BP/CC/MF/KEGG) → DataFrame (12 columns)
- This is our ground truth for Mode A validation

**Validation test:**
Run `pp.tl.enrich(adata, mode="plasmodb_api")` on 24h DEGs.
Compare results to `plasmoDB_functional_analysis_pb["24hpi"]["BP"]`.
Results should match (same pathways, same direction, similar p-values).

This validation is reported in the paper as evidence that plasmopack reproduces manual PlasmoDB analysis, and adds Mode B and Mode C on top.

**Known gap in PSA data:** The PlasmoDB release version used for manual submission was not recorded in the file. Our tool records this automatically in `adata.uns["plasmopack"]["history"]`.

---

## How plasmoRUtils (R) Does It — Comparison

| Feature | plasmoRUtils | plasmopack |
|---|---|---|
| GO enrichment method | topGO (local, elim/weight algorithm) | Mode B (Fisher's exact) + Mode C (fGSEA) |
| MetaCyc/KEGG | Scrapes MPMP website (fragile) | PlasmoDB pathway mapping file (stable) |
| PlasmoDB API mode | Not available | Mode A |
| Multi-mode comparison | Not available | Yes |
| Results in AnnData | Not available (R object) | Yes — adata.uns["enrichment"] |
| Version recording | Not available | Yes — provenance in adata.uns |
| GMT file support | Not available | Yes — flagship feature |
| Offline after first run | Yes (topGO local) | Yes (Mode B and C) |

**Key difference from topGO:** topGO's "elim" algorithm accounts for GO hierarchy (avoids double-counting parent terms when child is already significant). Our Mode B uses simple Fisher's exact like PlasmoDB does. For v0.2 we can implement GO hierarchy correction. This is documented as a known difference.

---

## AnnData Storage Convention

```
adata.uns["enrichment"]
    ├── "24hpi"
    │     ├── "plasmodb_api"      # Mode A results
    │     │     ├── "GO_BP"       → DataFrame (ID, Name, pval, padj, fold_enrichment, genes)
    │     │     ├── "GO_CC"       → DataFrame
    │     │     ├── "GO_MF"       → DataFrame
    │     │     └── "MetaCyc"     → DataFrame
    │     ├── "local_ora"         # Mode B results (same structure)
    │     ├── "local_gsea"        # Mode C results (NES, pval, padj, leading_edge)
    │     └── "comparison"        → DataFrame (pathway, found_in_modes, confidence_tier)
    ├── "48hpi"
    │     └── ... (same structure)
    └── "62hpi"
          └── ... (same structure)

adata.uns["plasmopack"]["history"]
    └── [{"function": "tl.enrich", "date": "2026-06-23", 
           "gmt_version": "PlasmoDB-68", "modes": ["local_ora","local_gsea"],
           "organism": "pberghei", "package_version": "0.1.0"}]
```

---

## Paper Talking Points

1. **Reproducibility gap in current practice:** Manual PlasmoDB submission loses database version, is not scriptable, and cannot be integrated into analysis pipelines. plasmopack solves this.

2. **Method comparison adds confidence:** No existing Plasmodium tool runs ORA and GSEA and compares them. The confidence-tier output helps users identify the most robust findings.

3. **GMT as a reusable artifact:** The GMT file built by plasmopack can be shared, version-controlled, and used with any enrichment tool outside plasmopack. It is a community resource, not just an internal file.

4. **Backward compatibility with PlasmoDB results:** Mode A reproduces existing manually obtained results exactly, allowing plasmopack to be used alongside existing datasets like PSA_minimal.h5ad.

5. **Extension to all VEuPathDB organisms:** The same GMT-building and enrichment pipeline works for any organism in VEuPathDB by changing `organism=` parameter.

---

## Open Questions (Not Yet Decided)

- GO hierarchy correction (topGO elim equivalent): v0.1 does not include it; document clearly; add in v0.2
- Custom background gene list: allow users to specify background (e.g. only expressed genes, not all annotated genes) — important for single-cell where not all genes are detected
- Minimum gene set size filter: default min=5, max=500 (standard GSEA practice)
