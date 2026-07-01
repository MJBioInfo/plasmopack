# Transcription Factor Analysis Strategy — Design Decision

**Date:** 2026-06-23
**Status:** Decided
**Author:** Majeed Jamakhani + Claude

---

## The Problem We Solve

Identifying and analysing transcription factors (TFs) in Plasmodium is uniquely difficult:
- Standard TF databases (JASPAR, TRANSFAC, ENCODE) are built for human/mouse — they have no Plasmodium entries
- Generic tools (SCENIC, decoupler) have no Plasmodium regulon data
- plasmoRUtils offers only a lookup (which genes are TFs) — no activity scoring, no AnnData integration
- The Plasmodium TF family is biologically unusual — dominated by AP2 domain proteins not found in mammals

plasmopack provides a complete pipeline: identify TFs → annotate in AnnData → score activity per cell → find differentially active TFs across stages/conditions.

---

## The Biology (Essential Background for the Paper)

### The AP2 / ApiAP2 TF Family

In Plasmodium, the dominant transcription factor family is the **ApiAP2 (Apicomplexan AP2) proteins**.
- Named after the Apetala2 (AP2) DNA-binding domain (Pfam PF00847)
- ~27 ApiAP2 factors in *P. falciparum* (Pf3D7)
- ~24 ApiAP2 factors in *P. berghei* (PBANKA_)
- They control stage-specific gene expression across the entire life cycle
- Other TF families (e.g. BDP1, a bromodomain protein) exist but are minor
- Standard mammalian TF families (zinc fingers like SP1, nuclear receptors, bHLH) are absent or rare

**Why this matters for the package:** any tool that uses a human/mouse TF database will find essentially zero Plasmodium TFs. We must use Plasmodium-specific resources.

### Regulon = TF → Target Gene Links

A regulon is the set of genes that a specific TF controls. To score TF activity in a cell, you need both:
1. Which genes are TFs (identification)
2. Which genes each TF regulates (regulon)

For human data, rich regulon databases exist (DoRothEA, CollecTRI cover hundreds of TFs with thousands of targets). For Plasmodium, regulon data is sparse:
- ~8–10 of 27 P. falciparum AP2 factors have published ChIP-seq data
- For P. berghei (our PSA dataset organism), even fewer ChIP-seq datasets exist
- For other VEuPathDB organisms (Toxoplasma, Leishmania), some AP2 ortholog data exists

This is the honest limitation. We document it and work around it (see Three-Step Pipeline below).

---

## Databases and Data Sources

### Primary Source: ApicoTFdb

- URL: bioinfo.icgeb.res.in/ApicoTFdb/
- What it contains: curated list of all apicomplexan TFs, classified by domain family, organism, evidence
- Coverage: P. falciparum, P. berghei, T. gondii, C. parvum, B. bovis and others
- Quality: manually curated, peer-reviewed (published paper), gold standard for the community
- Access: web interface only (no REST API) → we ship a versioned parquet snapshot
- Update frequency: infrequent (major updates with new publications)

### Secondary Sources

| Source | What it adds | Access | Quality |
|---|---|---|---|
| PlasmoDB GO:0003700 | Genes annotated as "DNA-binding TF activity" | REST API | Derives from ApicoTFdb + literature |
| UniProt Pfam PF00847 | All proteins with AP2 domain (catches any missed by ApicoTFdb) | REST API | Stable, EBI-maintained |
| Published ChIP-seq (GEO) | TF → target gene links for ~8-10 AP2 factors | Static files, PlasmoDB downloads | Highest quality; incomplete |
| PlasmoDB gene expression | Co-expression for inferring missing regulons | REST API | Lower confidence; clearly flagged |

### Published ChIP-seq Datasets Available (for regulon building)

These are real published datasets that provide TF-target links. This list should be updated
as new papers are published:

- AP2-G (gametocyte commitment) — Sinha et al., Nature 2014; López-Rubio lab
- AP2-I (invasion genes) — Santos et al., Cell Host Microbe 2017
- AP2-EXP (asexual stages) — Campagnac & Bhatt labs
- PfBDP1 — Josling et al., Cell Host Microbe 2015
- AP2-G2 — Yuda lab
- AP2-MG — Campagnac lab
- AP2-O — Yuda lab (ookinete stage)
- AP2-SP — sporogony stage

For P. berghei specifically: some orthologous ChIP-seq data exists but is sparser. Ortholog
mapping from Pf to Pb is used as fallback with explicit confidence flag.

---

## Three-Step Pipeline: pp.tl.tf_*

### Step 1 — Identify TFs: pp.tl.tf_identify(adata)

```
Sources queried (in order):
  1. ApicoTFdb snapshot (primary)
  2. PlasmoDB GO:0003700 annotation (confirms and extends)
  3. UniProt PF00847 domain (catches any missed)
  4. Merge → deduplicate → confidence score per gene

Results written to:
  adata.var["is_tf"]          True/False
  adata.var["tf_family"]      "AP2", "BDP", "other", NaN
  adata.var["tf_name"]        "AP2-G", "AP2-I", etc. (where named)
  adata.var["tf_confidence"]  "curated" / "domain_only" / "go_only"
  adata.uns["plasmopack"]["history"] → records ApicoTFdb snapshot version
```

This step is equivalent to what plasmoRUtils does with searchApicoTFdb, but:
- Uses a snapshot (does not scrape live — more robust)
- Writes results into AnnData
- Adds confidence levels
- Cross-validates across three sources

### Step 2 — Build Regulons: pp.tl.tf_build_regulons(adata)

```
For each TF in adata.var["is_tf"]:

  IF ChIP-seq data available (known ~8-10 Pf AP2s):
    → load target gene list from bundled dataset (GEO-derived, versioned)
    → regulon_confidence = "experimental"

  ELIF ortholog in Pf with ChIP-seq (for Pb AP2s with Pf ChIP-seq ortholog):
    → map targets via OrthoMCL Pf→Pb
    → regulon_confidence = "ortholog_inferred"

  ELSE (most AP2 factors):
    → compute co-expression network (Spearman, top N correlated genes)
    → regulon_confidence = "coexpression_inferred"
    → clearly flagged as hypothesis, not experimental

Results written to:
  adata.uns["regulons"][tf_name] = {
      "targets": [gene_id, ...],
      "n_targets": int,
      "confidence": "experimental" / "ortholog_inferred" / "coexpression_inferred",
      "source": "ChIP-seq paper DOI" / "OrthoMCL Pf→Pb" / "coexpression"
  }
```

**Why co-expression as fallback is scientifically defensible:**
Co-expression-based regulon inference is standard practice in single-cell analysis
(SCENIC, pySCENIC, decoupler all use this). We are transparent about confidence level.
Reviewers accept this if the confidence flag is shown clearly in results and plots.

### Step 3 — Score TF Activity Per Cell: pp.tl.tf_activity(adata)

```
Method: mean normalised expression of regulon target genes per cell
(same principle as VIPER/decoupler's mean method — simple, interpretable, reproducible)

For each cell and each TF with a regulon:
  activity_score = mean(normalised_expression[target_genes])

Results written to:
  adata.obsm["X_tf_activity"]    shape: (n_cells, n_tfs_with_regulon)
  Columns = TF names; rows = cells

Differential TF activity across timepoints / clusters:
  pp.tl.tf_differential(adata, groupby="time")
  → adata.uns["tf_activity"]["differential"]["time"]
     DataFrame: tf_name, mean_24h, mean_48h, mean_62h, padj, confidence
```

---

## How We Compare to plasmoRUtils

| Capability | plasmoRUtils | plasmopack |
|---|---|---|
| Identify TFs | ✓ (live scrape) | ✓ (snapshot, more robust) |
| Cross-validate across sources | ✗ | ✓ (3 sources) |
| Write TF labels to AnnData | ✗ | ✓ adata.var["is_tf"] |
| Build regulons | ✗ | ✓ (experimental + inferred) |
| Score TF activity per cell | ✗ | ✓ adata.obsm["X_tf_activity"] |
| Differential TF activity | ✗ | ✓ |
| Confidence levels | ✗ | ✓ at every step |
| Works offline | ✗ | ✓ (after snapshot download) |

---

## Honest Limitations for Reviewers

1. **Regulon coverage is incomplete.** Only ~8-10 AP2 factors have experimental
   ChIP-seq regulons. The rest use co-expression inference, which is lower confidence.
   We do not hide this — every result carries a confidence label.

2. **P. berghei regulon data is sparser than P. falciparum.** The PSA dataset uses
   Pb. Regulons are either inferred from Pf orthologs or from co-expression.

3. **We do not do motif scanning.** SCENIC/pySCENIC does motif-based regulon
   inference (finds TF binding sites in promoters). This requires genome sequence +
   motif database + BLAST — too many dependencies. We use expression-based activity
   scoring, which is well-established and simpler.

4. **We do not claim to discover new TFs.** All identifications come from ApicoTFdb
   and GO annotations. Novel TF discovery is out of scope.

---

## AnnData Storage Convention (Complete)

```
adata.var
  ├── "is_tf"              bool
  ├── "tf_family"          str or NaN
  ├── "tf_name"            str or NaN
  └── "tf_confidence"      "curated" / "domain_only" / "go_only" / NaN

adata.obsm
  └── "X_tf_activity"      ndarray (n_cells × n_tfs_with_regulon)

adata.uns
  ├── "regulons"
  │     └── "AP2-G" → {targets, n_targets, confidence, source}
  │     └── "AP2-I" → {targets, n_targets, confidence, source}
  │     └── ...
  ├── "tf_activity"
  │     └── "differential"
  │           └── "time" → DataFrame (tf, mean_per_group, padj, confidence)
  └── "plasmopack"
        └── "history" → records ApicoTFdb version, ChIP-seq dataset versions, date
```

---

## Paper Talking Points

1. **First AnnData-native TF pipeline for Plasmodium.** No existing Python tool
   integrates ApicoTFdb, ChIP-seq regulons, and single-cell activity scoring for
   Plasmodium in one workflow.

2. **Transparent confidence tiers.** Every TF identification and every activity score
   carries an explicit confidence label. Reviewers and users know exactly what is
   experimental and what is inferred.

3. **Comparison to plasmoRUtils:** plasmoRUtils provides only TF identification via
   web scraping. plasmopack adds regulon-based activity scoring and differential analysis
   — a qualitative improvement in analytical depth.

4. **Validation strategy:** Using the PSA dataset (P. berghei liver stage, 3 timepoints),
   we identify which AP2 factors are differentially active between 24h, 48h, and 62hpi.
   Results are cross-referenced against published stage-specific AP2 expression patterns
   (Yuda, Llinas lab papers).

---

## Version / Release Plan

- **v0.1.0:** tf_identify() only (Step 1). Safe, solid, validated.
- **v0.2.0:** tf_build_regulons() + tf_activity() (Steps 2-3). Adds activity scoring.
- **v0.3.0:** tf_differential() + plotting. Full pipeline complete.
