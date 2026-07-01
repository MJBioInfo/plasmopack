# Cell-Cell Interaction (CCI) Strategy — Design Decision

**Date:** 2026-06-23
**Status:** Decided — v0.2 feature (not v0.1)
**Author:** Majeed Jamakhani + Claude

---

## Why This Exists — The Gap

Cell-cell interaction (CCI) analysis is well-established for mammalian single-cell data.
Tools like CellPhoneDB, LIANA, CellChat, NicheNet are widely used. They work by scoring
co-expression of ligand-receptor (LR) pairs across cell types.

**None of these tools work properly for Plasmodium because:**

1. Their LR databases contain only human/mouse proteins. Plasmodium proteins are absent.
2. They cannot model parasite-to-host interactions — only cell-to-cell within one species.
3. The parasite actively exports proteins into the host cell (a unique biology with no
   mammalian equivalent), which no CCI framework currently models.

**The gap this fills:**
- No existing tool scores known parasite-host molecular interactions in scRNA-seq data
- No structured, versioned, AnnData-compatible parasite-host LR database exists
- The Malaria Cell Atlas and emerging parasite+host co-profiling datasets have no
  dedicated analytical tool for interaction scoring

**Important honest statement:** We are building the FIRST version of this. The field does
not have this tool yet. That is the contribution. We are not claiming completeness — we
are claiming to be first and transparent.

---

## What "Cell-Cell Interaction" Means in Plasmodium Context

Unlike mammalian CCI (e.g. macrophage→T cell via cytokine), Plasmodium CCI involves
three biologically distinct interaction types:

### Type 1 — Parasite invasion of host cell
The parasite surface proteins bind to RBC (or hepatocyte) surface receptors and enable
entry. These are the best-characterised Plasmodium-host molecular interactions.

| Parasite ligand | Host receptor | Stage | Evidence |
|---|---|---|---|
| PfEBA-175 | Glycophorin A (GYPA) | merozoite | Experimental — high confidence |
| PfRH5 | Basigin (BSG/CD147) | merozoite | Experimental — crystal structure |
| AMA1 | RON2 (parasite-derived, in host membrane) | merozoite | Experimental |
| PfRH1 | Complement receptor 1 (CR1) | merozoite | Experimental |
| MSP1 | Band 3 (SLC4A1) | merozoite | Experimental |
| PbCIRP | HSP70 (host surface) | liver stage | Experimental (Pb) |

### Type 2 — Infected RBC to host endothelium/immune cell (cytoadherence)
The parasite exports proteins that appear on the surface of the infected RBC, which then
bind to blood vessel walls or immune cells. This causes cerebral malaria, severe disease.

| Exported protein | Host receptor | Disease relevance |
|---|---|---|
| PfEMP1 (var gene family) | ICAM-1, PECAM-1, CD36, EPCR | Cytoadherence, cerebral malaria |
| RIFIN | LILRB1 (on NK/T cells) | Immune evasion |
| STEVOR | Unknown | RBC deformability |
| PfMC-2TM | Unknown | Knob formation |

### Type 3 — Host immune response to infected RBC
Standard mammalian CCI (macrophage, neutrophil, dendritic cell, T cell responses).
Existing tools (LIANA, CellPhoneDB) handle this side. We complement them by adding the
parasite side and providing a unified analysis framework.

---

## The Parasite-Host LR Database

This is the central artifact of pp.tl.cci. It does not exist elsewhere.
We build it by merging four sources:

### Source 1 — IntAct (EBI Molecular Interaction Database)
- URL: ebi.ac.uk/intact
- What it has: experimentally validated molecular interactions including cross-species;
  Plasmodium entries curated from primary literature with PSI-MI evidence codes
- Access: REST API (stable, EBI-maintained)
- Filter we apply: only interactions with experimental evidence codes
  (two-hybrid, pull-down, co-immunoprecipitation, surface plasmon resonance)
  — exclude predicted and text-mined entries

### Source 2 — PHI-base (Pathogen-Host Interaction Database)
- URL: phi-base.org
- What it has: curated pathogen-host interactions classified by phenotype outcome
  (loss of pathogenicity, reduced virulence, etc.); has Plasmodium entries
- Access: bulk download only (no REST API) → we ship versioned parquet snapshot
- Quality: manually curated, peer-reviewed, respected in the pathogen community

### Source 3 — BioGRID
- URL: thebiogrid.org
- What it has: biological interactions including intra-Plasmodium PPI; some
  cross-species Plasmodium-human interaction entries
- Access: REST API (requires free registration key)
- Note: BioGRID is stronger for intra-Plasmodium PPI than cross-species.
  Cross-species entries are fewer but real.

### Source 4 — CellPhoneDB / LIANA database (host LR pairs)
- What it has: curated human ligand-receptor pairs, complexes, co-factors
- Access: static download (versioned)
- Use: provides the host-to-host and immune-cell interaction pairs that complement
  the parasite-specific pairs from Sources 1-3

### Merging and Evidence Tiers

```
All four sources → merge on (parasite_gene, host_gene) → deduplicate

Evidence tier assigned:
  "experimental"      → direct protein interaction, primary assay, published paper
                        (from IntAct PSI-MI experimental codes, PHI-base, BioGRID experimental)
  "text_mined"        → mentioned in papers but not directly assayed
  "ortholog_inferred" → interaction known for ortholog species, mapped to Plasmodium
  "predicted"         → computational prediction only

Default filter: only "experimental" tier shown unless user sets evidence_level="all"
```

### What the Database Actually Contains (Honest Size Estimate)

The Plasmodium-host experimental interaction database is small:
- **Well-characterised invasion pairs:** ~15–20 (EBA, RH, AMA1, MSP families + receptors)
- **Cytoadherence pairs (PfEMP1/var genes):** ~10–15 well-validated receptor-ligand pairs
- **Immune evasion pairs (RIFIN/STEVOR):** ~5–8
- **Liver stage invasion (Pb, Pf sporozoites):** ~5–10
- **Total experimental tier:** approximately 50–80 unique parasite-host pairs

This is small compared to the thousands of pairs in human CCI databases. We are explicit
about this. The contribution is organising these known pairs into a queryable, versioned,
AnnData-compatible database — not claiming to have a complete LR atlas.

Intra-host immune pairs (from CellPhoneDB/LIANA): thousands — well-covered.

---

## The Pipeline: pp.tl.cci.*

### Step 1 — Build/Load the LR Database

```python
pp.tl.cci.load_lr_database(
    organism="pfalciparum",    # or "pberghei"
    host="human",              # or "mouse" for Pb liver stage studies
    evidence_level="experimental",  # default; "all" to include predicted
    version="2024-01"          # pinned version of the bundled snapshot
)
# → adata.uns["cci"]["lr_database"] = DataFrame:
#   parasite_gene, host_gene, interaction_type, evidence_tier, pmid, source_db, notes
```

### Step 2 — Score Interactions Across Cell Types

```python
pp.tl.cci.score_interactions(
    adata,
    groupby="leiden",        # cell type / cluster column in obs
    parasite_key="species",  # column in obs identifying parasite vs host cells
    n_permutations=1000      # for p-value estimation
)
# Method: for each LR pair in the database, and each (parasite_cluster, host_cluster) pair:
#   score = mean_expression(ligand in parasite_cluster) × mean_expression(receptor in host_cluster)
#   p-value by permutation test (shuffle cell labels, recompute score)
#   — same method as CellPhoneDB, well-established
#
# → adata.uns["cci"]["scores"] = DataFrame:
#   ligand, receptor, parasite_cluster, host_cluster, score, pval, padj,
#   evidence_tier, interaction_type
```

### Step 3 — Var Gene / PfEMP1 Activity (Plasmodium-specific)

```python
pp.tl.cci.var_gene_activity(adata)
# P. falciparum has ~60 var genes (PfEMP1 variants) — only 1-2 expressed per cell
# This function:
#   → identifies which var genes are expressed in each cell cluster
#   → predicts which host receptors may be engaged based on the LR database
#   → writes adata.obs["dominant_var_gene"]
#   → writes adata.uns["cci"]["var_gene_activity"]
```

---

## Relationship to Existing CCI Tools

We do not replace CellPhoneDB or LIANA. We extend them:

```
CellPhoneDB / LIANA
  → excellent for host-only immune interactions (human-human, mouse-mouse)
  → no parasite LR pairs

plasmopack pp.tl.cci
  → adds parasite-host LR pairs from IntAct, PHI-base, BioGRID
  → uses same scoring method as CellPhoneDB (compatible, comparable)
  → AnnData-native — results stay in AnnData
  → versioned database — reproducible

Future integration path: export pp.tl.cci results in CellPhoneDB-compatible format
so users can visualise in existing tools.
```

---

## How We Compare to plasmoRUtils

| Capability | plasmoRUtils | plasmopack |
|---|---|---|
| Any CCI functionality | ✗ None | ✓ First Plasmodium-specific CCI tool |
| Parasite-host LR database | ✗ | ✓ from IntAct + PHI-base + BioGRID |
| Interaction scoring | ✗ | ✓ permutation-based, CellPhoneDB method |
| Var gene / PfEMP1 analysis | ✗ | ✓ |
| AnnData integration | ✗ | ✓ |
| Evidence tiers | ✗ | ✓ experimental vs predicted |

---

## Honest Limitations for Reviewers

1. **The parasite-host experimental LR database is small (~50-80 pairs).** This reflects
   the state of the field, not a gap in our approach. We document exactly which pairs are
   included and which papers they come from.

2. **Scoring requires single-cell data that captures both parasite and host cells.**
   Most existing Plasmodium scRNA-seq datasets profile only the parasite (Malaria Cell
   Atlas) or only the host. Dual-species profiling datasets exist but are fewer. We
   document which datasets are suitable and provide example data via pp.datasets.

3. **Permutation-based p-values are sensitive to cell number per cluster.** Small
   clusters (< 20 cells) produce unreliable p-values. We warn the user explicitly
   when clusters are too small and recommend minimum cell counts.

4. **Var gene expression is low and noisy in scRNA-seq.** PfEMP1/var genes are
   transcribed at low levels and are difficult to detect in single cells. We
   report raw expression values with appropriate caveats.

5. **We do not do spatial CCI.** Tools like Squidpy model physical proximity.
   Our scoring assumes cell types are in the same biological compartment. Physical
   proximity data is out of scope for v0.1 of this module.

---

## AnnData Storage Convention (Complete)

```
adata.uns
  ├── "cci"
  │     ├── "lr_database"       DataFrame — the full LR pair table with evidence
  │     ├── "scores"            DataFrame — (ligand, receptor, clusters, score, pval, padj, tier)
  │     ├── "var_gene_activity" DataFrame — var gene expression per cluster
  │     └── "metadata"         dict — database version, organism, host, date built
  └── "plasmopack"
        └── "history"          records LR database version + scoring parameters

adata.obs
  └── "dominant_var_gene"      str or NaN (from var_gene_activity step)

adata.var
  └── "is_parasite_lr_ligand"  bool — is this gene a ligand in our LR database
```

---

## Version / Release Plan

- **v0.1.0:** NOT included. Too sparse and complex to validate without dual-species data.
- **v0.2.0:** load_lr_database() only — users can inspect and query the LR database.
              No scoring yet. This alone is useful (lookup tool) and safe to release.
- **v0.3.0:** score_interactions() — full CCI scoring with permutation test.
- **v0.4.0:** var_gene_activity() + plotting (pp.pl.cci_dotplot, pp.pl.cci_network).
- **Future:** Integration with spatial transcriptomics (Squidpy compatibility).

**Why delay CCI to v0.2+:** The LR database needs to be validated carefully. Releasing
an unvalidated database and presenting it as authoritative would be scientifically
irresponsible and a reviewer risk. Better to release the database in v0.2 with clear
documentation and let the community validate it before adding automated scoring.

---

## Paper Talking Points

1. **Novel contribution:** First structured, versioned, AnnData-native parasite-host
   interaction database and scoring tool for Plasmodium. Nothing comparable exists.

2. **Conservative approach is a strength:** By restricting the default database to
   experimentally validated pairs (IntAct experimental codes, PHI-base curated entries),
   we trade completeness for accuracy. This is the correct scientific choice.

3. **Complementary to existing tools:** We use the same scoring method as CellPhoneDB
   (mean expression × permutation p-value), making our results directly comparable to
   host-only CCI analyses in the same dataset.

4. **The field is ready for this:** Dual-species Plasmodium+host scRNA-seq datasets now
   exist (Chua et al. 2023, Poran et al. 2017) and the Malaria Cell Atlas provides a
   reference for parasite cell states. The data infrastructure exists; the analysis tool
   is the missing piece.

5. **Transparent about gaps:** We clearly list which interaction types are covered and
   which are not. Reviewers appreciate honesty over overclaiming.

---

## Validation Strategy

Validate against published biology, not against another tool (no equivalent tool exists):

1. **Load PSA_minimal.h5ad** (P. berghei liver stage) — check that known Pb liver-stage
   invasion proteins (PbCIRP, TRAP, CSP) appear in the database and are detected as
   expressed at the relevant timepoints.

2. **Use Poran et al. 2017 dataset** (P. berghei + mouse liver cells, GEO) — score
   interactions and check that known Pb sporozoite→hepatocyte interactions are recovered.

3. **Correlation with known biology:** At 24hpi (early liver stage), invasion-related LR
   pairs should score high. At 62hpi (late liver stage / merozoite formation), egress-
   related proteins should appear. This is the biological validation.
