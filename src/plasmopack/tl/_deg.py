"""Differential gene expression with user-selectable methods.

The user chooses the statistical method via ``method=``. The choice reflects a
real trade-off documented in the single-cell DE literature:

- **Cell-level tests** (``wilcoxon``, ``t-test``, ``logreg``) treat every cell
  as an observation. Fast and fine for exploration, but they suffer from
  *pseudo-replication*: with few biological samples they overstate significance
  because cells from one sample are not independent.
- **Pseudobulk / replicate-aware** methods (``pseudobulk``, ``deseq2``)
  aggregate counts per biological replicate first, then test across replicates.
  This is the recommended approach for multi-sample designs and the one that
  controls false positives (Squair et al. 2021; sc-best-practices).

Results are written to ``adata.uns["deg"][groupby][method]`` as a tidy
DataFrame and provenance is stamped. Use :func:`list_deg_methods` to see which
methods are available in the current install.
"""

from __future__ import annotations

import importlib.util
import warnings
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from plasmopack._utils.provenance import record_provenance

_EPS = 1e-9


# ---------------------------------------------------------------------------
# Method registry (what exists + what it needs)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class DEGMethod:
    name: str
    level: str  # "cell" or "pseudobulk"
    extra: str | None  # pip extra needed, or None for core
    needs_module: str | None  # importable module that must be present
    needs_replicates: bool
    expects: str  # what input data the method expects
    description: str


# What each data level expects as input (so users pick the right layer).
_CELL_EXPECTS = (
    "log-normalized expression (e.g. run pp.normalize, uses layer='lognorm')"
)
_PB_EXPECTS = (
    "raw integer counts (uses layer='counts'); needs sample_key for replicates"
)

_METHODS: dict[str, DEGMethod] = {
    "wilcoxon": DEGMethod(
        name="wilcoxon",
        level="cell",
        extra=None,
        needs_module=None,
        needs_replicates=False,
        expects=_CELL_EXPECTS,
        description="Mann-Whitney/Wilcoxon rank-sum per gene (non-parametric, core).",
    ),
    "t-test": DEGMethod(
        name="t-test",
        level="cell",
        extra=None,
        needs_module=None,
        needs_replicates=False,
        expects=_CELL_EXPECTS,
        description="Welch's t-test per gene (core).",
    ),
    "pseudobulk": DEGMethod(
        name="pseudobulk",
        level="pseudobulk",
        extra=None,
        needs_module=None,
        needs_replicates=True,
        expects=_PB_EXPECTS,
        description="Aggregate counts per replicate, then Welch's t-test across "
        "replicates (replicate-aware, core). Requires sample_key.",
    ),
    "deseq2": DEGMethod(
        name="deseq2",
        level="pseudobulk",
        extra="deg",
        needs_module="pydeseq2",
        needs_replicates=True,
        expects=_PB_EXPECTS,
        description="Pseudobulk negative-binomial GLM via pydeseq2 (gold standard "
        "for multi-sample designs). Requires sample_key and plasmopack[deg].",
    ),
    "logreg": DEGMethod(
        name="logreg",
        level="cell",
        extra="ml",
        needs_module="sklearn",
        needs_replicates=False,
        expects=_CELL_EXPECTS,
        description="Logistic regression per gene via scikit-learn. "
        "Requires plasmopack[ml].",
    ),
    "mast": DEGMethod(
        name="mast",
        level="cell",
        extra=None,
        needs_module="rpy2",
        needs_replicates=False,
        expects=_CELL_EXPECTS,
        description="MAST hurdle model (zero-inflation aware). Needs R + MAST "
        "via rpy2; not yet implemented.",
    ),
}


def _available(m: DEGMethod) -> bool:
    if m.name == "mast":
        return False  # not implemented yet
    if m.needs_module is None:
        return True
    return importlib.util.find_spec(m.needs_module) is not None


def list_deg_methods() -> pd.DataFrame:
    """Return a table of DEG methods and whether each is available to use."""
    rows = []
    for m in _METHODS.values():
        rows.append(
            {
                "method": m.name,
                "level": m.level,
                "available": _available(m),
                "extra": m.extra or "(core)",
                "needs_replicates": m.needs_replicates,
                "expects": m.expects,
                "description": m.description,
            }
        )
    return pd.DataFrame(rows).set_index("method")


def _looks_like_counts(x: Any) -> bool:
    """Heuristic: does this matrix look like raw integer counts?"""
    sample = _to_dense(x[: min(50, x.shape[0])])
    if sample.size == 0:
        return False
    finite = sample[np.isfinite(sample)]
    if finite.size == 0:
        return False
    nonneg = bool((finite >= 0).all())
    # log-normalized data is non-integer; all-integer non-negative values with
    # some magnitude are the signature of raw counts. The integer test is the
    # reliable discriminator; max>1 just rules out boolean/degenerate matrices.
    integral = bool(np.allclose(finite, np.round(finite)))
    return nonneg and integral and float(finite.max()) > 1.0


def _warn_input_mismatch(x: Any, spec: DEGMethod, layer: str | None) -> None:
    """Warn if the chosen data looks wrong for the method's level."""
    is_counts = _looks_like_counts(x)
    where = f"layer={layer!r}" if layer is not None else "adata.X"
    if spec.level == "cell" and is_counts:
        warnings.warn(
            f"method={spec.name!r} expects {spec.expects}, but {where} looks "
            f"like raw counts. Run pp.normalize(adata) first and/or pass "
            f"layer='lognorm'. Proceeding, but results may be unreliable.",
            stacklevel=3,
        )
    elif spec.level == "pseudobulk" and not is_counts:
        warnings.warn(
            f"method={spec.name!r} expects {spec.expects}, but {where} does not "
            f"look like raw counts. Pass layer='counts' (raw). Proceeding, but "
            f"results may be unreliable.",
            stacklevel=3,
        )


# ---------------------------------------------------------------------------
# Small numeric helpers
# ---------------------------------------------------------------------------
def _to_dense(x: Any) -> np.ndarray:
    return np.asarray(x.toarray() if sp.issparse(x) else x, dtype=np.float64)


def _bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment."""
    p = np.asarray(pvals, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out


def _pct_nonzero(mat: np.ndarray) -> np.ndarray:
    return (mat > 0).mean(axis=0) if mat.shape[0] else np.zeros(mat.shape[1])


def _log2fc_from_log(mean_g: np.ndarray, mean_r: np.ndarray) -> np.ndarray:
    """log2 fold change from log1p-space means (un-log, then ratio)."""
    return np.log2((np.expm1(mean_g) + _EPS) / (np.expm1(mean_r) + _EPS))


def _log2fc_linear(mean_g: np.ndarray, mean_r: np.ndarray) -> np.ndarray:
    return np.log2((mean_g + _EPS) / (mean_r + _EPS))


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def deg(
    adata: AnnData,
    groupby: str,
    *,
    method: str = "wilcoxon",
    groups: list[str] | None = None,
    reference: str = "rest",
    layer: str | None = None,
    sample_key: str | None = None,
) -> AnnData:
    """Differential expression, one group vs a reference, with a chosen method.

    Parameters
    ----------
    adata
        Annotated data; result written to ``adata.uns["deg"]``.
    groupby
        A categorical column in ``adata.obs`` defining the groups to compare.
    method
        One of :func:`list_deg_methods` (default ``"wilcoxon"``). Optional
        methods raise a clear error naming the pip extra to install.
    groups
        Which groups to test (default: all categories in ``groupby``).
    reference
        ``"rest"`` (each group vs all other cells) or a specific group name.
    layer
        Expression layer. Cell-level methods default to ``"lognorm"`` if
        present else ``X`` (expect log-normalized data). Pseudobulk methods
        default to ``"counts"`` if present else ``X`` (expect raw counts).
    sample_key
        Biological-replicate column in ``adata.obs``. **Required** for
        ``pseudobulk`` and ``deseq2``.

    Returns
    -------
    AnnData
        The same object; results in ``adata.uns["deg"][groupby][method]``.
    """
    if method not in _METHODS:
        raise ValueError(
            f"Unknown DEG method {method!r}. Available: "
            f"{', '.join(_METHODS)}. See plasmopack.tl.list_deg_methods()."
        )
    spec = _METHODS[method]
    if not _available(spec):
        if method == "mast":
            raise NotImplementedError(
                "method='mast' is not implemented yet (needs R + MAST via "
                "rpy2). Use 'wilcoxon' or, for replicate-aware testing, "
                "'pseudobulk'/'deseq2'."
            )
        raise ImportError(
            f"method={method!r} needs an optional dependency. Install it with: "
            f"pip install 'plasmopack[{spec.extra}]'."
        )
    if spec.needs_replicates and sample_key is None:
        raise ValueError(
            f"method={method!r} is replicate-aware and requires "
            f"sample_key=<obs column identifying biological replicates>."
        )

    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs.")
    labels = adata.obs[groupby].astype(str).to_numpy()
    all_groups = list(dict.fromkeys(labels))
    test_groups = groups if groups is not None else all_groups

    # choose default layer per method level (documented in list_deg_methods)
    if layer is None:
        if spec.level == "pseudobulk":
            layer = "counts" if "counts" in adata.layers else None
        else:
            layer = "lognorm" if "lognorm" in adata.layers else None
    if layer is not None and layer not in adata.layers:
        raise KeyError(
            f"layer={layer!r} not found in adata.layers "
            f"(available: {list(adata.layers)}). method={method!r} expects "
            f"{spec.expects}."
        )
    x = adata.layers[layer] if layer is not None else adata.X

    # help the user: warn if the data looks wrong for this method's level
    _warn_input_mismatch(x, spec, layer)

    var_names = np.asarray(adata.var_names)
    frames: list[pd.DataFrame] = []
    for g in test_groups:
        in_g = labels == g
        in_r = ~in_g if reference == "rest" else labels == reference
        if in_g.sum() == 0 or in_r.sum() == 0:
            continue

        if method in ("wilcoxon", "t-test"):
            res = _cell_test(x, in_g, in_r, method)
        elif method == "logreg":
            res = _logreg(x, in_g, in_r)
        elif method == "pseudobulk":
            res = _pseudobulk_ttest(
                x, in_g, in_r, adata.obs[sample_key].astype(str).to_numpy()
            )
        elif method == "deseq2":
            res = _deseq2(x, in_g, in_r, adata.obs[sample_key].astype(str).to_numpy())
        else:  # pragma: no cover - guarded above
            raise AssertionError(method)

        logfc, pval, pct_g, pct_r = res
        frame = pd.DataFrame(
            {
                "group": g,
                "gene": var_names,
                "logFC": logfc,
                "pval": pval,
                "padj": _bh(pval),
                "pct_grp": pct_g,
                "pct_rest": pct_r,
                "method": method,
            }
        )
        frames.append(frame)

    result = (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(
            columns=[
                "group",
                "gene",
                "logFC",
                "pval",
                "padj",
                "pct_grp",
                "pct_rest",
                "method",
            ]
        )
    )

    store = adata.uns.setdefault("deg", {})
    if not isinstance(store, dict):
        store = {}
        adata.uns["deg"] = store
    by = store.setdefault(groupby, {})
    by[method] = result

    record_provenance(
        adata,
        "tl.deg",
        params={
            "groupby": groupby,
            "method": method,
            "reference": reference,
            "layer": layer,
            "sample_key": sample_key,
            "n_groups": len(frames),
        },
    )
    return adata


# ---------------------------------------------------------------------------
# Per-method implementations
# ---------------------------------------------------------------------------
def _cell_test(
    x: Any, in_g: np.ndarray, in_r: np.ndarray, method: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    from scipy import stats

    xg = _to_dense(x[in_g])
    xr = _to_dense(x[in_r])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if method == "wilcoxon":
            stat = stats.mannwhitneyu(xg, xr, axis=0, alternative="two-sided")
        else:
            stat = stats.ttest_ind(xg, xr, axis=0, equal_var=False)
        pval = np.asarray(stat.pvalue, dtype=float)
    pval = np.where(np.isfinite(pval), pval, 1.0)
    logfc = _log2fc_from_log(xg.mean(axis=0), xr.mean(axis=0))
    return logfc, pval, _pct_nonzero(xg), _pct_nonzero(xr)


def _logreg(
    x: Any, in_g: np.ndarray, in_r: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Per-gene univariate logistic regression Wald-style p-value approximation
    # via the likelihood-ratio is heavy; we use sklearn coefficients + a simple
    # permutation-free score. For a first version we fit one multivariate model
    # and report per-gene coefficients with a normal-approx p-value.
    from scipy import stats
    from sklearn.linear_model import LogisticRegression

    xg = _to_dense(x[in_g])
    xr = _to_dense(x[in_r])
    xx = np.vstack([xg, xr])
    yy = np.concatenate([np.ones(xg.shape[0]), np.zeros(xr.shape[0])])
    model = LogisticRegression(max_iter=1000, C=1.0)
    model.fit(xx, yy)
    coef = model.coef_.ravel()
    # normal approximation around coefficient spread
    z = coef / (np.std(coef) + _EPS)
    pval = 2 * stats.norm.sf(np.abs(z))
    logfc = _log2fc_from_log(xg.mean(axis=0), xr.mean(axis=0))
    return logfc, pval, _pct_nonzero(xg), _pct_nonzero(xr)


def _aggregate(
    x: Any, mask: np.ndarray, samples: np.ndarray
) -> tuple[np.ndarray, list[str]]:
    """Sum counts per sample within a cell mask -> (n_samples x n_genes)."""
    xs = _to_dense(x[mask])
    s = samples[mask]
    uniq = list(dict.fromkeys(s))
    out = (
        np.vstack([xs[s == u].sum(axis=0) for u in uniq])
        if uniq
        else np.empty((0, xs.shape[1]))
    )
    return out, uniq


def _cpm(mat: np.ndarray) -> np.ndarray:
    lib = mat.sum(axis=1, keepdims=True)
    lib = np.where(lib > 0, lib, 1.0)
    return np.asarray(mat / lib * 1e6, dtype=np.float64)


def _pseudobulk_ttest(
    x: Any, in_g: np.ndarray, in_r: np.ndarray, samples: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    from scipy import stats

    pb_g, sg = _aggregate(x, in_g, samples)
    pb_r, sr = _aggregate(x, in_r, samples)
    if pb_g.shape[0] < 2 or pb_r.shape[0] < 2:
        raise ValueError(
            "pseudobulk needs at least 2 replicates per side; got "
            f"{pb_g.shape[0]} vs {pb_r.shape[0]}. Check sample_key."
        )
    cpm_g = np.log2(_cpm(pb_g) + 1.0)
    cpm_r = np.log2(_cpm(pb_r) + 1.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stat = stats.ttest_ind(cpm_g, cpm_r, axis=0, equal_var=False)
        pval = np.asarray(stat.pvalue, dtype=float)
    pval = np.where(np.isfinite(pval), pval, 1.0)
    logfc = cpm_g.mean(axis=0) - cpm_r.mean(axis=0)
    # pct across original cells
    return (
        logfc,
        pval,
        _pct_nonzero(_to_dense(x[in_g])),
        _pct_nonzero(_to_dense(x[in_r])),
    )


def _deseq2(
    x: Any, in_g: np.ndarray, in_r: np.ndarray, samples: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    pb_g, sg = _aggregate(x, in_g, samples)
    pb_r, sr = _aggregate(x, in_r, samples)
    if pb_g.shape[0] < 2 or pb_r.shape[0] < 2:
        raise ValueError(
            "deseq2 needs at least 2 replicates per side; got "
            f"{pb_g.shape[0]} vs {pb_r.shape[0]}. Check sample_key."
        )
    counts = pd.DataFrame(
        np.vstack([pb_g, pb_r]).astype(int),
        index=[f"g_{s}" for s in sg] + [f"r_{s}" for s in sr],
    )
    meta = pd.DataFrame(
        {"condition": ["group"] * len(sg) + ["rest"] * len(sr)},
        index=counts.index,
    )
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "group", "rest"], quiet=True)
    stat.summary()
    res = stat.results_df
    logfc = res["log2FoldChange"].to_numpy()
    pval = np.where(
        np.isfinite(res["pvalue"].to_numpy()), res["pvalue"].to_numpy(), 1.0
    )
    return (
        logfc,
        pval,
        _pct_nonzero(_to_dense(x[in_g])),
        _pct_nonzero(_to_dense(x[in_r])),
    )
