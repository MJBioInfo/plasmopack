"""Circos-style pathway <-> gene chord diagram.

Draws enriched pathways as arcs on one side of a circle and their member genes
as arcs on the other, linked by ribbons. Genes are coloured by fold change
(up vs down), so one picture shows *which pathways are enriched* and *which
direction their genes move*.

Built on matplotlib only (no circos dependency). The user controls how many
pathways and genes to show.
"""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path

_R_OUTER = 1.0
_R_INNER = 0.90  # arcs are a thin ring between inner and outer
_R_RIBBON = 0.88  # ribbons start just inside the ring
_RGBA = tuple[float, float, float, float]
_PATHWAY_COLORS = [plt.get_cmap("tab10")(i) for i in range(10)]


def _polar(theta_deg: float, r: float) -> tuple[float, float]:
    t = np.deg2rad(theta_deg)
    return r * np.cos(t), r * np.sin(t)


def _ribbon(a0: float, a1: float, b0: float, b1: float, color: _RGBA) -> PathPatch:
    """A filled ribbon between arc span [a0,a1] and [b0,b1] (degrees)."""
    c = 0.25  # curvature: pull control points toward the centre
    pa0, pa1 = _polar(a0, _R_RIBBON), _polar(a1, _R_RIBBON)
    pb0, pb1 = _polar(b0, _R_RIBBON), _polar(b1, _R_RIBBON)
    verts = [
        pa0,
        (pa0[0] * c, pa0[1] * c),
        (pb1[0] * c, pb1[1] * c),
        pb1,
        pb0,  # straight chord across the (small) gene/pathway arc
        (pb0[0] * c, pb0[1] * c),
        (pa1[0] * c, pa1[1] * c),
        pa1,
        pa0,
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    return PathPatch(Path(verts, codes), facecolor=color, edgecolor="none", zorder=1)


def _gene_color(lfc: float, vmax: float, up_color: str, down_color: str) -> _RGBA:
    base = up_color if lfc >= 0 else down_color
    alpha = 0.35 + 0.65 * min(abs(lfc) / vmax, 1.0) if vmax > 0 else 0.6
    return to_rgba(base, alpha)


def _label_map(adata: AnnData, gene_label: str | None) -> dict[str, str]:
    """gene id -> display label (symbol where available, else the id).

    ``gene_label`` names a ``var`` column; when ``None``, ``"symbol"`` is used
    if present (e.g. from ``pp.pp.add_gene_symbols``), otherwise the gene id.
    """
    import numpy as np

    col = gene_label or ("symbol" if "symbol" in adata.var.columns else None)
    ids = [str(v) for v in adata.var_names]
    if col is None:
        return {gid: gid for gid in ids}
    if col not in adata.var.columns:
        raise KeyError(f"gene_label={col!r} not found in adata.var.")
    out: dict[str, str] = {}
    for gid, val in zip(ids, adata.var[col], strict=True):
        s = (
            ""
            if val is None or (isinstance(val, float) and np.isnan(val))
            else str(val).strip()
        )
        out[gid] = s if s and s.lower() not in {"nan", "none", "na", "-"} else gid
    return out


def pathway_chord(
    adata: AnnData,
    enrichment_key: str,
    *,
    groupby: str,
    group: str,
    deg_method: str = "wilcoxon",
    n_pathways: int = 3,
    n_genes: int = 20,
    gene_rank: str = "logfc",
    gene_label: str | None = None,
    up_color: str = "#d62728",
    down_color: str = "#1f77b4",
    figsize: tuple[float, float] = (8.0, 8.0),
    title: str | None = None,
    ax: Axes | None = None,
) -> Figure:
    """Draw a pathway <-> gene chord (circos-style) diagram.

    Parameters
    ----------
    adata
        Must contain enrichment results in ``adata.uns["enrichment"]
        [enrichment_key]`` (from :func:`plasmopack.tl.enrich`) and DEG results
        in ``adata.uns["deg"][groupby][deg_method]`` (for gene fold changes).
    enrichment_key
        Key of the enrichment result to plot.
    groupby, group, deg_method
        Locate the DEG table that supplies each gene's ``logFC``.
    n_pathways
        How many top pathways (by padj) to show.
    n_genes
        How many genes to show (top-ranked across the shown pathways).
    gene_rank
        Rank genes by ``"logfc"`` (|log fold change|), ``"padj"``, or
        ``"combined"`` (``-log10(padj) * |logFC|``).
    up_color, down_color
        Colours for up- and down-regulated genes.
    figsize, title, ax
        Standard matplotlib options; ``ax`` must be a square axes if given.

    Returns
    -------
    matplotlib.figure.Figure
        The figure (also drawn into ``ax`` if provided).
    """
    enr = adata.uns.get("enrichment", {}).get(enrichment_key)
    if enr is None or len(enr) == 0:
        raise KeyError(
            f"No enrichment result at adata.uns['enrichment'][{enrichment_key!r}]."
        )
    deg = adata.uns.get("deg", {}).get(groupby, {}).get(deg_method)
    if deg is None:
        raise KeyError(
            f"No DEG result at adata.uns['deg'][{groupby!r}][{deg_method!r}]."
        )
    deg_g = deg[deg["group"].astype(str) == str(group)]
    lfc_map = dict(zip(deg_g["gene"].astype(str), deg_g["logFC"], strict=False))
    padj_map = dict(zip(deg_g["gene"].astype(str), deg_g["padj"], strict=False))

    # --- choose pathways and their member (overlap) genes ---
    top = enr.sort_values("padj", kind="stable").head(n_pathways)
    pathways = []
    for _, row in top.iterrows():
        members = [g for g in str(row["overlap_genes"]).split(",") if g]
        pathways.append((str(row["name"] or row["set_id"]), members))

    # candidate genes = union, ranked, top n_genes
    cand = {g for _, members in pathways for g in members if g in lfc_map}

    def _score(g: str) -> float:
        lfc = abs(float(lfc_map.get(g, 0.0)))
        nlp = -float(np.log10(max(float(padj_map.get(g, 1.0)), 1e-300)))
        return float(
            {"logfc": lfc, "padj": nlp, "combined": lfc * nlp}.get(gene_rank, lfc)
        )

    genes = sorted(cand, key=_score, reverse=True)[:n_genes]
    if not genes:
        raise ValueError("No genes to plot (check enrichment/DEG overlap).")
    gene_set = set(genes)

    links = [
        (pi, g)
        for pi, (_, members) in enumerate(pathways)
        for g in members
        if g in gene_set
    ]

    # --- geometry: pathways on left arc, genes on right arc ---
    gap = 8.0
    # genes on the right half: from -90+gap up to 90-gap
    g_lo, g_hi = -90 + gap, 90 - gap
    gene_edges = np.linspace(g_lo, g_hi, len(genes) + 1)
    gene_span = {g: (gene_edges[i], gene_edges[i + 1]) for i, g in enumerate(genes)}
    gene_mid = {g: 0.5 * (a + b) for g, (a, b) in gene_span.items()}
    # pathways on the left half: from 90+gap to 270-gap
    p_lo, p_hi = 90 + gap, 270 - gap
    # pathway arc proportional to number of shown links
    counts = np.array(
        [sum(1 for lp, _ in links if lp == pi) for pi in range(len(pathways))]
    )
    counts = np.where(counts > 0, counts, 1)
    widths = (p_hi - p_lo) * counts / counts.sum()
    p_edges = p_lo + np.concatenate([[0], np.cumsum(widths)])
    pathway_span = {pi: (p_edges[pi], p_edges[pi + 1]) for pi in range(len(pathways))}

    # sub-divide each pathway arc among its genes (ordered by gene angle)
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    else:
        parent = ax.get_figure()
        fig = parent if isinstance(parent, Figure) else plt.figure(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.55, 1.55)

    vmax = max((abs(float(lfc_map.get(g, 0.0))) for g in genes), default=1.0)
    labels = _label_map(adata, gene_label)  # gene id -> symbol-or-id

    def _radial(mid: float) -> tuple[float, str]:
        """Rotation + horizontal alignment so a label reads outward, upright."""
        m = mid % 360
        if 90 < m < 270:  # left half -> flip so text isn't upside-down
            return mid - 180, "right"
        return mid, "left"

    # draw gene wedges + labels (coloured by fold change)
    for g in genes:
        a0, a1 = gene_span[g]
        col = _gene_color(float(lfc_map.get(g, 0.0)), vmax, up_color, down_color)
        ax.add_patch(
            Wedge(
                (0, 0),
                _R_OUTER,
                a0,
                a1,
                width=_R_OUTER - _R_INNER,
                facecolor=col,
                edgecolor="white",
                linewidth=0.5,
                zorder=3,
            )
        )
        mid = gene_mid[g]
        rot, ha = _radial(mid)
        lx, ly = _polar(mid, _R_OUTER + 0.03)
        text = labels.get(g, g)
        label = text if len(text) <= 16 else text[:15] + "…"
        ax.text(
            lx,
            ly,
            label,
            rotation=rot,
            rotation_mode="anchor",
            ha=ha,
            va="center",
            fontsize=6.5,
        )

    # draw pathway wedges + labels
    for pi, (name, _) in enumerate(pathways):
        a0, a1 = pathway_span[pi]
        ax.add_patch(
            Wedge(
                (0, 0),
                _R_OUTER,
                a0,
                a1,
                width=_R_OUTER - _R_INNER,
                facecolor=_PATHWAY_COLORS[pi % len(_PATHWAY_COLORS)],
                edgecolor="white",
                linewidth=0.5,
                zorder=3,
            )
        )
        mid = 0.5 * (a0 + a1)
        rot, ha = _radial(mid)
        lx, ly = _polar(mid, _R_OUTER + 0.05)
        disp = name if len(name) <= 28 else name[:27] + "…"
        ax.text(
            lx,
            ly,
            disp,
            rotation=rot,
            rotation_mode="anchor",
            ha=ha,
            va="center",
            fontsize=8.5,
            fontweight="bold",
        )

    # draw ribbons: subdivide each pathway arc among its links and connect
    for pi in range(len(pathways)):
        p_links = [g for lp, g in links if lp == pi]
        p_links.sort(key=lambda g: gene_mid[g])  # order to reduce crossings
        if not p_links:
            continue
        a0, a1 = pathway_span[pi]
        sub = np.linspace(a0, a1, len(p_links) + 1)
        for j, g in enumerate(p_links):
            gb0, gb1 = gene_span[g]
            r, gr, b, al = _gene_color(
                float(lfc_map.get(g, 0.0)), vmax, up_color, down_color
            )
            ribbon_col = (r, gr, b, al * 0.55)  # translucent so overlaps read
            ax.add_patch(_ribbon(sub[j], sub[j + 1], gb0, gb1, ribbon_col))

    # legend for direction
    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(facecolor=up_color, label="up-regulated"),
            Patch(facecolor=down_color, label="down-regulated"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.06),
        ncol=2,
        frameon=False,
        fontsize=8,
    )
    ax.set_title(title or f"Pathways ↔ genes  ({group})", fontsize=12, y=1.04)
    return fig
