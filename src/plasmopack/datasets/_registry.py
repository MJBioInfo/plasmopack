"""Organism registry.

A small, bundled, version-controlled table mapping each supported organism to
its VEuPathDB site, gene-ID pattern, and NCBI taxonomy id. Every plasmopack
function accepts ``organism=`` and looks the details up here, so adding a new
organism is a one-entry change rather than a code change.

This is deliberately plain data (a dict of dataclasses) — no network, no files,
no external dependency. See Phase 1 in ``publish/ROADMAP.md``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Organism:
    """Static metadata for one supported organism.

    Attributes
    ----------
    key
        Short lowercase identifier used as the ``organism=`` argument.
    scientific_name
        Full species (and strain) name.
    common_name
        Human-friendly label.
    veupathdb_site
        The VEuPathDB project database this organism belongs to.
    site_url
        Base URL of that project database.
    gene_id_prefix
        The literal prefix real gene IDs start with (e.g. ``"PF3D7_"``).
    gene_id_regex
        Regex matching a full valid gene ID for this organism.
    taxonomy_id
        NCBI taxonomy identifier.
    """

    key: str
    scientific_name: str
    common_name: str
    veupathdb_site: str
    site_url: str
    gene_id_prefix: str
    gene_id_regex: str
    taxonomy_id: int
    aliases: tuple[str, ...] = field(default_factory=tuple)

    def matches_gene_id(self, gene_id: str) -> bool:
        """Return True if ``gene_id`` looks like an ID for this organism."""
        return re.fullmatch(self.gene_id_regex, gene_id) is not None


# ---------------------------------------------------------------------------
# The registry. Plasmodium first; other VEuPathDB organisms added as needed.
# Gene-ID regexes follow VEuPathDB conventions (verified against PlasmoDB).
# ---------------------------------------------------------------------------
_REGISTRY: dict[str, Organism] = {
    "pfalciparum": Organism(
        key="pfalciparum",
        scientific_name="Plasmodium falciparum 3D7",
        common_name="P. falciparum (3D7)",
        veupathdb_site="PlasmoDB",
        site_url="https://plasmodb.org",
        gene_id_prefix="PF3D7_",
        gene_id_regex=r"PF3D7_\d{7}",
        taxonomy_id=36329,
        aliases=("pf", "pf3d7", "falciparum"),
    ),
    "pberghei": Organism(
        key="pberghei",
        scientific_name="Plasmodium berghei ANKA",
        common_name="P. berghei (ANKA)",
        veupathdb_site="PlasmoDB",
        site_url="https://plasmodb.org",
        gene_id_prefix="PBANKA_",
        gene_id_regex=r"PBANKA_\d{6,7}",
        taxonomy_id=5823,
        aliases=("pb", "pbanka", "berghei"),
    ),
    "pvivax": Organism(
        key="pvivax",
        scientific_name="Plasmodium vivax P01",
        common_name="P. vivax (P01)",
        veupathdb_site="PlasmoDB",
        site_url="https://plasmodb.org",
        gene_id_prefix="PVP01_",
        gene_id_regex=r"PVP01_\d{7}",
        taxonomy_id=126793,
        aliases=("pv", "vivax"),
    ),
    "pknowlesi": Organism(
        key="pknowlesi",
        scientific_name="Plasmodium knowlesi H",
        common_name="P. knowlesi (H)",
        veupathdb_site="PlasmoDB",
        site_url="https://plasmodb.org",
        gene_id_prefix="PKNH_",
        gene_id_regex=r"PKNH_\d{7}",
        taxonomy_id=5851,
        aliases=("pk", "knowlesi"),
    ),
    "tgondii": Organism(
        key="tgondii",
        scientific_name="Toxoplasma gondii ME49",
        common_name="T. gondii (ME49)",
        veupathdb_site="ToxoDB",
        site_url="https://toxodb.org",
        gene_id_prefix="TGME49_",
        gene_id_regex=r"TGME49_\d{6,9}",
        taxonomy_id=508771,
        aliases=("tg", "tgme49", "gondii", "toxo"),
    ),
}

# Build an alias lookup once, at import time.
_ALIAS_INDEX: dict[str, str] = {}
for _org in _REGISTRY.values():
    _ALIAS_INDEX[_org.key.lower()] = _org.key
    for _alias in _org.aliases:
        _ALIAS_INDEX[_alias.lower()] = _org.key


def get_organism(organism: str) -> Organism:
    """Look up an :class:`Organism` by key or alias (case-insensitive).

    Parameters
    ----------
    organism
        Registry key (e.g. ``"pberghei"``) or a known alias (e.g. ``"pb"``).

    Raises
    ------
    KeyError
        If the organism is not in the registry, with a helpful message
        listing what is available.
    """
    resolved = _ALIAS_INDEX.get(organism.lower())
    if resolved is None:
        available = ", ".join(sorted(_REGISTRY))
        raise KeyError(
            f"Unknown organism {organism!r}. "
            f"Available organisms: {available}. "
            f"To add one, extend plasmopack.datasets._registry."
        )
    return _REGISTRY[resolved]


def list_organisms() -> list[str]:
    """Return the sorted list of supported organism keys."""
    return sorted(_REGISTRY)
