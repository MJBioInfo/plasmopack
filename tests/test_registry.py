"""Phase 1 tests — the organism registry."""

from __future__ import annotations

import pytest

from plasmopack.datasets import get_organism, list_organisms


def test_list_organisms_includes_plasmodium() -> None:
    orgs = list_organisms()
    assert "pfalciparum" in orgs
    assert "pberghei" in orgs


def test_lookup_by_key() -> None:
    org = get_organism("pberghei")
    assert org.veupathdb_site == "PlasmoDB"
    assert org.gene_id_prefix == "PBANKA_"
    assert org.taxonomy_id == 5823


def test_strain_and_species_taxa_differ() -> None:
    # Strain taxon matches the reference genome; species taxon is what UniProt
    # indexes proteins under. They are deliberately different.
    pf = get_organism("pfalciparum")
    assert pf.taxonomy_id == 36329  # P. falciparum 3D7 (strain)
    assert pf.species_taxonomy_id == 5833  # P. falciparum (species)


@pytest.mark.parametrize("alias", ["pb", "PB", "berghei", "pbanka"])
def test_lookup_by_alias(alias: str) -> None:
    assert get_organism(alias).key == "pberghei"


def test_unknown_organism_raises_helpful_error() -> None:
    with pytest.raises(KeyError, match="Available organisms"):
        get_organism("nonesuch")


def test_gene_id_matching() -> None:
    pb = get_organism("pberghei")
    assert pb.matches_gene_id("PBANKA_010010")
    assert pb.matches_gene_id("PBANKA_0100100")
    assert not pb.matches_gene_id("PF3D7_0100100")

    pf = get_organism("pfalciparum")
    assert pf.matches_gene_id("PF3D7_0100100")
    assert not pf.matches_gene_id("PBANKA_010010")


def test_cross_organism_ids_do_not_collide() -> None:
    tg = get_organism("tgondii")
    assert tg.matches_gene_id("TGME49_200010")
    assert not tg.matches_gene_id("PF3D7_0100100")
