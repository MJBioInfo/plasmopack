"""Phase 2 tests — GMT read/write and the provenance sidecar."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from plasmopack.io import GeneSet, GeneSets, read_gmt, write_gmt


def _example() -> GeneSets:
    return GeneSets(
        sets=[
            GeneSet("GO:0006412", "translation", ["g1", "g2", "g3"], "BP"),
            GeneSet("GO:0003735", "ribosome", ["g1", "g2"], "MF"),
        ],
        metadata={"source": "test", "version": "PlasmoDB-68", "organism": "pberghei"},
    )


def test_write_then_read_roundtrip(tmp_path: Path) -> None:
    gs = _example()
    path = write_gmt(gs, tmp_path / "sets.gmt")
    assert path.exists()

    back = read_gmt(path)
    assert back.ids == ["GO:0006412", "GO:0003735"]
    assert back["GO:0006412"].genes == ["g1", "g2", "g3"]
    assert back["GO:0003735"].description == "MF"


def test_sidecar_written_and_loaded(tmp_path: Path) -> None:
    gs = _example()
    path = write_gmt(gs, tmp_path / "sets.gmt")

    sidecar = tmp_path / "sets.gmt.meta.json"
    assert sidecar.exists()
    meta = json.loads(sidecar.read_text())
    assert meta["version"] == "PlasmoDB-68"
    assert meta["_summary"]["n_sets"] == 2
    assert meta["_summary"]["n_genes_total"] == 3  # g1,g2,g3 unique

    # read_gmt should pull the sidecar back into metadata
    back = read_gmt(path)
    assert back.metadata["version"] == "PlasmoDB-68"


def test_read_without_sidecar_records_source(tmp_path: Path) -> None:
    gs = _example()
    path = write_gmt(gs, tmp_path / "sets.gmt", write_sidecar=False)
    assert not (tmp_path / "sets.gmt.meta.json").exists()

    back = read_gmt(path)
    assert back.metadata["source"] == str(path)


def test_gmt_is_standard_tab_format(tmp_path: Path) -> None:
    path = write_gmt(_example(), tmp_path / "sets.gmt")
    first = path.read_text().splitlines()[0].split("\t")
    assert first[0] == "GO:0006412"  # set id
    assert first[1] == "BP"  # description
    assert first[2:] == ["g1", "g2", "g3"]  # genes


def test_malformed_line_raises(tmp_path: Path) -> None:
    bad = tmp_path / "bad.gmt"
    bad.write_text("only_one_field\n")
    with pytest.raises(ValueError, match="at least 2"):
        read_gmt(bad)
