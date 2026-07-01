"""Shared pytest fixtures."""

from __future__ import annotations

import pytest
from anndata import AnnData

from tests.synthetic import make_dummy_pberghei


@pytest.fixture
def dummy_pberghei() -> AnnData:
    """A small synthetic P. berghei AnnData with planted timepoint markers."""
    return make_dummy_pberghei(n_cells=300, n_genes=800, seed=0)
