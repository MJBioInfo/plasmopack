"""Tools — analytical methods.

DEG (Phase 5), compositional analysis (Phase 6), enrichment (Phase 7),
transcription-factor analysis (Phase 8), orthologs (Phase 9), and the
cell-cell interaction submodule ``tl.cci`` (Phase 11).
"""

from __future__ import annotations

from plasmopack.tl import cci
from plasmopack.tl._deg import deg, list_deg_methods
from plasmopack.tl._enrich import enrich

__all__ = ["cci", "deg", "enrich", "list_deg_methods"]
