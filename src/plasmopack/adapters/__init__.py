"""Private adapter layer — one module per external service.

Not part of the public API. Each adapter isolates a single external service
so that an upstream change touches exactly one file. See Phase 3 in
``publish/ROADMAP.md``.
"""

from __future__ import annotations

__all__: list[str] = []
