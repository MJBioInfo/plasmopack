"""A tiny on-disk cache for external-API responses.

Stdlib only: cached values are JSON files named by a SHA-256 of the cache key.
Once a response is cached, subsequent lookups are offline and reproducible —
the same principle that lets the test suite avoid the network.

The cache key includes a schema version so that changing how we store or parse
a response can invalidate old entries cleanly (bump ``SCHEMA_VERSION``).
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
_ENV_VAR = "PLASMOPACK_CACHE_DIR"


def default_cache_dir() -> Path:
    """Return the cache root, honouring ``$PLASMOPACK_CACHE_DIR`` if set."""
    override = os.environ.get(_ENV_VAR)
    if override:
        return Path(override)
    return Path.home() / ".plasmopack_cache"


def make_key(*parts: str) -> str:
    """Build a stable cache key from string parts + the schema version."""
    joined = "\x1f".join((str(SCHEMA_VERSION), *parts))
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()


class DiskCache:
    """A minimal filesystem cache of JSON-serialisable values.

    Parameters
    ----------
    root
        Directory to store cache files in. Defaults to
        :func:`default_cache_dir`. Created on first write.
    namespace
        Sub-directory grouping (e.g. an adapter name) to keep sources tidy.
    """

    def __init__(self, root: str | Path | None = None, *, namespace: str = "") -> None:
        base = Path(root) if root is not None else default_cache_dir()
        self.root = base / namespace if namespace else base

    def _path_for(self, key: str) -> Path:
        return self.root / f"{key}.json"

    def get(self, key: str) -> Any | None:
        """Return the cached value for ``key``, or ``None`` if absent."""
        path = self._path_for(key)
        if not path.exists():
            return None
        return json.loads(path.read_text(encoding="utf-8"))

    def set(self, key: str, value: Any) -> None:
        """Store ``value`` under ``key`` (must be JSON-serialisable)."""
        self.root.mkdir(parents=True, exist_ok=True)
        self._path_for(key).write_text(
            json.dumps(value, ensure_ascii=False), encoding="utf-8"
        )

    def clear(self) -> None:
        """Delete every cached file in this namespace."""
        if not self.root.exists():
            return
        for path in self.root.glob("*.json"):
            path.unlink()
