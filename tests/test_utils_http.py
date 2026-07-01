"""Phase 3 tests — the urllib HTTP helper (retry logic, no real network).

We patch ``urllib.request.urlopen`` so no sockets are opened.
"""

from __future__ import annotations

import io
import json
import urllib.error

import pytest

from plasmopack._utils import http
from plasmopack._utils.http import HTTPError, http_get_json


class _FakeResponse(io.BytesIO):
    def __enter__(self) -> _FakeResponse:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


def _ok_body(payload: dict) -> _FakeResponse:
    return _FakeResponse(json.dumps(payload).encode("utf-8"))


def test_returns_parsed_json(monkeypatch) -> None:
    monkeypatch.setattr(
        http.urllib.request, "urlopen", lambda *a, **k: _ok_body({"hello": "world"})
    )
    assert http_get_json("https://example.org/x") == {"hello": "world"}


def test_client_error_fails_fast(monkeypatch) -> None:
    calls = {"n": 0}

    def boom(*a: object, **k: object) -> None:
        calls["n"] += 1
        raise urllib.error.HTTPError("url", 404, "Not Found", {}, None)  # type: ignore[arg-type]

    monkeypatch.setattr(http.urllib.request, "urlopen", boom)
    with pytest.raises(HTTPError) as exc:
        http_get_json("https://example.org/missing", retries=3)
    assert exc.value.status == 404
    assert calls["n"] == 1  # no retries on 4xx


def test_transient_error_retries_then_succeeds(monkeypatch) -> None:
    calls = {"n": 0}

    def flaky(*a: object, **k: object) -> _FakeResponse:
        calls["n"] += 1
        if calls["n"] < 3:
            raise urllib.error.URLError("temporary network glitch")
        return _ok_body({"ok": True})

    monkeypatch.setattr(http.urllib.request, "urlopen", flaky)
    result = http_get_json("https://example.org/x", retries=3, backoff=0)
    assert result == {"ok": True}
    assert calls["n"] == 3


def test_gives_up_after_retries(monkeypatch) -> None:
    def always_5xx(*a: object, **k: object) -> None:
        raise urllib.error.HTTPError("url", 503, "Unavailable", {}, None)  # type: ignore[arg-type]

    monkeypatch.setattr(http.urllib.request, "urlopen", always_5xx)
    with pytest.raises(HTTPError):
        http_get_json("https://example.org/x", retries=2, backoff=0)
