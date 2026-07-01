"""Minimal HTTP helper built on the standard library only.

We deliberately avoid ``requests``/``httpx`` to honour the minimum-dependency
contract: ``urllib`` ships with Python and never breaks because a third-party
HTTP library changed. This module provides just what the adapter layer needs —
a JSON GET with sensible retries and a package-identifying User-Agent.

Retry policy: transient failures (network errors, HTTP 5xx, HTTP 429) are
retried with exponential backoff; client errors (other HTTP 4xx) fail fast.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from typing import Any

DEFAULT_TIMEOUT = 30.0
DEFAULT_RETRIES = 3
_RETRYABLE_STATUS = {429, 500, 502, 503, 504}


class HTTPError(RuntimeError):
    """Raised when an HTTP request ultimately fails."""

    def __init__(self, message: str, *, status: int | None = None) -> None:
        super().__init__(message)
        self.status = status


def _user_agent() -> str:
    from plasmopack import __version__

    return f"plasmopack/{__version__} (https://github.com/MJBioInfo/plasmopack)"


def http_get_json(
    url: str,
    *,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    headers: dict[str, str] | None = None,
    backoff: float = 0.5,
) -> Any:
    """GET ``url`` and parse the response body as JSON.

    Parameters
    ----------
    url
        The full request URL (including any query string).
    timeout
        Per-attempt socket timeout in seconds.
    retries
        Total number of attempts for transient failures.
    headers
        Extra request headers. ``Accept: application/json`` and a
        ``User-Agent`` are added automatically.
    backoff
        Base seconds for exponential backoff between retries.

    Returns
    -------
    Any
        The parsed JSON (dict or list).

    Raises
    ------
    HTTPError
        On a non-retryable HTTP status, or after exhausting retries.
    """
    req_headers = {
        "Accept": "application/json",
        "User-Agent": _user_agent(),
    }
    if headers:
        req_headers.update(headers)

    request = urllib.request.Request(url, headers=req_headers, method="GET")

    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as resp:
                raw = resp.read().decode("utf-8")
            return json.loads(raw)
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code not in _RETRYABLE_STATUS:
                raise HTTPError(
                    f"GET {url} failed with HTTP {exc.code}", status=exc.code
                ) from exc
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
        # transient: back off unless this was the last attempt
        if attempt < retries - 1:
            time.sleep(backoff * (2**attempt))

    raise HTTPError(f"GET {url} failed after {retries} attempts: {last_error}")
