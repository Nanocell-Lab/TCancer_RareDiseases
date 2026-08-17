"""Small dependency-free JSON HTTP client with retries and rate limiting."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any


class ApiError(RuntimeError):
    """Raised when a remote biomedical API cannot return a valid response."""


class JsonHttpClient:
    def __init__(
        self,
        *,
        user_agent: str,
        timeout: float = 45.0,
        retries: int = 3,
        min_interval: float = 0.0,
    ) -> None:
        """Role: Compute init for the http workflow."""
        self.user_agent = user_agent
        self.timeout = timeout
        self.retries = retries
        self.min_interval = min_interval
        self._last_request_at = 0.0

    def get_json(self, url: str, params: dict[str, Any]) -> dict[str, Any]:
        """Role: Issue an HTTP GET request and return a validated JSON object."""
        query = urllib.parse.urlencode(params)
        return self._request(f"{url}?{query}")

    def post_json(self, url: str, payload: dict[str, Any]) -> dict[str, Any]:
        """Role: Issue an HTTP POST request and return a validated JSON object."""
        body = json.dumps(payload).encode("utf-8")
        return self._request(url, data=body)

    def _request(self, url: str, data: bytes | None = None) -> dict[str, Any]:
        """Role: Compute request for the http workflow."""
        safe_url = _redact_url(url)
        headers = {"Accept": "application/json", "User-Agent": self.user_agent}
        if data is not None:
            headers["Content-Type"] = "application/json"

        for attempt in range(self.retries + 1):
            self._respect_rate_limit()
            request = urllib.request.Request(url, data=data, headers=headers)
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    payload = json.load(response)
                if isinstance(payload, dict) and payload.get("error"):
                    raise ApiError(str(payload["error"]))
                if not isinstance(payload, dict):
                    raise ApiError(f"Unexpected JSON response from {safe_url}")
                return payload
            except ApiError:
                raise
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                if attempt >= self.retries:
                    detail = (
                        f"HTTP {exc.code}"
                        if isinstance(exc, urllib.error.HTTPError)
                        else type(exc).__name__
                    )
                    raise ApiError(f"Fallo consultando {safe_url}: {detail}") from exc
                time.sleep(2**attempt)

        raise ApiError(f"Fallo consultando {safe_url}")

    def _respect_rate_limit(self) -> None:
        """Role: Compute respect rate limit for the http workflow."""
        elapsed = time.monotonic() - self._last_request_at
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last_request_at = time.monotonic()


def _redact_url(url: str) -> str:
    """Role: Compute redact url for the http workflow."""
    parsed = urllib.parse.urlsplit(url)
    sensitive = {"api_key", "access_token", "token"}
    query = urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
    redacted = [
        (key, "REDACTED" if key.casefold() in sensitive else value)
        for key, value in query
    ]
    return urllib.parse.urlunsplit(
        (
            parsed.scheme,
            parsed.netloc,
            parsed.path,
            urllib.parse.urlencode(redacted),
            parsed.fragment,
        )
    )
