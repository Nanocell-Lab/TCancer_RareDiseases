"""Safe browser-download exports for dashboard tables."""

from __future__ import annotations

import csv
import io
import json
import zipfile
from collections.abc import Mapping, Sequence
from typing import Any


def records_to_csv(records: Sequence[Mapping[str, Any]]) -> bytes:
    """Role: Provide the records to csv operation for the export workflow."""
    if not records:
        return b""
    fields = list(dict.fromkeys(key for record in records for key in record))
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    for record in records:
        writer.writerow({key: _safe_cell(record.get(key)) for key in fields})
    return buffer.getvalue().encode("utf-8-sig")


def payload_to_json(payload: Mapping[str, Any]) -> bytes:
    """Role: Provide the payload to json operation for the export workflow."""
    return json.dumps(
        payload,
        ensure_ascii=False,
        indent=2,
        allow_nan=False,
    ).encode("utf-8")


def files_to_zip(files: Sequence[tuple[str, bytes]]) -> bytes:
    """Create a flat in-memory ZIP while rejecting unsafe archive names."""

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
        for filename, content in files:
            safe_name = filename.replace("\\", "/").rsplit("/", 1)[-1]
            if not safe_name or safe_name in {".", ".."}:
                raise ValueError("Invalid download filename")
            archive.writestr(safe_name, content)
    return buffer.getvalue()


def _safe_cell(value: Any) -> Any:
    """Role: Sanitize cell for safe serialization."""
    if value is None or isinstance(value, (int, float)):
        return value
    if isinstance(value, (list, tuple, dict)):
        value = json.dumps(value, ensure_ascii=False, allow_nan=False)
    text = str(value)
    stripped = text.lstrip(" \t\r\n")
    if stripped.startswith(("=", "+", "-", "@")) or text[:1] in {"\t", "\r", "\n"}:
        return "'" + text
    return text
