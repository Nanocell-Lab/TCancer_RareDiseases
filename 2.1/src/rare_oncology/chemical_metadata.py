"""Chemical metadata enrichment for therapies identified in ChEMBL."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Mapping, Protocol

from .http import JsonHttpClient


_CHEMBL_ID = re.compile(r"^CHEMBL\d+$", re.IGNORECASE)


@dataclass(frozen=True)
class ChemicalMetadata:
    """Structure metadata returned for one ChEMBL molecule."""

    chembl_id: str
    molecule_type: str | None
    canonical_smiles: str | None


class ChemicalMetadataResolver(Protocol):
    """Resolve stable molecule identifiers without interpreting therapy evidence."""

    def cache_identity(self) -> Mapping[str, Any]:
        """Return configuration that participates in therapy cache keys."""
        ...

    def resolve(self, chembl_ids: tuple[str, ...]) -> dict[str, ChemicalMetadata]:
        """Return available metadata keyed by normalized ChEMBL identifier."""
        ...


class ChemicalGetJsonClient(Protocol):
    """Minimal HTTP surface required by the ChEMBL resolver."""

    def get_json(self, url: str, params: dict[str, Any]) -> dict[str, Any]:
        """Retrieve a JSON payload from the molecule endpoint."""
        ...


class ChEMBLChemicalMetadataClient:
    """Dependency-free client for the official ChEMBL molecule endpoint."""

    url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

    def __init__(
        self,
        *,
        http: ChemicalGetJsonClient | None = None,
        batch_size: int = 50,
    ) -> None:
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        self.http = http or JsonHttpClient(
            user_agent="rare-oncology/2.1 (ChEMBL structure lookup)",
            retries=2,
            min_interval=0.1,
        )
        self.batch_size = batch_size
        self._cache: dict[str, ChemicalMetadata | None] = {}

    def cache_identity(self) -> dict[str, Any]:
        return {
            "url": self.url,
            "batch_size": self.batch_size,
            "fields": "molecule_type,canonical_smiles",
        }

    def resolve(self, chembl_ids: tuple[str, ...]) -> dict[str, ChemicalMetadata]:
        normalized = tuple(
            sorted(
                {
                    value.strip().upper()
                    for value in chembl_ids
                    if _CHEMBL_ID.fullmatch(value.strip())
                }
            )
        )
        missing = [value for value in normalized if value not in self._cache]
        for start in range(0, len(missing), self.batch_size):
            batch = missing[start : start + self.batch_size]
            payload = self.http.get_json(
                self.url,
                {
                    "molecule_chembl_id__in": ",".join(batch),
                    "limit": len(batch),
                },
            )
            found: set[str] = set()
            for item in payload.get("molecules") or []:
                if not isinstance(item, dict):
                    continue
                chembl_id = str(item.get("molecule_chembl_id") or "").upper()
                if chembl_id not in batch:
                    continue
                structures = item.get("molecule_structures") or {}
                if not isinstance(structures, dict):
                    structures = {}
                metadata = ChemicalMetadata(
                    chembl_id=chembl_id,
                    molecule_type=_optional_text(item.get("molecule_type")),
                    canonical_smiles=_optional_text(
                        structures.get("canonical_smiles")
                    ),
                )
                self._cache[chembl_id] = metadata
                found.add(chembl_id)
            for chembl_id in set(batch) - found:
                self._cache[chembl_id] = None
        return {
            chembl_id: metadata
            for chembl_id in normalized
            if (metadata := self._cache.get(chembl_id)) is not None
        }


def _optional_text(value: object) -> str | None:
    text = str(value).strip() if value is not None else ""
    return text or None
