"""Clients and parsers for ClinVar/Orphanet and GDC/TCGA."""

from __future__ import annotations

import os
from collections.abc import Iterable
from typing import Any

from .http import ApiError, JsonHttpClient
from .models import CancerObservation, RareVariant
from .normalization import normalize_protein_change, primary_protein_change


CLINVAR_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GDC_SSM_OCCURRENCES = "https://api.gdc.cancer.gov/ssm_occurrences"
GDC_STATUS = "https://api.gdc.cancer.gov/status"


class NcbiClinVarClient:
    """Retrieve gene variants whose ClinVar germline traits reference Orphanet."""

    def __init__(
        self,
        *,
        email: str | None = None,
        api_key: str | None = None,
        http: JsonHttpClient | None = None,
        page_size: int = 500,
    ) -> None:
        """Role: Compute init for the clients workflow."""
        self.email = email or os.getenv("NCBI_EMAIL")
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.page_size = page_size
        interval = 0.11 if self.api_key else 0.34
        self.http = http or JsonHttpClient(
            user_agent="rare-oncology/2.0",
            min_interval=interval,
        )

    def fetch_rare_variants(
        self, gene: str, *, max_records: int | None = None
    ) -> list[RareVariant]:
        """Role: Retrieve and normalize rare variants data."""
        gene = _validate_gene(gene)
        search = self.http.get_json(
            CLINVAR_ESEARCH,
            self._ncbi_params(
                db="clinvar",
                term=f"{gene}[gene] AND single_gene[prop]",
                usehistory="y",
                retmode="json",
                retmax="0",
            ),
        )
        result = search.get("esearchresult", {})
        count = int(result.get("count", 0))
        if max_records is not None:
            if max_records < 1:
                raise ValueError("max_records debe ser positivo")
            count = min(count, max_records)
        if count == 0:
            return []

        query_key = result.get("querykey")
        webenv = result.get("webenv")
        if not query_key or not webenv:
            raise ApiError("ClinVar no devolvio querykey/webenv para paginar")

        variants: dict[tuple[str, str], RareVariant] = {}
        for start in range(0, count, self.page_size):
            batch_size = min(self.page_size, count - start)
            payload = self.http.get_json(
                CLINVAR_ESUMMARY,
                self._ncbi_params(
                    db="clinvar",
                    query_key=query_key,
                    WebEnv=webenv,
                    retstart=str(start),
                    retmax=str(batch_size),
                    retmode="json",
                ),
            )
            for summary in _summary_documents(payload):
                parsed = _parse_clinvar_summary(gene, summary)
                if parsed is not None:
                    variants[(parsed.clinvar_id, parsed.protein_change)] = parsed
        return sorted(variants.values(), key=lambda item: (item.protein_change, item.clinvar_id))

    def _ncbi_params(self, **params: str) -> dict[str, str]:
        """Role: Compute ncbi params for the clients workflow."""
        params["tool"] = "rare_oncology"
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        return params


class GdcClient:
    """Retrieve canonical somatic protein changes from GDC occurrence records."""

    def __init__(
        self,
        *,
        http: JsonHttpClient | None = None,
        page_size: int = 1000,
    ) -> None:
        """Role: Compute init for the clients workflow."""
        self.page_size = page_size
        self.http = http or JsonHttpClient(user_agent="rare-oncology/2.0")
        self.release_metadata: dict[str, Any] = {}

    def fetch_cancer_observations(
        self,
        gene: str,
        *,
        program: str | None = "TCGA",
        projects: Iterable[str] | None = None,
    ) -> list[CancerObservation]:
        """Role: Retrieve and normalize cancer observations data."""
        gene = _validate_gene(gene)
        try:
            self.release_metadata = self.http.get_json(GDC_STATUS, {})
        except ApiError:
            self.release_metadata = {}
        selected_projects = tuple(sorted(set(projects or ())))
        filters: list[dict[str, Any]] = [
            _in_filter("ssm.consequence.transcript.gene.symbol", [gene])
        ]
        if program:
            filters.append(_in_filter("case.project.program.name", [program]))
        if selected_projects:
            filters.append(_in_filter("case.project.project_id", list(selected_projects)))

        fields = ",".join(
            [
                "ssm.ssm_id",
                "ssm.chromosome",
                "ssm.genomic_dna_change",
                "ssm.consequence.transcript.annotation.hgvsp_short",
                "ssm.consequence.transcript.annotation.hgvsp",
                "ssm.consequence.transcript.gene.symbol",
                "ssm.consequence.transcript.is_canonical",
                "ssm.consequence.transcript.transcript_id",
                "case.case_id",
                "case.submitter_id",
                "case.project.project_id",
            ]
        )

        observations: dict[tuple[str, str, str], CancerObservation] = {}
        start = 0
        while True:
            payload = self.http.post_json(
                GDC_SSM_OCCURRENCES,
                {
                    "filters": {"op": "and", "content": filters},
                    "fields": fields,
                    "format": "JSON",
                    "size": self.page_size,
                    "from": start,
                },
            )
            data = payload.get("data", {})
            hits = data.get("hits", [])
            for hit in hits:
                parsed = _parse_gdc_occurrence(gene, hit)
                if parsed is not None:
                    key = (parsed.ssm_id, parsed.case_id, parsed.protein_change)
                    observations[key] = parsed

            pagination = data.get("pagination", {})
            total = int(pagination.get("total", len(hits)))
            start += len(hits)
            if not hits or start >= total:
                break

        return sorted(
            observations.values(),
            key=lambda item: (item.protein_change, item.project_id, item.case_id),
        )


def _summary_documents(payload: dict[str, Any]) -> Iterable[dict[str, Any]]:
    """Role: Compute summary documents for the clients workflow."""
    result = payload.get("result", {})
    for uid in result.get("uids", []):
        document = result.get(str(uid))
        if isinstance(document, dict):
            yield document


def _parse_clinvar_summary(gene: str, summary: dict[str, Any]) -> RareVariant | None:
    """Role: Parse and validate clinvar summary from its source representation."""
    gene_symbols = {str(item.get("symbol", "")).upper() for item in summary.get("genes", [])}
    if gene not in gene_symbols:
        return None

    classification = summary.get("germline_classification") or {}
    diseases: set[str] = set()
    orphanet_codes: set[str] = set()
    condition_orphanet_pairs: set[tuple[str, str]] = set()
    for trait in classification.get("trait_set", []):
        codes = {
            str(xref.get("db_id"))
            for xref in trait.get("trait_xrefs", [])
            if str(xref.get("db_source", "")).lower() == "orphanet" and xref.get("db_id")
        }
        if codes:
            orphanet_codes.update(codes)
            disease = str(trait.get("trait_name") or "Condicion ClinVar sin nombre")
            diseases.add(disease)
            condition_orphanet_pairs.update((disease, code) for code in codes)
    if not orphanet_codes:
        return None

    protein = primary_protein_change(summary)
    if protein is None:
        return None
    normalized, raw = protein

    canonical_spdi = None
    for variation in summary.get("variation_set", []):
        if variation.get("canonical_spdi"):
            canonical_spdi = str(variation["canonical_spdi"])
            break

    return RareVariant(
        gene=gene,
        protein_change=normalized,
        raw_protein_change=raw,
        clinvar_id=str(summary.get("uid", "")),
        clinvar_accession=str(
            summary.get("accession_version") or summary.get("accession") or ""
        ),
        diseases=tuple(sorted(diseases)),
        orphanet_codes=tuple(sorted(orphanet_codes)),
        condition_orphanet_pairs=tuple(sorted(condition_orphanet_pairs)),
        germline_classification=str(classification.get("description", "")),
        molecular_consequences=tuple(
            sorted(str(value) for value in summary.get("molecular_consequence_list", []))
        ),
        canonical_spdi=canonical_spdi,
    )


def _parse_gdc_occurrence(gene: str, hit: dict[str, Any]) -> CancerObservation | None:
    """Role: Parse and validate gdc occurrence from its source representation."""
    ssm = hit.get("ssm") or {}
    candidates: list[dict[str, Any]] = []
    for consequence in ssm.get("consequence", []):
        transcript = consequence.get("transcript") or {}
        symbol = str((transcript.get("gene") or {}).get("symbol", "")).upper()
        annotation = transcript.get("annotation") or {}
        if symbol == gene and annotation.get("hgvsp_short"):
            candidates.append(transcript)
    if not candidates:
        return None

    canonical = [item for item in candidates if item.get("is_canonical") is True]
    transcript = canonical[0] if canonical else candidates[0]
    annotation = transcript.get("annotation") or {}
    raw_change = str(annotation.get("hgvsp_short", ""))
    protein_change = normalize_protein_change(raw_change)
    if protein_change is None or protein_change.endswith("="):
        return None

    case = hit.get("case") or {}
    return CancerObservation(
        gene=gene,
        protein_change=protein_change,
        raw_protein_change=raw_change,
        ssm_id=str(ssm.get("ssm_id", "")),
        case_id=str(case.get("case_id", "")),
        patient_id=str(case.get("submitter_id", "")),
        project_id=str((case.get("project") or {}).get("project_id", "")),
        transcript_id=(
            str(transcript["transcript_id"])
            if transcript.get("transcript_id")
            else None
        ),
        genomic_dna_change=(
            str(ssm["genomic_dna_change"]) if ssm.get("genomic_dna_change") else None
        ),
    )


def _validate_gene(gene: str) -> str:
    """Role: Validate gene constraints before execution."""
    normalized = gene.strip().upper()
    if not normalized or not normalized.replace("-", "").isalnum():
        raise ValueError(f"Invalid gene symbol: {gene!r}")
    return normalized


def _in_filter(field: str, values: list[str]) -> dict[str, Any]:
    """Role: Compute in filter for the clients workflow."""
    return {"op": "in", "content": {"field": field, "value": values}}
