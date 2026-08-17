"""Complementary therapy lookup with three deliberately separate evidence layers."""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Protocol, Sequence

from .http import ApiError, JsonHttpClient
from .models import (
    AnalysisResult,
    SourceStatus,
    TherapyAssociation,
    TherapyResult,
)
from .normalization import normalize_protein_change


VARIANT_THERAPY_EVIDENCE = "variant_therapy_evidence"
REGULATORY_THERAPY = "regulatory_therapy"
TRIAL_AND_TARGET_DISCOVERY = "trial_and_target_discovery"


@dataclass(frozen=True)
class TherapyContext:
    gene: str
    protein_changes: tuple[str, ...]
    disease_ids: tuple[str, ...]
    disease_names: tuple[str, ...]

    @classmethod
    def from_analysis(cls, result: AnalysisResult) -> "TherapyContext":
        """Role: Build a therapy-query context from a completed gene analysis."""
        return cls(
            gene=result.gene.upper(),
            protein_changes=tuple(
                sorted({item.protein_change for item in result.shared_mutations})
            ),
            disease_ids=tuple(
                sorted(
                    {
                        f"ORPHA:{code}"
                        for item in result.shared_mutations
                        for code in item.orphanet_codes
                    }
                )
            ),
            disease_names=tuple(
                sorted(
                    {
                        disease
                        for item in result.shared_mutations
                        for disease in item.diseases
                    }
                )
            ),
        )


class TherapyProvider(Protocol):
    @property
    def name(self) -> str:
        """Role: Expose the stable therapy-provider name."""
        ...

    @property
    def layer(self) -> str:
        """Role: Expose the independent therapy evidence layer."""
        ...

    @property
    def source_version(self) -> str | None:
        """Role: Expose the observed source release or version."""
        ...

    @property
    def url(self) -> str | None:
        """Role: Expose the configured source endpoint."""
        ...

    def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
        """Role: Collect therapy associations for one gene context."""
        ...


class PostJsonClient(Protocol):
    def post_json(self, url: str, payload: dict[str, Any]) -> dict[str, Any]:
        """Role: Submit a JSON query to a therapy endpoint."""
        ...


class GetJsonClient(Protocol):
    def get_json(self, url: str, params: dict[str, Any]) -> dict[str, Any]:
        """Role: Retrieve JSON from a public therapy endpoint."""
        ...


class OncoKBGetJsonClient(Protocol):
    def get_json(
        self,
        url: str,
        params: dict[str, Any],
        *,
        headers: dict[str, str],
    ) -> dict[str, Any]:
        """Role: Retrieve authorized JSON from the OncoKB endpoint."""
        ...


class CivicTherapyProvider:
    """Accepted predictive CIViC evidence for the matched protein changes."""

    name = "civic"
    layer = VARIANT_THERAPY_EVIDENCE
    source_version: str | None = None
    url = "https://civicdb.org/api/graphql"
    license_name = "CC0-1.0"

    _VARIANTS_QUERY = """
      query($gene: String!, $after: String) {
        gene(entrezSymbol: $gene) {
          variants(first: 100, after: $after) {
            pageInfo { hasNextPage endCursor }
            nodes { id name }
          }
        }
      }
    """
    _EVIDENCE_QUERY = """
      query($variantId: Int!, $after: String) {
        evidenceItems(
          variantId: $variantId,
          status: ACCEPTED,
          evidenceType: PREDICTIVE,
          first: 100,
          after: $after
        ) {
          pageInfo { hasNextPage endCursor }
          nodes {
            id status evidenceType evidenceLevel evidenceDirection significance
            evidenceRating therapyInteractionType description link
            molecularProfile { id name isMultiVariant }
            disease { id name doid link }
            therapies { id name ncitId link }
            source { citationId sourceType link }
          }
        }
      }
    """

    def __init__(
        self,
        *,
        http: PostJsonClient | None = None,
        max_variants: int = 25,
    ) -> None:
        """Role: Compute init for the therapy workflow."""
        if max_variants < 1:
            raise ValueError("max_variants must be positive")
        self.http = http or JsonHttpClient(
            user_agent="rare-oncology/2.1 (CIViC therapy lookup)",
            retries=2,
            min_interval=0.34,
        )
        self.max_variants = max_variants

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {"url": self.url, "max_variants": self.max_variants}

    def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
        """Role: Collect source-specific records for the current analysis context."""
        requested = set(context.protein_changes)
        if not requested:
            return []
        matching_variants: list[tuple[int, str]] = []
        cursor: str | None = None
        seen_cursors: set[str] = set()
        while True:
            payload = self._graphql(
                self._VARIANTS_QUERY,
                {"gene": context.gene, "after": cursor},
            )
            gene = payload.get("data", {}).get("gene")
            if not gene:
                return []
            connection = gene.get("variants") or {}
            for node in connection.get("nodes") or []:
                normalized = normalize_protein_change(str(node.get("name") or ""))
                if normalized in requested:
                    matching_variants.append((int(node["id"]), normalized))
            page = connection.get("pageInfo") or {}
            cursor = _next_graphql_cursor(
                page, cursor, seen_cursors, "CIViC variants"
            )
            if cursor is None:
                break

        retrieved_at = _utc_now()
        records: list[TherapyAssociation] = []
        for variant_id, protein_change in matching_variants[: self.max_variants]:
            cursor = None
            seen_cursors = set()
            while True:
                payload = self._graphql(
                    self._EVIDENCE_QUERY,
                    {"variantId": variant_id, "after": cursor},
                )
                connection = payload.get("data", {}).get("evidenceItems") or {}
                for node in connection.get("nodes") or []:
                    record = self._record(node, context, protein_change, retrieved_at)
                    if record is not None:
                        records.append(record)
                page = connection.get("pageInfo") or {}
                cursor = _next_graphql_cursor(
                    page,
                    cursor,
                    seen_cursors,
                    f"CIViC evidence variant {variant_id}",
                )
                if cursor is None:
                    break
        return _deduplicate(records)

    def _graphql(self, query: str, variables: dict[str, Any]) -> dict[str, Any]:
        """Role: Compute graphql for the therapy workflow."""
        payload = self.http.post_json(
            self.url, {"query": query, "variables": variables}
        )
        errors = payload.get("errors")
        if errors:
            message = "; ".join(
                str(item.get("message", item)) if isinstance(item, dict) else str(item)
                for item in errors
            )
            raise ApiError(f"CIViC GraphQL: {message}")
        return payload

    def _record(
        self,
        node: dict[str, Any],
        context: TherapyContext,
        protein_change: str,
        retrieved_at: str,
    ) -> TherapyAssociation | None:
        """Role: Compute record for the therapy workflow."""
        therapies = [item for item in node.get("therapies") or [] if item.get("name")]
        if not therapies:
            return None
        names = tuple(str(item["name"]) for item in therapies)
        therapy_ids = tuple(
            str(item.get("ncitId") or f"CIViC:{item.get('id')}") for item in therapies
        )
        disease = node.get("disease") or {}
        profile = node.get("molecularProfile") or {}
        source = node.get("source") or {}
        evidence_id = str(node.get("id"))
        interaction = _optional_text(node.get("therapyInteractionType"))
        return TherapyAssociation(
            layer=self.layer,
            source=self.name,
            source_record_id=evidence_id,
            gene=context.gene,
            association_scope=(
                "molecular_profile"
                if profile.get("isMultiVariant")
                else "candidate_protein_variant"
            ),
            therapy_name=" + ".join(names),
            therapy_components=names,
            therapy_status="curated_predictive_evidence",
            therapy_id=" + ".join(therapy_ids),
            protein_change=protein_change,
            molecular_profile=_optional_text(profile.get("name")),
            disease_ids=(
                (f"DOID:{disease['doid']}",) if disease.get("doid") else ()
            ),
            disease_names=((str(disease["name"]),) if disease.get("name") else ()),
            effect=_civic_effect(node.get("significance")),
            evidence_type=_snake_case(node.get("evidenceType")),
            evidence_level=_optional_text(node.get("evidenceLevel")),
            evidence_direction=_snake_case(node.get("evidenceDirection")),
            evidence_rating=_optional_float(node.get("evidenceRating")),
            citation_id=_optional_text(source.get("citationId")),
            description=_optional_text(node.get("description")),
            source_url=_absolute_civic_url(node.get("link")),
            retrieved_at=retrieved_at,
            source_version=self.source_version,
            license=self.license_name,
            details=(("therapy_interaction_type", interaction),) if interaction else (),
        )


class OncoKBTherapyProvider:
    """Licensed OncoKB protein-variant treatment implications."""

    name = "oncokb"
    layer = VARIANT_THERAPY_EVIDENCE
    url = "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange"
    license_name = "OncoKB licensed API; use subject to account terms"

    def __init__(
        self,
        token: str,
        *,
        http: OncoKBGetJsonClient | None = None,
        max_variants: int = 25,
        reference_genome: str = "GRCh38",
    ) -> None:
        """Role: Compute init for the therapy workflow."""
        clean_token = token.strip()
        if not clean_token:
            raise ValueError("OncoKB requires a non-empty token")
        if max_variants < 1:
            raise ValueError("max_variants must be positive")
        if reference_genome not in {"GRCh37", "GRCh38"}:
            raise ValueError("reference_genome must be GRCh37 or GRCh38")
        self._token = clean_token
        self.http = http or JsonHttpClient(
            user_agent="rare-oncology/2.1 (OncoKB licensed lookup)",
            retries=2,
            min_interval=0.1,
        )
        self.max_variants = max_variants
        self.reference_genome = reference_genome
        self.source_version: str | None = None
        self.last_status = "ok"
        self.last_message: str | None = None

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {
            "url": self.url,
            "max_variants": self.max_variants,
            "reference_genome": self.reference_genome,
        }

    def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
        """Role: Collect source-specific records for the current analysis context."""
        self.last_status = "ok"
        self.last_message = None
        changes = context.protein_changes[: self.max_variants]
        if len(context.protein_changes) > len(changes):
            self.last_status = "partial"
            self.last_message = (
                f"variants queried {len(changes)}/{len(context.protein_changes)}"
            )
        versions: set[str] = set()
        retrieved_at = _utc_now()
        records: list[TherapyAssociation] = []
        for protein_change in changes:
            payload = self.http.get_json(
                self.url,
                {
                    "referenceGenome": self.reference_genome,
                    "hugoSymbol": context.gene,
                    "alteration": protein_change,
                },
                headers={"Authorization": f"Bearer {self._token}"},
            )
            version = _optional_text(payload.get("dataVersion"))
            if version:
                versions.add(version)
            records.extend(
                self._treatment_records(
                    payload, context.gene, protein_change, retrieved_at, version
                )
            )
        self.source_version = ", ".join(sorted(versions)) or None
        return _deduplicate(records)

    def _treatment_records(
        self,
        payload: dict[str, Any],
        gene: str,
        protein_change: str,
        retrieved_at: str,
        source_version: str | None,
    ) -> list[TherapyAssociation]:
        """Role: Compute treatment records for the therapy workflow."""
        records: list[TherapyAssociation] = []
        for treatment in payload.get("treatments") or []:
            if not isinstance(treatment, dict):
                continue
            drugs = [
                item
                for item in treatment.get("drugs") or []
                if isinstance(item, dict) and _optional_text(item.get("drugName"))
            ]
            components = tuple(
                str(item["drugName"]).strip() for item in drugs
            )
            if not components:
                continue
            drug_ids = tuple(
                str(item["ncitCode"]).strip()
                for item in drugs
                if _optional_text(item.get("ncitCode"))
            )
            level = _optional_text(treatment.get("level"))
            cancer = treatment.get("levelAssociatedCancerType") or {}
            if not isinstance(cancer, dict):
                cancer = {}
            cancer_code = _optional_text(cancer.get("code"))
            cancer_id = _optional_text(cancer.get("id"))
            main_type = cancer.get("mainType") or {}
            if not isinstance(main_type, dict):
                main_type = {}
            cancer_name = _optional_text(cancer.get("name")) or _optional_text(
                main_type.get("name")
            )
            disease_id = (
                f"ONCOTREE:{cancer_code}"
                if cancer_code
                else (f"ONCOTREE_ID:{cancer_id}" if cancer_id else None)
            )
            pmids = tuple(
                str(value).strip()
                for value in treatment.get("pmids") or []
                if str(value).strip()
            )
            effect = (
                "resistance"
                if level and level.upper().startswith("LEVEL_R")
                else "sensitivity_response"
            )
            identity = "|".join(
                (
                    gene,
                    protein_change,
                    level or "unlevelled",
                    disease_id or cancer_name or "unspecified_cancer",
                    "+".join(drug_ids or components),
                )
            )
            details = {
                "fda_level": _optional_text(treatment.get("fdaLevel")) or "",
                "oncogenic": _optional_text(payload.get("oncogenic")) or "",
                "highest_sensitive_level": (
                    _optional_text(payload.get("highestSensitiveLevel")) or ""
                ),
                "highest_resistance_level": (
                    _optional_text(payload.get("highestResistanceLevel")) or ""
                ),
                "pmids": "|".join(pmids),
                "last_update": _optional_text(payload.get("lastUpdate")) or "",
            }
            records.append(
                TherapyAssociation(
                    layer=self.layer,
                    source=self.name,
                    source_record_id=identity,
                    gene=gene,
                    association_scope="oncokb_protein_variant",
                    therapy_name=" + ".join(components),
                    therapy_components=components,
                    therapy_status="oncokb_curated_evidence",
                    retrieved_at=retrieved_at,
                    therapy_id=" + ".join(f"NCIT:{value}" for value in drug_ids)
                    or None,
                    protein_change=protein_change,
                    molecular_profile=f"{gene} {protein_change}",
                    disease_ids=(disease_id,) if disease_id else (),
                    disease_names=(cancer_name,) if cancer_name else (),
                    effect=effect,
                    evidence_type="predictive",
                    evidence_level=level,
                    evidence_direction="supports",
                    citation_id=(
                        "|".join(f"PMID:{value}" for value in pmids) or None
                    ),
                    source_url=f"https://www.oncokb.org/gene/{gene}/{protein_change}",
                    source_version=source_version,
                    license=self.license_name,
                    details=tuple(sorted(details.items())),
                )
            )
        return records


class ClinicalTrialsTherapyProvider:
    """Gene/variant mentions in ClinicalTrials.gov, never interpreted as efficacy."""

    name = "clinicaltrials"
    layer = TRIAL_AND_TARGET_DISCOVERY
    source_version = "API v2"
    url = "https://clinicaltrials.gov/api/v2/studies"
    license_name = "Public domain / ClinicalTrials.gov terms"

    def __init__(
        self,
        *,
        http: GetJsonClient | None = None,
        max_studies: int = 100,
    ) -> None:
        """Role: Compute init for the therapy workflow."""
        if max_studies < 1:
            raise ValueError("max_studies must be positive")
        self.http = http or JsonHttpClient(
            user_agent="rare-oncology/2.1 (ClinicalTrials.gov therapy lookup)",
            retries=2,
        )
        self.max_studies = max_studies

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {"url": self.url, "max_studies": self.max_studies}

    def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
        """Role: Collect source-specific records for the current analysis context."""
        self.last_status = "ok"
        self.last_message: str | None = None
        warnings: list[str] = []
        changes = context.protein_changes[:10]
        if len(context.protein_changes) > len(changes):
            warnings.append(
                f"variants queried 10/{len(context.protein_changes)}"
            )
        query_terms = tuple(changes) + (f"{context.gene} gene",)
        studies_by_id: dict[str, dict[str, Any]] = {}
        for query_term in dict.fromkeys(query_terms):
            token: str | None = None
            seen_tokens: set[str] = set()
            scanned = 0
            while scanned < self.max_studies:
                params: dict[str, Any] = {
                    "query.term": query_term,
                    "pageSize": min(1000, self.max_studies - scanned),
                    "format": "json",
                }
                if token:
                    params["pageToken"] = token
                payload = self.http.get_json(self.url, params)
                page_studies = payload.get("studies") or []
                scanned += len(page_studies)
                for study in page_studies:
                    nct_id = _study_nct_id(study)
                    if nct_id:
                        studies_by_id.setdefault(nct_id, study)
                next_token = _optional_text(payload.get("nextPageToken"))
                if not next_token:
                    break
                if not page_studies:
                    raise ApiError(
                        "ClinicalTrials.gov returned an empty page with a continuation token"
                    )
                if next_token == token or next_token in seen_tokens:
                    raise ApiError(
                        "ClinicalTrials.gov repeated the pagination token"
                    )
                if scanned >= self.max_studies:
                    warnings.append(
                        f"query '{query_term}' truncated after {self.max_studies} studies"
                    )
                    break
                seen_tokens.add(next_token)
                token = next_token

        retrieved_at = _utc_now()
        records: list[TherapyAssociation] = []
        for study in studies_by_id.values():
            records.extend(self._study_records(study, context, retrieved_at))
        if warnings:
            self.last_status = "partial"
            self.last_message = "; ".join(dict.fromkeys(warnings))
        return _deduplicate(records)

    def _study_records(
        self,
        study: dict[str, Any],
        context: TherapyContext,
        retrieved_at: str,
    ) -> list[TherapyAssociation]:
        """Role: Compute study records for the therapy workflow."""
        protocol = study.get("protocolSection") or {}
        identification = protocol.get("identificationModule") or {}
        status = protocol.get("statusModule") or {}
        design = protocol.get("designModule") or {}
        conditions = protocol.get("conditionsModule") or {}
        arms = protocol.get("armsInterventionsModule") or {}
        eligibility = protocol.get("eligibilityModule") or {}
        description = protocol.get("descriptionModule") or {}
        nct_id = _optional_text(identification.get("nctId"))
        if not nct_id:
            return []
        sections = _study_text_sections(
            identification, description, eligibility, conditions, arms
        )
        gene_locations = tuple(
            name
            for name, text in sections.items()
            if _gene_symbol_mentioned(context.gene, text)
        )
        if not gene_locations:
            return []
        change_locations = {
            change: tuple(
                name
                for name, text in sections.items()
                if _protein_change_mentioned(change, text)
            )
            for change in context.protein_changes
        }
        matched_changes = tuple(
            change for change, locations in change_locations.items() if locations
        )
        match_locations = tuple(
            sorted(
                set(gene_locations).union(
                    location
                    for change in matched_changes
                    for location in change_locations[change]
                )
            )
        )
        phases = tuple(str(value) for value in design.get("phases") or [])
        disease_names = tuple(
            sorted({str(value) for value in conditions.get("conditions") or [] if value})
        )
        arm_roles = _trial_intervention_roles(arms.get("armGroups") or [])
        interventions: list[tuple[dict[str, Any], tuple[str, ...]]] = []
        for item in arms.get("interventions") or []:
            if (
                str(item.get("type") or "").upper()
                not in {"DRUG", "BIOLOGICAL", "COMBINATION_PRODUCT", "GENETIC"}
                or not item.get("name")
            ):
                continue
            therapy_name = str(item["name"])
            roles = arm_roles.get(therapy_name.casefold(), ("unspecified",))
            if (
                any(
                    role in {"placebo_comparator", "sham_comparator", "no_intervention"}
                    for role in roles
                )
                or "placebo" in therapy_name.casefold()
            ):
                continue
            interventions.append((item, roles))
        records: list[TherapyAssociation] = []
        for intervention, roles in interventions:
            therapy_name = str(intervention["name"])
            records.append(
                TherapyAssociation(
                    layer=self.layer,
                    source=self.name,
                    source_record_id=f"{nct_id}:{therapy_name}",
                    gene=context.gene,
                    association_scope=(
                        "trial_variant_mention"
                        if matched_changes
                        else "trial_gene_mention"
                    ),
                    therapy_name=therapy_name,
                    therapy_components=(therapy_name,),
                    therapy_status="clinical_trial",
                    protein_change=(matched_changes[0] if len(matched_changes) == 1 else None),
                    molecular_profile=(
                        " + ".join(matched_changes) if len(matched_changes) > 1 else None
                    ),
                    disease_names=disease_names,
                    evidence_type="clinical_trial_registry",
                    evidence_level=" | ".join(phases) or None,
                    trial_id=nct_id,
                    trial_phase=" | ".join(phases) or None,
                    trial_status=_optional_text(status.get("overallStatus")),
                    trial_role=" | ".join(roles),
                    description=_optional_text(identification.get("briefTitle")),
                    source_url=f"https://clinicaltrials.gov/study/{nct_id}",
                    retrieved_at=retrieved_at,
                    source_version=self.source_version,
                    license=self.license_name,
                    details=tuple(
                        item
                        for item in (
                            ("match_locations", " | ".join(match_locations)),
                            (
                                "matched_variants",
                                " | ".join(matched_changes),
                            ),
                        )
                        if item[1]
                    ),
                )
            )
        return records


class OpenTargetsTherapyProvider:
    """Gene-target drug landscape; results remain mechanism/discovery evidence."""

    name = "opentargets"
    layer = TRIAL_AND_TARGET_DISCOVERY
    source_version: str | None = None
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    license_name = "CC0-1.0 (input-source licences may also apply)"

    _SEARCH_QUERY = """
      query($query: String!) {
        search(queryString: $query, entityNames: ["target"], page: {index: 0, size: 10}) {
          hits { id name entity }
        }
      }
    """
    _DRUG_QUERY = """
      query($id: String!) {
        target(ensemblId: $id) {
          id approvedSymbol approvedName
          drugAndClinicalCandidates {
            count
            rows {
              id maxClinicalStage
              drug { id name drugType maximumClinicalStage }
              diseases {
                diseaseFromSource
                disease { id name }
              }
            }
          }
        }
      }
    """

    def __init__(
        self,
        *,
        http: PostJsonClient | None = None,
        max_drugs: int = 100,
    ) -> None:
        """Role: Compute init for the therapy workflow."""
        if max_drugs < 1:
            raise ValueError("max_drugs must be positive")
        self.http = http or JsonHttpClient(
            user_agent="rare-oncology/2.1 (Open Targets therapy lookup)",
            retries=2,
        )
        self.max_drugs = max_drugs

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {"url": self.url, "max_drugs": self.max_drugs}

    def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
        """Role: Collect source-specific records for the current analysis context."""
        search = self._graphql(self._SEARCH_QUERY, {"query": context.gene})
        hits = search.get("data", {}).get("search", {}).get("hits") or []
        exact = next(
            (
                hit
                for hit in hits
                if str(hit.get("name") or "").upper() == context.gene
                and hit.get("entity") == "target"
            ),
            None,
        )
        if not exact:
            return []
        payload = self._graphql(self._DRUG_QUERY, {"id": str(exact["id"])})
        target = payload.get("data", {}).get("target")
        if not target:
            return []
        rows = (target.get("drugAndClinicalCandidates") or {}).get("rows") or []
        retrieved_at = _utc_now()
        records: list[TherapyAssociation] = []
        for row in rows[: self.max_drugs]:
            drug = row.get("drug") or {}
            if not drug.get("name"):
                continue
            disease_ids: set[str] = set()
            disease_names: set[str] = set()
            for item in row.get("diseases") or []:
                disease = item.get("disease") or {}
                if disease.get("id"):
                    disease_ids.add(str(disease["id"]))
                name = disease.get("name") or item.get("diseaseFromSource")
                if name:
                    disease_names.add(str(name))
            stage = _optional_text(row.get("maxClinicalStage"))
            drug_id = str(drug.get("id") or row.get("id"))
            records.append(
                TherapyAssociation(
                    layer=self.layer,
                    source=self.name,
                    source_record_id=str(row.get("id") or drug_id),
                    gene=context.gene,
                    association_scope="gene_target",
                    therapy_name=str(drug["name"]),
                    therapy_components=(str(drug["name"]),),
                    therapy_status="mechanism_target",
                    therapy_id=drug_id,
                    disease_ids=tuple(sorted(disease_ids)),
                    disease_names=tuple(sorted(disease_names)),
                    evidence_type="target_drug_clinical_landscape",
                    evidence_level=stage,
                    description=_optional_text(drug.get("drugType")),
                    source_url=f"https://platform.opentargets.org/target/{target['id']}",
                    retrieved_at=retrieved_at,
                    source_version=self.source_version,
                    license=self.license_name,
                    details=(("maximum_clinical_stage", stage),) if stage else (),
                )
            )
        return _deduplicate(records)

    def _graphql(self, query: str, variables: dict[str, Any]) -> dict[str, Any]:
        """Role: Compute graphql for the therapy workflow."""
        payload = self.http.post_json(
            self.url, {"query": query, "variables": variables}
        )
        errors = payload.get("errors")
        if errors:
            message = "; ".join(
                str(item.get("message", item)) if isinstance(item, dict) else str(item)
                for item in errors
            )
            raise ApiError(f"Open Targets GraphQL: {message}")
        return payload


class RegulatoryTableTherapyProvider:
    """Import an official FDA/EMA-style table without scraping or inferring approval."""

    layer = REGULATORY_THERAPY

    def __init__(
        self,
        source: str,
        path: str | Path,
        *,
        jurisdiction: str | None = None,
        source_version: str | None = None,
        license_name: str | None = None,
    ) -> None:
        """Role: Compute init for the therapy workflow."""
        self.name = source.strip().lower()
        if not self.name:
            raise ValueError("source cannot be empty")
        self.path = Path(path)
        self.jurisdiction = jurisdiction
        self.source_version = source_version
        self.license_name = license_name
        self.url = self.path.resolve().as_uri()

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {
            "source": self.name,
            "path": self.path,
            "jurisdiction": self.jurisdiction,
            "source_version": self.source_version,
            "license": self.license_name,
        }

    def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
        """Role: Collect source-specific records for the current analysis context."""
        rows = _read_table(self.path)
        retrieved_at = _utc_now()
        records: list[TherapyAssociation] = []
        for index, raw in enumerate(rows, start=1):
            row = {_column_key(key): value for key, value in raw.items()}
            gene = _first(row, "gene", "gene_symbol", "symbol", "hgnc_symbol")
            if not gene:
                raise ValueError(f"Row {index} has no gene in {self.path}")
            if gene.upper() != context.gene:
                continue
            therapy = _first(row, "therapy_name", "drug_name", "drug", "product")
            if not therapy:
                raise ValueError(f"Row {index} has no therapy in {self.path}")
            protein_raw = _first(
                row, "protein_change", "hgvsp", "hgvs_p", "variant"
            )
            protein = normalize_protein_change(protein_raw)
            biomarker = _first(row, "biomarker_text", "biomarker", "alteration")
            biomarker_profile = biomarker or (protein_raw if not protein else None)
            record_id = _first(
                row,
                "source_record_id",
                "record_id",
                "application_number",
                "id",
            )
            source_url = _first(row, "source_url", "label_url", "url")
            if not record_id and not source_url:
                raise ValueError(
                    f"Row {index} has neither an identifier nor a regulatory URL in {self.path}"
                )
            record_id = record_id or f"url:{index}"
            status = _first(
                row, "therapy_status", "regulatory_status", "approval_status", "status"
            ) or "regulatory_status_not_provided"
            components = tuple(
                value.strip()
                for value in re.split(r"\s*[+|;]\s*", therapy)
                if value.strip()
            )
            records.append(
                TherapyAssociation(
                    layer=self.layer,
                    source=self.name,
                    source_record_id=record_id,
                    gene=context.gene,
                    association_scope=(
                        "variant_biomarker"
                        if protein
                        else "biomarker"
                        if biomarker_profile
                        else "indication"
                    ),
                    therapy_name=therapy,
                    therapy_components=components or (therapy,),
                    therapy_status=_snake_case(status) or "label_biomarker",
                    therapy_id=_first(row, "therapy_id", "drug_id", "chembl_id"),
                    protein_change=protein,
                    molecular_profile=biomarker_profile,
                    disease_ids=_split_values(
                        _first(row, "disease_id", "condition_id", "orpha_code")
                    ),
                    disease_names=_split_values(
                        _first(row, "disease_name", "disease", "indication")
                    ),
                    effect=_snake_case(
                        _first(row, "effect", "clinical_effect", "response")
                    ),
                    evidence_type=_first(row, "evidence_type", "labeling_section"),
                    evidence_level=_first(row, "evidence_level"),
                    jurisdiction=(
                        _first(row, "jurisdiction", "region") or self.jurisdiction
                    ),
                    citation_id=_first(row, "citation_id", "label_id"),
                    description=_first(row, "description", "biomarker_text", "notes"),
                    source_url=source_url,
                    retrieved_at=retrieved_at,
                    source_version=(
                        _first(row, "source_version", "version") or self.source_version
                    ),
                    license=_first(row, "license") or self.license_name,
                )
            )
        return _deduplicate(records)


def find_therapy_associations(
    result: AnalysisResult,
    *,
    providers: Iterable[TherapyProvider] | None = None,
    strict: bool = False,
) -> TherapyResult:
    """Enrich a mutation analysis without changing its molecular conclusions."""

    context = TherapyContext.from_analysis(result)
    selected: tuple[TherapyProvider, ...] = tuple(
        providers
        if providers is not None
        else (
            CivicTherapyProvider(),
            ClinicalTrialsTherapyProvider(),
            OpenTargetsTherapyProvider(),
        )
    )
    records, statuses = collect_therapy_associations(context, selected, strict=strict)
    by_layer = {
        VARIANT_THERAPY_EVIDENCE: tuple(
            record for record in records if record.layer == VARIANT_THERAPY_EVIDENCE
        ),
        REGULATORY_THERAPY: tuple(
            record for record in records if record.layer == REGULATORY_THERAPY
        ),
        TRIAL_AND_TARGET_DISCOVERY: tuple(
            record
            for record in records
            if record.layer == TRIAL_AND_TARGET_DISCOVERY
        ),
    }
    return TherapyResult(
        gene=context.gene,
        generated_at=_utc_now(),
        protein_changes=context.protein_changes,
        variant_therapy_evidence=by_layer[VARIANT_THERAPY_EVIDENCE],
        regulatory_therapy=by_layer[REGULATORY_THERAPY],
        trial_and_target_discovery=by_layer[TRIAL_AND_TARGET_DISCOVERY],
        source_statuses=statuses,
    )


def collect_therapy_associations(
    context: TherapyContext,
    providers: Sequence[TherapyProvider],
    *,
    strict: bool = False,
) -> tuple[tuple[TherapyAssociation, ...], tuple[SourceStatus, ...]]:
    """Role: Collect and consolidate therapy associations records."""
    records: list[TherapyAssociation] = []
    statuses: list[SourceStatus] = []
    allowed_layers = {
        VARIANT_THERAPY_EVIDENCE,
        REGULATORY_THERAPY,
        TRIAL_AND_TARGET_DISCOVERY,
    }
    for provider in providers:
        try:
            if provider.layer not in allowed_layers:
                raise ValueError(
                    f"{provider.name} declared an unknown therapy layer: "
                    f"{provider.layer}"
                )
            collected = provider.collect(context)
            wrong_layer = [record for record in collected if record.layer != provider.layer]
            if wrong_layer:
                raise ValueError(
                    f"{provider.name} returned records outside {provider.layer}"
                )
            records.extend(collected)
            statuses.append(
                SourceStatus(
                    source=provider.name,
                    status=getattr(provider, "last_status", "ok"),
                    evidence_layer=provider.layer,
                    record_count=len(collected),
                    message=getattr(provider, "last_message", None),
                    source_version=provider.source_version,
                    url=provider.url,
                )
            )
        except (ApiError, OSError, ValueError, KeyError, TypeError) as exc:
            if strict:
                raise
            statuses.append(
                SourceStatus(
                    source=provider.name,
                    status="error",
                    evidence_layer=provider.layer,
                    record_count=0,
                    message=str(exc),
                    source_version=provider.source_version,
                    url=provider.url,
                )
            )
    return tuple(_deduplicate(records)), tuple(statuses)


def _deduplicate(records: Iterable[TherapyAssociation]) -> list[TherapyAssociation]:
    """Role: Deduplicate therapy associations while preserving provenance."""
    unique: dict[tuple[object, ...], TherapyAssociation] = {}
    for record in records:
        key = (
            record.layer,
            record.source,
            record.source_record_id,
            record.therapy_name,
            record.protein_change,
            record.molecular_profile,
            record.disease_ids,
            record.disease_names,
            record.effect,
            record.therapy_status,
            record.jurisdiction,
        )
        unique.setdefault(key, record)
    return sorted(
        unique.values(),
        key=lambda item: (
            item.layer,
            item.source,
            item.protein_change or "",
            item.therapy_name.casefold(),
            item.source_record_id,
        ),
    )


def _read_table(path: Path) -> list[dict[str, Any]]:
    """Role: Compute read table for the therapy workflow."""
    if not path.exists():
        raise ValueError(f"File does not exist: {path}")
    if path.suffix.casefold() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            payload = payload.get("records", payload.get("data"))
        if not isinstance(payload, list) or not all(
            isinstance(item, dict) for item in payload
        ):
            raise ValueError("Therapy JSON must contain a list of objects")
        return payload
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel_tab if path.suffix.casefold() == ".tsv" else csv.excel
        return list(csv.DictReader(handle, dialect=dialect))


def _column_key(value: str) -> str:
    """Role: Compute column key for the therapy workflow."""
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().casefold()).strip("_")


def _first(row: dict[str, Any], *keys: str) -> str | None:
    """Role: Compute first for the therapy workflow."""
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def _split_values(value: str | None) -> tuple[str, ...]:
    """Role: Compute split values for the therapy workflow."""
    if not value:
        return ()
    return tuple(
        sorted({item.strip() for item in re.split(r"\s*[|;]\s*", value) if item.strip()})
    )


def _optional_text(value: Any) -> str | None:
    """Role: Convert optional text input into its typed value."""
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _optional_float(value: Any) -> float | None:
    """Role: Convert optional float input into its typed value."""
    if value is None or value == "":
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _snake_case(value: Any) -> str | None:
    """Role: Compute snake case for the therapy workflow."""
    text = _optional_text(value)
    if not text:
        return None
    text = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", text)
    return re.sub(r"[^a-z0-9]+", "_", text.casefold()).strip("_")


def _civic_effect(value: Any) -> str | None:
    """Role: Compute civic effect for the therapy workflow."""
    text = _optional_text(value)
    if not text:
        return None
    normalized = re.sub(r"[^A-Z]", "", text.upper())
    return {
        "SENSITIVITYRESPONSE": "sensitivity_response",
        "RESISTANCE": "resistance",
        "REDUCEDSENSITIVITY": "reduced_sensitivity",
        "ADVERSERESPONSE": "adverse_response",
        "NORESPONSE": "no_response",
    }.get(normalized, _snake_case(text))


def _absolute_civic_url(value: Any) -> str | None:
    """Role: Compute absolute civic url for the therapy workflow."""
    text = _optional_text(value)
    if not text:
        return None
    return text if text.startswith(("http://", "https://")) else f"https://civicdb.org{text}"


def _protein_change_mentioned(change: str, text: str) -> bool:
    """Role: Compute protein change mentioned for the therapy workflow."""
    return bool(
        re.search(
            rf"(?<![A-Z0-9])(?:P\.)?{re.escape(change)}(?![A-Z0-9])",
            text,
            flags=re.IGNORECASE,
        )
    )


def _study_nct_id(study: dict[str, Any]) -> str | None:
    """Role: Compute study nct id for the therapy workflow."""
    return _optional_text(
        (study.get("protocolSection") or {})
        .get("identificationModule", {})
        .get("nctId")
    )


def _study_text_sections(
    identification: dict[str, Any],
    description: dict[str, Any],
    eligibility: dict[str, Any],
    conditions: dict[str, Any],
    arms: dict[str, Any],
) -> dict[str, str]:
    """Role: Compute study text sections for the therapy workflow."""
    return {
        "title": "\n".join(
            str(value or "")
            for value in (
                identification.get("briefTitle"),
                identification.get("officialTitle"),
            )
        ),
        "summary": "\n".join(
            str(value or "")
            for value in (
                description.get("briefSummary"),
                description.get("detailedDescription"),
            )
        ),
        "eligibility": str(eligibility.get("eligibilityCriteria") or ""),
        "conditions": "\n".join(
            str(value or "")
            for value in (
                *(conditions.get("conditions") or []),
                *(conditions.get("keywords") or []),
            )
        ),
        "interventions": "\n".join(
            str(value or "")
            for intervention in arms.get("interventions") or []
            for value in (
                intervention.get("name"),
                intervention.get("description"),
                *(intervention.get("otherNames") or []),
            )
        ),
    }


def _next_graphql_cursor(
    page: dict[str, Any],
    current: str | None,
    seen: set[str],
    label: str,
) -> str | None:
    """Role: Compute next graphql cursor for the therapy workflow."""
    if not page.get("hasNextPage"):
        return None
    cursor = _optional_text(page.get("endCursor"))
    if not cursor:
        raise ApiError(f"{label}: missing cursor for the next page")
    if cursor == current or cursor in seen:
        raise ApiError(f"{label}: repeated pagination cursor")
    seen.add(cursor)
    return cursor


def _gene_symbol_mentioned(gene: str, text: str) -> bool:
    """Require the canonical uppercase token to avoid short-word false matches."""

    return bool(
        re.search(
            rf"(?<![A-Z0-9]){re.escape(gene.upper())}(?![A-Z0-9])",
            text,
        )
    )


def _trial_intervention_roles(
    arm_groups: Sequence[dict[str, Any]],
) -> dict[str, tuple[str, ...]]:
    """Role: Compute trial intervention roles for the therapy workflow."""
    roles: dict[str, set[str]] = {}
    for arm in arm_groups:
        role = _snake_case(arm.get("type")) or "unspecified"
        for value in arm.get("interventionNames") or []:
            name = str(value).strip()
            if ":" in name:
                name = name.split(":", 1)[1].strip()
            if name:
                roles.setdefault(name.casefold(), set()).add(role)
    return {key: tuple(sorted(values)) for key, values in roles.items()}


def _utc_now() -> str:
    """Role: Return the current UTC timestamp in canonical form."""
    return datetime.now(timezone.utc).isoformat()
