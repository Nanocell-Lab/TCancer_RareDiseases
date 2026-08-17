"""Layered evidence adapters for rare-disease and oncology interpretation.

The module deliberately keeps variant assertions, case observations, population
frequency, gene-disease validity and oncology interpretation as different
evidence layers.  Combining those categories into one pathogenicity flag would
lose provenance and overstate weaker sources.
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Protocol
from xml.etree import ElementTree

from .http import ApiError, JsonHttpClient
from .models import CancerObservation, EvidenceRecord, RareVariant, SourceStatus
from .normalization import normalize_protein_change


CLINVAR_URL = "https://www.ncbi.nlm.nih.gov/clinvar/"
GDC_URL = "https://api.gdc.cancer.gov/ssm_occurrences"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
PANELAPP_GENE = "https://panelapp.genomicsengland.co.uk/api/v1/entities/{gene}/"

VARIANT_ASSERTION = "variant_assertion"
CASE_OBSERVATION = "case_observation"
POPULATION_FREQUENCY = "population_frequency"
GENE_DISEASE_VALIDITY = "gene_disease_validity"
DISEASE_ONTOLOGY = "disease_ontology"
ONCOLOGY_EVIDENCE = "oncology_evidence"


@dataclass(frozen=True)
class EvidenceContext:
    gene: str
    rare_variants: tuple[RareVariant, ...]
    cancer_observations: tuple[CancerObservation, ...]


class EvidenceProvider(Protocol):
    @property
    def name(self) -> str:
        """Role: Expose the stable evidence-provider name."""
        ...

    @property
    def evidence_layer(self) -> str:
        """Role: Expose the evidence layer produced by the provider."""
        ...

    @property
    def url(self) -> str | None:
        """Role: Expose the provider endpoint when one is configured."""
        ...

    @property
    def source_version(self) -> str | None:
        """Role: Expose the observed source release or version."""
        ...

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        """Role: Collect normalized evidence records for one gene context."""
        ...


class ClinVarEvidenceProvider:
    name: str = "clinvar"
    evidence_layer: str = VARIANT_ASSERTION
    url: str | None = CLINVAR_URL
    source_version: str | None = None

    def __init__(
        self,
        *,
        details_queried: bool = True,
        unresolved_condition_count: int = 0,
    ) -> None:
        """Role: Compute init for the evidence workflow."""
        self.details_queried = details_queried
        self.unresolved_condition_count = unresolved_condition_count
        self.last_status = "ok"
        self.last_message: str | None = None

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {
            "url": CLINVAR_URL,
            "details_queried": self.details_queried,
            "unresolved_condition_count": self.unresolved_condition_count,
        }

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        """Role: Collect source-specific records for the current analysis context."""
        unresolved = self.unresolved_condition_count + sum(
            1
            for variant in context.rare_variants
            if self.details_queried and not variant.condition_assertions
        )
        if not self.details_queried:
            self.last_status = "ok"
            self.last_message = "VCV/RCV/SCV details were not queried (--no-clinvar-details)"
        elif unresolved:
            self.last_status = "partial"
            self.last_message = (
                f"{unresolved} VCV records lack a resolvable RCV-Orphanet relationship; "
                "aggregate classification was not assigned to the condition"
            )
        else:
            self.last_status = "ok"
            self.last_message = None
        retrieved_at = _utc_now()
        records: list[EvidenceRecord] = []
        for variant in context.rare_variants:
            identity = spdi_identity(variant.canonical_spdi)
            condition_rows: list[
                tuple[str, str, str, str | None, str | None, str | None]
            ] = []
            if variant.condition_assertions:
                for assertion in variant.condition_assertions:
                    condition_rows.extend(
                        (
                            assertion.disease,
                            code,
                            assertion.rcv_accession,
                            assertion.classification,
                            assertion.review_status,
                            assertion.last_evaluated,
                        )
                        for code in assertion.orphanet_codes
                    )
            else:
                condition_rows.extend(
                    (disease, code, variant.clinvar_accession, None, None, None)
                    for disease, code in variant.condition_orphanet_pairs
                )
            for (
                disease,
                orphanet_code,
                condition_accession,
                classification,
                review_status,
                last_evaluated,
            ) in condition_rows:
                details = {
                    "clinvar_id": variant.clinvar_id,
                    "last_evaluated": last_evaluated or "",
                    "detail_query": "queried" if self.details_queried else "not_queried",
                }
                records.append(
                    EvidenceRecord(
                        source=self.name,
                        source_record_id=(
                            f"{condition_accession}|ORPHA:{orphanet_code}"
                            if orphanet_code
                            else condition_accession
                        ),
                        evidence_layer=self.evidence_layer,
                        gene=context.gene,
                        retrieved_at=retrieved_at,
                        license="NCBI data usage policies",
                        variant_key=identity.get("variant_key"),
                        assembly=identity.get("assembly"),
                        chromosome=identity.get("chromosome"),
                        position=identity.get("position"),
                        reference=identity.get("reference"),
                        alternate=identity.get("alternate"),
                        protein_change=variant.protein_change,
                        disease_id=(
                            f"ORPHA:{orphanet_code}" if orphanet_code else None
                        ),
                        disease_name=disease or None,
                        classification=classification,
                        review_level=review_status,
                        details=_details(details),
                    )
                )
            # Aggregate VCV classification and flattened SCV submitters remain
            # explicitly variant-level and are never copied to a disease.
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=(
                        f"{variant.clinvar_accession}|variant_provenance"
                    ),
                    evidence_layer=self.evidence_layer,
                    gene=context.gene,
                    retrieved_at=retrieved_at,
                    license="NCBI data usage policies",
                    variant_key=identity.get("variant_key"),
                    assembly=identity.get("assembly"),
                    chromosome=identity.get("chromosome"),
                    position=identity.get("position"),
                    reference=identity.get("reference"),
                    alternate=identity.get("alternate"),
                    protein_change=variant.protein_change,
                    classification=variant.germline_classification or None,
                    review_level=variant.review_status,
                    details=_details(
                        {
                            "clinvar_id": variant.clinvar_id,
                            "rcv_accessions": "|".join(variant.rcv_accessions),
                            "scv_accessions": "|".join(variant.scv_accessions),
                            "submitters": "|".join(variant.submitters),
                            "number_of_submissions": _text(
                                variant.number_of_submissions
                            ),
                            "number_of_submitters": _text(
                                variant.number_of_submitters
                            ),
                            "last_evaluated": variant.last_evaluated or "",
                            "date_last_updated": variant.date_last_updated or "",
                            "variation_type": variant.variation_type or "",
                            "detail_query": (
                                "queried" if self.details_queried else "not_queried"
                            ),
                            "scope": "variant_level_not_condition_attributed",
                        }
                    ),
                )
            )
        return records


class GdcEvidenceProvider:
    name: str = "gdc"
    evidence_layer: str = CASE_OBSERVATION
    url: str | None = GDC_URL

    def __init__(self, source_version: str | None = None) -> None:
        """Role: Compute init for the evidence workflow."""
        self.source_version = source_version

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {"url": GDC_URL, "source_version": self.source_version}

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        """Role: Collect source-specific records for the current analysis context."""
        retrieved_at = _utc_now()
        records: list[EvidenceRecord] = []
        for observation in context.cancer_observations:
            identity = gdc_identity(observation.genomic_dna_change)
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=(
                        f"{observation.ssm_id}|{observation.case_id}"
                    ),
                    evidence_layer=self.evidence_layer,
                    gene=context.gene,
                    retrieved_at=retrieved_at,
                    source_version=self.source_version,
                    license="GDC data usage policies",
                    variant_key=identity.get("variant_key"),
                    assembly=identity.get("assembly"),
                    chromosome=identity.get("chromosome"),
                    position=identity.get("position"),
                    reference=identity.get("reference"),
                    alternate=identity.get("alternate"),
                    transcript=observation.transcript_id,
                    protein_change=observation.protein_change,
                    details=_details(
                        {
                            "case_id": observation.case_id,
                            "patient_id": observation.patient_id,
                            "project_id": observation.project_id,
                            "ssm_id": observation.ssm_id,
                        }
                    ),
                )
            )
        return records


class GnomadEvidenceProvider:
    """Annotate candidate GRCh38 SPDI alleles with gnomAD joint frequency."""

    name: str = "gnomad"
    evidence_layer: str = POPULATION_FREQUENCY
    url: str | None = GNOMAD_API
    source_version: str | None = "gnomAD v4 (latest v4 API release)"

    def __init__(
        self,
        *,
        http: JsonHttpClient | None = None,
        max_variants: int = 25,
    ) -> None:
        """Role: Compute init for the evidence workflow."""
        if max_variants < 1 or max_variants > 25:
            raise ValueError("max_variants must be between 1 and 25")
        self.max_variants = max_variants
        self.http = http or JsonHttpClient(user_agent="rare-oncology/2.1")
        self.last_message: str | None = None
        self.last_status = "ok"

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {"url": GNOMAD_API, "max_variants": self.max_variants}

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        """Role: Collect source-specific records for the current analysis context."""
        self.last_message = None
        self.last_status = "ok"
        warnings: list[str] = []
        identities: dict[str, dict[str, Any]] = {}
        unmappable = 0
        for variant in context.rare_variants:
            identity = spdi_identity(variant.canonical_spdi)
            variant_id = identity.get("gnomad_variant_id")
            if variant_id:
                identities[str(variant_id)] = identity
            elif variant.canonical_spdi:
                unmappable += 1
        if unmappable:
            warnings.append(
                f"{unmappable} alleles skipped: gnomAD queries only unambiguous "
                "GRCh38 A/C/G/T SNVs"
            )
        selected = sorted(identities)[: self.max_variants]
        if len(identities) > len(selected):
            warnings.append(
                f"Queried {len(selected)} of {len(identities)} alleles; "
                "increasing the limit requires batch processing."
            )
        if not selected:
            self.last_message = " ".join(warnings) or None
            return []

        aliases = []
        for index, variant_id in enumerate(selected):
            aliases.append(
                f'v{index}: variant(variantId: "{variant_id}", dataset: gnomad_r4) '
                "{ variant_id reference_genome chrom pos ref alt "
                "joint { ac an } exome { ac an } genome { ac an } }"
            )
        response = self.http.post_json(GNOMAD_API, {"query": "query {" + " ".join(aliases) + "}"})
        errors = response.get("errors") or []
        unexpected_errors = [
            error
            for error in errors
            if str(error.get("message", "")).casefold() != "variant not found"
        ]
        if unexpected_errors:
            raise ApiError(f"gnomAD GraphQL: {unexpected_errors}")
        data = response.get("data") or {}
        retrieved_at = _utc_now()
        records: list[EvidenceRecord] = []
        missing = 0
        for index, variant_id in enumerate(selected):
            item = data.get(f"v{index}")
            if not isinstance(item, dict):
                missing += 1
                identity = identities[variant_id]
                records.append(
                    EvidenceRecord(
                        source=self.name,
                        source_record_id=variant_id,
                        evidence_layer=self.evidence_layer,
                        gene=context.gene,
                        retrieved_at=retrieved_at,
                        source_version=self.source_version,
                        license="gnomAD terms of use",
                        variant_key=identity.get("variant_key"),
                        assembly=identity.get("assembly"),
                        chromosome=identity.get("chromosome"),
                        position=identity.get("position"),
                        reference=identity.get("reference"),
                        alternate=identity.get("alternate"),
                        details=(("observed", "false"),),
                    )
                )
                continue
            joint = item.get("joint") or {}
            ac = _optional_int(joint.get("ac"))
            an = _optional_int(joint.get("an"))
            frequency = _frequency_from_counts(ac, an)
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=str(item.get("variant_id") or variant_id),
                    evidence_layer=self.evidence_layer,
                    gene=context.gene,
                    retrieved_at=retrieved_at,
                    source_version=self.source_version,
                    license="gnomAD terms of use",
                    variant_key=_variant_key(
                        str(item.get("reference_genome") or "GRCh38"),
                        str(item.get("chrom") or identities[variant_id]["chromosome"]),
                        _optional_int(item.get("pos")),
                        str(item.get("ref") or identities[variant_id]["reference"]),
                        str(item.get("alt") or identities[variant_id]["alternate"]),
                    ),
                    assembly=str(item.get("reference_genome") or "GRCh38"),
                    chromosome=str(item.get("chrom") or identities[variant_id]["chromosome"]),
                    position=_optional_int(item.get("pos")),
                    reference=str(item.get("ref") or identities[variant_id]["reference"]),
                    alternate=str(item.get("alt") or identities[variant_id]["alternate"]),
                    population_frequency=frequency,
                    details=_details(
                        {
                            "observed": "true",
                            "joint_ac": _text(ac),
                            "joint_an": _text(an),
                            "exome_ac": _nested_number(item, "exome", "ac"),
                            "exome_an": _nested_number(item, "exome", "an"),
                            "genome_ac": _nested_number(item, "genome", "ac"),
                            "genome_an": _nested_number(item, "genome", "an"),
                        }
                    ),
                )
            )
        if missing:
            warnings.append(
                f"{missing} of {len(selected)} alleles were not observed in gnomAD"
            )
        self.last_message = " ".join(warnings) or None
        return records


class PanelAppEvidenceProvider:
    name: str = "panelapp"
    evidence_layer: str = GENE_DISEASE_VALIDITY
    source_version: str | None = None

    def __init__(self, *, http: JsonHttpClient | None = None) -> None:
        """Role: Compute init for the evidence workflow."""
        self.http = http or JsonHttpClient(user_agent="rare-oncology/2.1")
        self.url: str | None = None

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {"url_template": PANELAPP_GENE}

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        """Role: Collect source-specific records for the current analysis context."""
        self.url = PANELAPP_GENE.format(gene=context.gene)
        payload = self.http.get_json(self.url, {})
        retrieved_at = _utc_now()
        records: list[EvidenceRecord] = []
        for row in payload.get("results", []):
            gene_data = row.get("gene_data") or {}
            symbol = str(
                gene_data.get("gene_symbol") or gene_data.get("hgnc_symbol") or ""
            ).upper()
            if symbol != context.gene:
                continue
            panel = row.get("panel") or {}
            panel_id = str(panel.get("id") or "unknown")
            version = str(panel.get("version") or "unknown")
            phenotypes = tuple(str(value) for value in row.get("phenotypes", []) if value)
            confidence = str(row.get("confidence_level") or "")
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=f"panel:{panel_id}:v{version}:{context.gene}",
                    evidence_layer=self.evidence_layer,
                    gene=context.gene,
                    retrieved_at=retrieved_at,
                    source_version=version,
                    license="Genomics England PanelApp terms of use",
                    disease_name=" | ".join(phenotypes) or None,
                    mode_of_inheritance=(
                        str(row.get("mode_of_inheritance"))
                        if row.get("mode_of_inheritance")
                        else None
                    ),
                    classification=(
                        {"3": "Green", "2": "Amber", "1": "Red"}.get(
                            confidence, confidence or None
                        )
                    ),
                    evidence_level="panel_gene_confidence",
                    details=_details(
                        {
                            "panel_id": panel_id,
                            "panel_name": str(panel.get("name") or ""),
                            "panel_version": version,
                            "entity_type": str(row.get("entity_type") or "gene"),
                            "evidence": "|".join(
                                str(value) for value in row.get("evidence", [])
                            ),
                        }
                    ),
                )
            )
        return records


class OrphadataXmlEvidenceProvider:
    """Read the official Orphadata gene-association XML snapshot."""

    name: str = "orphadata"
    evidence_layer: str = GENE_DISEASE_VALIDITY

    def __init__(self, path: str | Path) -> None:
        """Role: Compute init for the evidence workflow."""
        self.path = Path(path)
        self.url = self.path.resolve().as_uri()
        self.source_version: str | None = None
        self.license_name: str | None = None

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {"path": self.path}

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        """Role: Collect source-specific records for the current analysis context."""
        if not self.path.is_file():
            raise ValueError(f"Orphadata XML does not exist: {self.path}")
        try:
            iterator = ElementTree.iterparse(self.path, events=("start", "end"))
            _event, root = next(iterator)
            self.source_version = root.get("version") or root.get("date")
            records: list[EvidenceRecord] = []
            retrieved_at = _utc_now()
            for event, element in iterator:
                if event != "end" or element.tag != "Disorder":
                    continue
                orpha = _xml_text(element.find("OrphaCode"))
                disease = _xml_text(element.find("Name"))
                for association in element.findall(
                    "./DisorderGeneAssociationList/DisorderGeneAssociation"
                ):
                    symbol = _xml_text(association.find("./Gene/Symbol"))
                    if symbol.upper() != context.gene:
                        continue
                    association_type = _xml_text(
                        association.find("./DisorderGeneAssociationType/Name")
                    )
                    association_status = _xml_text(
                        association.find("./DisorderGeneAssociationStatus/Name")
                    )
                    references = {
                        _xml_text(item.find("Source")): _xml_text(
                            item.find("Reference")
                        )
                        for item in association.findall(
                            "./Gene/ExternalReferenceList/ExternalReference"
                        )
                    }
                    records.append(
                        EvidenceRecord(
                            source=self.name,
                            source_record_id=f"ORPHA:{orpha}|{context.gene}",
                            evidence_layer=self.evidence_layer,
                            gene=context.gene,
                            retrieved_at=retrieved_at,
                            source_version=self.source_version,
                            license="CC BY 4.0",
                            disease_id=f"ORPHA:{orpha}" if orpha else None,
                            disease_name=disease or None,
                            classification=association_type or None,
                            evidence_level=association_status or None,
                            details=_details(
                                {
                                    "source_of_validation": _xml_text(
                                        association.find("SourceOfValidation")
                                    ),
                                    "hgnc": references.get("HGNC", ""),
                                    "ensembl": references.get("Ensembl", ""),
                                    "omim": references.get("OMIM", ""),
                                    "swissprot": references.get("SwissProt", ""),
                                    "snapshot_date": root.get("date") or "",
                                }
                            ),
                        )
                    )
                element.clear()
            return records
        except ElementTree.ParseError as exc:
            raise ValueError("Invalid Orphadata XML") from exc


class LocalTableEvidenceProvider:
    """Read a normalized or common-column CSV/TSV/JSON export.

    This is the safe integration point for ClinGen ERepo, DECIPHER, LOVD,
    CIViC, GenCC and Orphadata exports whose bulk access or licences differ.
    The original source remains explicit in every record.
    """

    def __init__(
        self,
        source: str,
        path: str | Path,
        *,
        evidence_layer: str | None = None,
        source_version: str | None = None,
        license_name: str | None = None,
        assume_gene: bool = False,
    ) -> None:
        """Role: Compute init for the evidence workflow."""
        self.name: str = source.strip().lower()
        if not self.name:
            raise ValueError("source cannot be empty")
        self.path = Path(path)
        self.evidence_layer: str = evidence_layer or _default_layer(self.name)
        self.source_version: str | None = source_version
        self.license_name = license_name
        self.assume_gene = assume_gene
        self.url: str | None = self.path.resolve().as_uri()

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {
            "source": self.name,
            "path": self.path,
            "evidence_layer": self.evidence_layer,
            "source_version": self.source_version,
            "license": self.license_name,
            "assume_gene": self.assume_gene,
        }

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        """Role: Collect source-specific records for the current analysis context."""
        rows = _read_table(self.path)
        retrieved_at = _utc_now()
        records: list[EvidenceRecord] = []
        for index, raw in enumerate(rows, start=1):
            row = {_column_key(key): value for key, value in raw.items()}
            gene = _first(row, "gene", "gene_symbol", "symbol", "hgnc_symbol")
            if not gene and not self.assume_gene:
                raise ValueError(
                    f"Row {index} has no gene in {self.path}; it was not assigned to {context.gene}"
                )
            if gene and gene.upper() != context.gene:
                continue
            protein_raw = _first(
                row, "protein_change", "hgvsp", "hgvs_p", "variant", "variant_name"
            )
            protein = normalize_protein_change(protein_raw) if protein_raw else None
            assembly = _first(row, "assembly", "reference_genome", "genome_build")
            chromosome = _first(row, "chromosome", "chrom", "chr")
            position = _optional_int(_first(row, "position", "pos", "start"))
            reference = _first(row, "reference", "ref")
            alternate = _first(row, "alternate", "alt")
            variant_key = _first(row, "variant_key", "variant_id", "caid")
            if not variant_key:
                variant_key = _variant_key(
                    assembly, chromosome, position, reference, alternate
                )
            record_id = _first(
                row,
                "source_record_id",
                "record_id",
                "id",
                "accession",
                "allele_registry_id",
                "clinvar_variation_id",
            ) or f"row:{index}"
            records.append(
                EvidenceRecord(
                    source=self.name,
                    source_record_id=record_id,
                    evidence_layer=(
                        _first(row, "evidence_layer") or self.evidence_layer
                    ),
                    gene=context.gene,
                    retrieved_at=retrieved_at,
                    source_version=(
                        _first(row, "source_version", "version")
                        or self.source_version
                    ),
                    license=_first(row, "license") or self.license_name,
                    variant_key=variant_key,
                    assembly=assembly,
                    chromosome=chromosome,
                    position=position,
                    reference=reference,
                    alternate=alternate,
                    transcript=_first(row, "transcript", "transcript_id"),
                    hgvs_c=_first(row, "hgvs_c", "hgvsc"),
                    protein_change=protein,
                    disease_id=_first(
                        row, "disease_id", "condition_id", "orphanet_code", "orphacode"
                    ),
                    disease_name=_first(
                        row, "disease_name", "disease", "condition", "phenotype"
                    ),
                    mode_of_inheritance=_first(
                        row, "mode_of_inheritance", "inheritance", "moi"
                    ),
                    classification=_first(
                        row, "classification", "clinical_significance", "assertion"
                    ),
                    review_level=_first(
                        row, "review_level", "review_status", "expert_panel"
                    ),
                    evidence_level=_first(
                        row, "evidence_level", "confidence", "confidence_level"
                    ),
                    population_frequency=_validated_frequency(
                        _first(row, "population_frequency", "allele_frequency", "af")
                    ),
                    details=(
                        (("gene_scope_assumption", context.gene),)
                        if not gene and self.assume_gene
                        else ()
                    ),
                )
            )
        return records


def collect_evidence(
    context: EvidenceContext,
    providers: Sequence[EvidenceProvider],
    *,
    strict: bool = False,
) -> tuple[tuple[EvidenceRecord, ...], tuple[SourceStatus, ...]]:
    """Role: Collect and consolidate evidence records."""
    records: list[EvidenceRecord] = []
    statuses: list[SourceStatus] = []
    for provider in providers:
        try:
            batch = provider.collect(context)
            records.extend(batch)
            statuses.append(
                SourceStatus(
                    source=provider.name,
                    status=getattr(provider, "last_status", "ok"),
                    evidence_layer=provider.evidence_layer,
                    record_count=len(batch),
                    message=getattr(provider, "last_message", None),
                    source_version=getattr(provider, "source_version", None),
                    url=getattr(provider, "url", None),
                )
            )
        except (ApiError, OSError, ValueError, KeyError, TypeError) as exc:
            if strict:
                raise
            statuses.append(
                SourceStatus(
                    source=provider.name,
                    status="error",
                    evidence_layer=provider.evidence_layer,
                    record_count=0,
                    message=str(exc),
                    source_version=getattr(provider, "source_version", None),
                    url=getattr(provider, "url", None),
                )
            )
    return _deduplicate_records(records), tuple(statuses)


def spdi_identity(spdi: str | None) -> dict[str, Any]:
    """Role: Provide the spdi identity operation for the evidence workflow."""
    if not spdi:
        return {}
    match = re.fullmatch(r"([^:]+):(\d+):([ACGT]+):([ACGT]*)", spdi.upper())
    if not match:
        return {}
    accession, zero_based, reference, alternate = match.groups()
    chromosome = GRCH38_ACCESSIONS.get(accession)
    if not chromosome:
        return {}
    position = int(zero_based) + 1
    identity = {
        "assembly": "GRCh38",
        "chromosome": chromosome,
        "position": position,
        "reference": reference,
        "alternate": alternate,
        "variant_key": _variant_key(
            "GRCh38", chromosome, position, reference, alternate
        ),
    }
    if len(reference) == 1 and len(alternate) == 1:
        identity["gnomad_variant_id"] = (
            f"{chromosome}-{position}-{reference}-{alternate}"
        )
    return identity


def gdc_identity(change: str | None) -> dict[str, Any]:
    """Role: Provide the gdc identity operation for the evidence workflow."""
    if not change:
        return {}
    match = re.fullmatch(
        r"(?:chr)?([0-9XYM]+):g\.(\d+)([ACGTN]+)>([ACGTN]+)",
        change,
        flags=re.IGNORECASE,
    )
    if not match:
        return {}
    chromosome, position, reference, alternate = match.groups()
    chromosome = chromosome.upper().replace("M", "MT")
    numeric_position = int(position)
    return {
        "assembly": "GRCh38",
        "chromosome": chromosome,
        "position": numeric_position,
        "reference": reference.upper(),
        "alternate": alternate.upper(),
        "variant_key": _variant_key(
            "GRCh38",
            chromosome,
            numeric_position,
            reference.upper(),
            alternate.upper(),
        ),
    }


def _deduplicate_records(records: Iterable[EvidenceRecord]) -> tuple[EvidenceRecord, ...]:
    """Role: Deduplicate records while preserving provenance."""
    unique = {
        (
            record.source,
            record.source_record_id,
            record.evidence_layer,
            record.variant_key,
            record.disease_id,
        ): record
        for record in records
    }
    return tuple(
        sorted(
            unique.values(),
            key=lambda item: (
                item.evidence_layer,
                item.source,
                item.variant_key or "",
                item.source_record_id,
            ),
        )
    )


def _read_table(path: Path) -> list[dict[str, Any]]:
    """Role: Compute read table for the evidence workflow."""
    if not path.is_file():
        raise ValueError(f"Evidence file does not exist: {path}")
    if path.suffix.casefold() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            payload = payload.get("records", payload.get("data"))
        if not isinstance(payload, list) or not all(
            isinstance(item, dict) for item in payload
        ):
            raise ValueError("JSON must be a list or contain records/data")
        return payload
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel_tab if path.suffix.casefold() == ".tsv" else csv.excel
        return list(csv.DictReader(handle, dialect=dialect))


def _default_layer(source: str) -> str:
    """Role: Compute default layer for the evidence workflow."""
    return {
        "clingen": VARIANT_ASSERTION,
        "clingen_erepo": VARIANT_ASSERTION,
        "decipher": CASE_OBSERVATION,
        "lovd": VARIANT_ASSERTION,
        "mseqdr": VARIANT_ASSERTION,
        "civic": ONCOLOGY_EVIDENCE,
        "cosmic": ONCOLOGY_EVIDENCE,
        "gencc": GENE_DISEASE_VALIDITY,
        "panelapp": GENE_DISEASE_VALIDITY,
        "orphadata": GENE_DISEASE_VALIDITY,
        "omim": GENE_DISEASE_VALIDITY,
    }.get(source, VARIANT_ASSERTION)


def _column_key(value: str) -> str:
    """Role: Compute column key for the evidence workflow."""
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().casefold()).strip("_")


def _first(row: dict[str, Any], *keys: str) -> str | None:
    """Role: Compute first for the evidence workflow."""
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def _variant_key(
    assembly: str | None,
    chromosome: str | None,
    position: int | None,
    reference: str | None,
    alternate: str | None,
) -> str | None:
    """Role: Compute variant key for the evidence workflow."""
    if not all((assembly, chromosome, position, reference, alternate)):
        return None
    return f"{assembly}:{chromosome}:{position}:{reference}>{alternate}"


def _optional_int(value: Any) -> int | None:
    """Role: Convert optional int input into its typed value."""
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _validated_frequency(value: Any) -> float | None:
    """Role: Compute validated frequency for the evidence workflow."""
    if value is None or value == "":
        return None
    try:
        frequency = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid population frequency: {value!r}") from None
    if not math.isfinite(frequency) or not 0.0 <= frequency <= 1.0:
        raise ValueError(f"Population frequency outside [0,1]: {value!r}")
    return frequency


def _frequency_from_counts(ac: int | None, an: int | None) -> float | None:
    """Role: Compute frequency from counts for the evidence workflow."""
    if ac is None and an is None:
        return None
    if ac is None or an is None or an <= 0 or ac < 0 or ac > an:
        raise ApiError(f"Invalid gnomAD counts: AC={ac}, AN={an}")
    return ac / an


def _nested_number(item: dict[str, Any], group: str, field: str) -> str:
    """Role: Compute nested number for the evidence workflow."""
    return _text((item.get(group) or {}).get(field))


def _text(value: Any) -> str:
    """Role: Compute text for the evidence workflow."""
    return "" if value is None else str(value)


def _details(values: dict[str, str]) -> tuple[tuple[str, str], ...]:
    """Role: Compute details for the evidence workflow."""
    return tuple(sorted((key, value) for key, value in values.items() if value))


def _xml_text(element: ElementTree.Element | None) -> str:
    """Role: Compute xml text for the evidence workflow."""
    return (element.text or "").strip() if element is not None else ""


def _utc_now() -> str:
    """Role: Return the current UTC timestamp in canonical form."""
    return datetime.now(timezone.utc).isoformat()


GRCH38_ACCESSIONS = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.10": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "MT",
}
