"""Typed records returned by the analysis pipeline."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class RareVariant:
    gene: str
    protein_change: str
    raw_protein_change: str
    clinvar_id: str
    clinvar_accession: str
    diseases: tuple[str, ...]
    orphanet_codes: tuple[str, ...]
    condition_orphanet_pairs: tuple[tuple[str, str], ...]
    germline_classification: str
    molecular_consequences: tuple[str, ...]
    canonical_spdi: str | None


@dataclass(frozen=True)
class CancerObservation:
    gene: str
    protein_change: str
    raw_protein_change: str
    ssm_id: str
    case_id: str
    patient_id: str
    project_id: str
    transcript_id: str | None
    genomic_dna_change: str | None


@dataclass(frozen=True)
class SharedMutation:
    gene: str
    protein_change: str
    diseases: tuple[str, ...]
    orphanet_codes: tuple[str, ...]
    clinvar_ids: tuple[str, ...]
    clinvar_accessions: tuple[str, ...]
    germline_classifications: tuple[str, ...]
    canonical_spdis: tuple[str, ...]
    tcga_projects: tuple[str, ...]
    cancer_patient_count: int
    cancer_event_count: int
    genomic_match_count: int
    genomic_confirmations: tuple[str, ...]
    ssm_ids: tuple[str, ...]
    match_level: str


@dataclass(frozen=True)
class MutationPlotRecord:
    gene: str
    protein_change: str
    amino_acid_position: int | None
    in_rare_disease: bool
    in_oncology: bool
    diseases: tuple[str, ...]
    orphanet_codes: tuple[str, ...]
    germline_classifications: tuple[str, ...]
    cancer_patient_count: int
    cancer_event_count: int
    project_patient_counts: tuple[tuple[str, int], ...]
    match_level: str | None


@dataclass(frozen=True)
class AnalysisResult:
    gene: str
    generated_at: str
    gdc_program: str | None
    gdc_projects: tuple[str, ...]
    min_cancer_patients: int
    excluded_conditions: tuple[str, ...]
    clinvar_record_limit: int | None
    allowed_germline_classifications: tuple[str, ...] | None
    gdc_data_release: str | None
    rare_variant_count: int
    cancer_observation_count: int
    shared_mutation_count: int
    shared_mutations: tuple[SharedMutation, ...]
    mutation_plot_data: tuple[MutationPlotRecord, ...]
    sources: tuple[str, ...] = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        "https://api.gdc.cancer.gov/ssm_occurrences",
    )
    caveat: str = (
        "Coincidencias candidatas por cambio proteico canonico; confirmar identidad "
        "with genomic coordinates/alleles or transcript HGVS before clinical interpretation."
    )

    def to_dict(self) -> dict[str, Any]:
        """Role: Convert the domain object into a JSON-safe dictionary."""
        return asdict(self)
