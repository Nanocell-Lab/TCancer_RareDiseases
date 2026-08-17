"""Typed records returned by the analysis pipeline."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class ClinVarConditionAssertion:
    disease: str
    orphanet_codes: tuple[str, ...]
    rcv_accession: str
    classification: str | None
    review_status: str | None
    last_evaluated: str | None


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
    variation_type: str | None = None
    review_status: str | None = None
    last_evaluated: str | None = None
    date_last_updated: str | None = None
    number_of_submissions: int | None = None
    number_of_submitters: int | None = None
    rcv_accessions: tuple[str, ...] = ()
    scv_accessions: tuple[str, ...] = ()
    submitters: tuple[str, ...] = ()
    hgvs_expressions: tuple[str, ...] = ()
    condition_assertions: tuple[ClinVarConditionAssertion, ...] = ()


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
    disease_variant_counts: tuple[tuple[str, int], ...] = ()


@dataclass(frozen=True)
class EvidenceRecord:
    """A source assertion kept separate from the normalized mutation match."""

    source: str
    source_record_id: str
    evidence_layer: str
    gene: str
    retrieved_at: str
    source_version: str | None = None
    license: str | None = None
    variant_key: str | None = None
    assembly: str | None = None
    chromosome: str | None = None
    position: int | None = None
    reference: str | None = None
    alternate: str | None = None
    transcript: str | None = None
    hgvs_c: str | None = None
    protein_change: str | None = None
    disease_id: str | None = None
    disease_name: str | None = None
    mode_of_inheritance: str | None = None
    classification: str | None = None
    review_level: str | None = None
    evidence_level: str | None = None
    population_frequency: float | None = None
    prediction_score: float | None = None
    prediction_class: str | None = None
    details: tuple[tuple[str, str], ...] = ()


@dataclass(frozen=True)
class SourceStatus:
    source: str
    status: str
    evidence_layer: str
    record_count: int
    message: str | None = None
    source_version: str | None = None
    url: str | None = None


@dataclass(frozen=True)
class TherapyAssociation:
    """A therapy assertion whose semantic layer and provenance are explicit."""

    layer: str
    source: str
    source_record_id: str
    gene: str
    association_scope: str
    therapy_name: str
    therapy_components: tuple[str, ...]
    therapy_status: str
    retrieved_at: str
    therapy_id: str | None = None
    protein_change: str | None = None
    molecular_profile: str | None = None
    disease_ids: tuple[str, ...] = ()
    disease_names: tuple[str, ...] = ()
    effect: str | None = None
    evidence_type: str | None = None
    evidence_level: str | None = None
    evidence_direction: str | None = None
    evidence_rating: float | None = None
    jurisdiction: str | None = None
    trial_id: str | None = None
    trial_phase: str | None = None
    trial_status: str | None = None
    trial_role: str | None = None
    citation_id: str | None = None
    description: str | None = None
    source_url: str | None = None
    source_version: str | None = None
    license: str | None = None
    details: tuple[tuple[str, str], ...] = ()


@dataclass(frozen=True)
class TherapyResult:
    """Complementary therapy report partitioned into non-interchangeable layers."""

    gene: str
    generated_at: str
    protein_changes: tuple[str, ...]
    variant_therapy_evidence: tuple[TherapyAssociation, ...]
    regulatory_therapy: tuple[TherapyAssociation, ...]
    trial_and_target_discovery: tuple[TherapyAssociation, ...]
    source_statuses: tuple[SourceStatus, ...]
    caveat: str = (
        "Informational result: a gene-drug association, orphan designation, or "
        "trial does not equal an approved therapy or clinical recommendation."
    )

    @property
    def association_count(self) -> int:
        """Role: Count all evidence associations attached to this mutation."""
        return (
            len(self.variant_therapy_evidence)
            + len(self.regulatory_therapy)
            + len(self.trial_and_target_discovery)
        )

    @property
    def associations(self) -> tuple[TherapyAssociation, ...]:
        """Role: Return the flattened evidence associations for this mutation."""
        return (
            self.variant_therapy_evidence
            + self.regulatory_therapy
            + self.trial_and_target_discovery
        )

    def to_dict(self) -> dict[str, Any]:
        """Role: Convert the domain object into a JSON-safe dictionary."""
        payload = asdict(self)
        payload["association_count"] = self.association_count
        return payload


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
    evidence_record_count: int = 0
    evidence_records: tuple[EvidenceRecord, ...] = ()
    source_statuses: tuple[SourceStatus, ...] = ()
    sources: tuple[str, ...] = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        "https://api.gdc.cancer.gov/ssm_occurrences",
    )
    caveat: str = (
        "Candidate matches use canonical protein changes; confirm identity with "
        "coordinates/alleles or transcript HGVS before clinical interpretation."
    )

    def to_dict(self) -> dict[str, Any]:
        """Role: Convert the domain object into a JSON-safe dictionary."""
        return asdict(self)
