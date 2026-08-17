"""Orchestrate and aggregate rare-disease/oncology mutation matches."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import replace
from datetime import datetime, timezone

from .clients import GdcClient, NcbiClinVarClient
from .models import (
    AnalysisResult,
    CancerObservation,
    MutationPlotRecord,
    RareVariant,
    SharedMutation,
)
from .normalization import clinvar_spdi_key, gdc_genomic_change_key, protein_position


DEFAULT_ALLOWED_CLASSIFICATIONS = (
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
)


def find_shared_mutations(
    gene: str,
    *,
    min_cancer_patients: int = 1,
    program: str | None = "TCGA",
    projects: Iterable[str] | None = None,
    max_clinvar_records: int | None = None,
    exclude_conditions: Iterable[str] | None = None,
    allowed_germline_classifications: Iterable[str] | None = (
        DEFAULT_ALLOWED_CLASSIFICATIONS
    ),
    clinvar_client: NcbiClinVarClient | None = None,
    gdc_client: GdcClient | None = None,
) -> AnalysisResult:
    """Find candidate protein-level matches between rare disease and oncology.

    Rare-disease evidence is defined objectively as a ClinVar germline trait with
    an Orphanet cross-reference. Oncology evidence comes from GDC somatic mutation
    occurrences. A match is candidate-level until genomic identity is confirmed.
    """

    normalized_gene = gene.strip().upper()
    selected_projects = tuple(sorted(set(projects or ())))
    selected_exclusions = tuple(sorted(set(exclude_conditions or ())))
    selected_classifications = (
        tuple(sorted(set(allowed_germline_classifications)))
        if allowed_germline_classifications is not None
        else None
    )
    if min_cancer_patients < 1:
        raise ValueError("min_cancer_patients must be at least 1")

    rare_variants = (clinvar_client or NcbiClinVarClient()).fetch_rare_variants(
        normalized_gene, max_records=max_clinvar_records
    )
    excluded = {condition.casefold() for condition in selected_exclusions}
    if excluded:
        filtered_variants: list[RareVariant] = []
        for variant in rare_variants:
            pairs = tuple(
                pair
                for pair in variant.condition_orphanet_pairs
                if pair[0].casefold() not in excluded
            )
            if not pairs:
                continue
            filtered_variants.append(
                replace(
                    variant,
                    diseases=tuple(sorted({disease for disease, _code in pairs})),
                    orphanet_codes=tuple(sorted({code for _disease, code in pairs})),
                    condition_orphanet_pairs=pairs,
                )
            )
        rare_variants = filtered_variants
    if selected_classifications is not None:
        allowed = {
            classification.casefold() for classification in selected_classifications
        }
        rare_variants = [
            variant
            for variant in rare_variants
            if variant.germline_classification.casefold() in allowed
        ]
    gdc = gdc_client or GdcClient()
    cancer_observations = gdc.fetch_cancer_observations(
        normalized_gene, program=program, projects=selected_projects
    )
    shared = _aggregate_matches(
        normalized_gene,
        rare_variants,
        cancer_observations,
        min_cancer_patients=min_cancer_patients,
    )
    plot_data = _build_plot_data(
        normalized_gene,
        rare_variants,
        cancer_observations,
        shared,
        min_cancer_patients=min_cancer_patients,
    )
    return AnalysisResult(
        gene=normalized_gene,
        generated_at=datetime.now(timezone.utc).isoformat(),
        gdc_program=program,
        gdc_projects=selected_projects,
        min_cancer_patients=min_cancer_patients,
        excluded_conditions=selected_exclusions,
        clinvar_record_limit=max_clinvar_records,
        allowed_germline_classifications=selected_classifications,
        gdc_data_release=(
            str(gdc.release_metadata.get("data_release"))
            if getattr(gdc, "release_metadata", {}).get("data_release")
            else None
        ),
        rare_variant_count=len(rare_variants),
        cancer_observation_count=len(cancer_observations),
        shared_mutation_count=len(shared),
        shared_mutations=tuple(shared),
        mutation_plot_data=tuple(plot_data),
    )


def _aggregate_matches(
    gene: str,
    rare_variants: Iterable[RareVariant],
    cancer_observations: Iterable[CancerObservation],
    *,
    min_cancer_patients: int,
) -> list[SharedMutation]:
    """Role: Compute aggregate matches for the analysis workflow."""
    rare_by_change: dict[str, list[RareVariant]] = defaultdict(list)
    cancer_by_change: dict[str, list[CancerObservation]] = defaultdict(list)
    for variant in rare_variants:
        rare_by_change[variant.protein_change].append(variant)
    for observation in cancer_observations:
        cancer_by_change[observation.protein_change].append(observation)

    matches: list[SharedMutation] = []
    for change in sorted(rare_by_change.keys() & cancer_by_change.keys()):
        rare = rare_by_change[change]
        cancer = cancer_by_change[change]
        patient_keys = {
            observation.case_id or observation.patient_id
            for observation in cancer
            if observation.case_id or observation.patient_id
        }
        if len(patient_keys) < min_cancer_patients:
            continue

        rare_keys = [clinvar_spdi_key(variant.canonical_spdi) for variant in rare]
        cancer_keys = [
            gdc_genomic_change_key(observation.genomic_dna_change)
            for observation in cancer
        ]
        genomic_confirmations = {
            (
                f"{variant.clinvar_accession}|{observation.ssm_id}|"
                f"{key[0]}:{key[1]}{key[2]}>{key[3]}"
            )
            for variant, rare_key in zip(rare, rare_keys, strict=True)
            for observation, cancer_key in zip(cancer, cancer_keys, strict=True)
            if rare_key is not None
            and cancer_key is not None
            and rare_key == cancer_key
            for key in (rare_key,)
        }
        confirmed_rare = {
            index
            for index, rare_key in enumerate(rare_keys)
            if rare_key is not None and rare_key in set(cancer_keys)
        }
        confirmed_cancer = {
            index
            for index, cancer_key in enumerate(cancer_keys)
            if cancer_key is not None and cancer_key in set(rare_keys)
        }
        all_records_confirmed = (
            bool(genomic_confirmations)
            and len(confirmed_rare) == len(rare)
            and len(confirmed_cancer) == len(cancer)
        )

        matches.append(
            SharedMutation(
                gene=gene,
                protein_change=change,
                diseases=_flatten(variant.diseases for variant in rare),
                orphanet_codes=_flatten(variant.orphanet_codes for variant in rare),
                clinvar_ids=tuple(sorted({variant.clinvar_id for variant in rare})),
                clinvar_accessions=tuple(
                    sorted(
                        {
                            variant.clinvar_accession
                            for variant in rare
                            if variant.clinvar_accession
                        }
                    )
                ),
                germline_classifications=tuple(
                    sorted(
                        {
                            variant.germline_classification
                            for variant in rare
                            if variant.germline_classification
                        }
                    )
                ),
                canonical_spdis=tuple(
                    sorted(
                        {
                            variant.canonical_spdi
                            for variant in rare
                            if variant.canonical_spdi
                        }
                    )
                ),
                tcga_projects=tuple(
                    sorted(
                        {
                            observation.project_id
                            for observation in cancer
                            if observation.project_id
                        }
                    )
                ),
                cancer_patient_count=len(patient_keys),
                cancer_event_count=len({observation.ssm_id for observation in cancer}),
                genomic_match_count=len(genomic_confirmations),
                genomic_confirmations=tuple(sorted(genomic_confirmations)),
                ssm_ids=tuple(sorted({observation.ssm_id for observation in cancer})),
                match_level=(
                    "confirmed_genomic_and_protein"
                    if all_records_confirmed
                    else (
                        "mixed_genomic_and_candidate_protein"
                        if genomic_confirmations
                        else "candidate_protein_match"
                    )
                ),
            )
        )
    return sorted(
        matches,
        key=lambda item: (-item.cancer_patient_count, item.protein_change),
    )


def _flatten(groups: Iterable[Iterable[str]]) -> tuple[str, ...]:
    """Role: Compute flatten for the analysis workflow."""
    return tuple(sorted({value for group in groups for value in group if value}))


def _build_plot_data(
    gene: str,
    rare_variants: Iterable[RareVariant],
    cancer_observations: Iterable[CancerObservation],
    shared_mutations: Iterable[SharedMutation],
    *,
    min_cancer_patients: int,
) -> list[MutationPlotRecord]:
    """Role: Compute build plot data for the analysis workflow."""
    rare_by_change: dict[str, list[RareVariant]] = defaultdict(list)
    cancer_by_change: dict[str, list[CancerObservation]] = defaultdict(list)
    for variant in rare_variants:
        rare_by_change[variant.protein_change].append(variant)
    for observation in cancer_observations:
        cancer_by_change[observation.protein_change].append(observation)
    shared_by_change = {
        mutation.protein_change: mutation for mutation in shared_mutations
    }

    records: list[MutationPlotRecord] = []
    for change in sorted(rare_by_change.keys() | cancer_by_change.keys()):
        rare = rare_by_change.get(change, [])
        cancer = cancer_by_change.get(change, [])
        patients = {
            observation.case_id or observation.patient_id
            for observation in cancer
            if observation.case_id or observation.patient_id
        }
        qualifies_as_oncology = len(patients) >= min_cancer_patients
        if not rare and not qualifies_as_oncology:
            continue
        project_patients: dict[str, set[str]] = defaultdict(set)
        if qualifies_as_oncology:
            for observation in cancer:
                patient = observation.case_id or observation.patient_id
                if observation.project_id and patient:
                    project_patients[observation.project_id].add(patient)

        records.append(
            MutationPlotRecord(
                gene=gene,
                protein_change=change,
                amino_acid_position=protein_position(change),
                in_rare_disease=bool(rare),
                in_oncology=qualifies_as_oncology,
                diseases=_flatten(variant.diseases for variant in rare),
                orphanet_codes=_flatten(variant.orphanet_codes for variant in rare),
                germline_classifications=_flatten(
                    (variant.germline_classification,) for variant in rare
                ),
                cancer_patient_count=(len(patients) if qualifies_as_oncology else 0),
                cancer_event_count=(
                    len({observation.ssm_id for observation in cancer})
                    if qualifies_as_oncology
                    else 0
                ),
                project_patient_counts=tuple(
                    sorted(
                        (project, len(project_cases))
                        for project, project_cases in project_patients.items()
                    )
                ),
                match_level=(
                    shared_by_change[change].match_level
                    if change in shared_by_change
                    else None
                ),
            )
        )
    return records
