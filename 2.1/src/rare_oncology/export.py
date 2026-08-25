"""Result exporters used by the command-line interface."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .models import AnalysisResult, EvidenceRecord, TherapyResult


def write_json(result: AnalysisResult, path: str | Path) -> Path:
    """Role: Serialize and write json output."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(
            result.to_dict(), ensure_ascii=False, indent=2, allow_nan=False
        )
        + "\n",
        encoding="utf-8",
    )
    return output


def write_csv(result: AnalysisResult, path: str | Path) -> Path:
    """Role: Serialize and write csv output."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "gene",
        "protein_change",
        "diseases",
        "hereditary_cancer_diseases",
        "other_inherited_diseases",
        "orphanet_codes",
        "clinvar_accessions",
        "germline_classifications",
        "tcga_projects",
        "cancer_patient_count",
        "cancer_event_count",
        "genomic_match_count",
        "genomic_confirmations",
        "match_level",
        "alphamissense_scores",
        "alphamissense_classes",
        "alphamissense_transcripts",
    ]
    alphamissense = _alphamissense_by_protein_change(result)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in result.shared_mutations:
            predictions = alphamissense.get(item.protein_change, ())
            writer.writerow(
                _safe_csv_row(
                    {
                        "gene": item.gene,
                        "protein_change": item.protein_change,
                        "diseases": " | ".join(item.diseases),
                        "hereditary_cancer_diseases": " | ".join(
                            item.hereditary_cancer_diseases
                        ),
                        "other_inherited_diseases": " | ".join(
                            item.other_inherited_diseases
                        ),
                        "orphanet_codes": " | ".join(item.orphanet_codes),
                        "clinvar_accessions": " | ".join(
                            item.clinvar_accessions
                        ),
                        "germline_classifications": " | ".join(
                            item.germline_classifications
                        ),
                        "tcga_projects": " | ".join(item.tcga_projects),
                        "cancer_patient_count": item.cancer_patient_count,
                        "cancer_event_count": item.cancer_event_count,
                        "genomic_match_count": item.genomic_match_count,
                        "genomic_confirmations": " | ".join(
                            item.genomic_confirmations
                        ),
                        "match_level": item.match_level,
                        "alphamissense_scores": " | ".join(
                            sorted(
                                {
                                    f"{record.prediction_score:.6g}"
                                    for record in predictions
                                    if record.prediction_score is not None
                                }
                            )
                        ),
                        "alphamissense_classes": " | ".join(
                            sorted(
                                {
                                    record.prediction_class
                                    for record in predictions
                                    if record.prediction_class
                                }
                            )
                        ),
                        "alphamissense_transcripts": " | ".join(
                            sorted(
                                {
                                    record.transcript
                                    for record in predictions
                                    if record.transcript
                                }
                            )
                        ),
                    }
                )
            )
    return output


def write_evidence_csv(result: AnalysisResult, path: str | Path) -> Path:
    """Role: Serialize and write evidence csv output."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "source",
        "source_record_id",
        "evidence_layer",
        "gene",
        "variant_key",
        "assembly",
        "chromosome",
        "position",
        "reference",
        "alternate",
        "transcript",
        "hgvs_c",
        "protein_change",
        "disease_id",
        "disease_name",
        "mode_of_inheritance",
        "classification",
        "review_level",
        "evidence_level",
        "population_frequency",
        "prediction_score",
        "prediction_class",
        "source_version",
        "retrieved_at",
        "license",
        "details",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in result.evidence_records:
            row = {field: getattr(record, field) for field in fieldnames[:-1]}
            row["details"] = json.dumps(
                dict(record.details), ensure_ascii=False, sort_keys=True
            )
            writer.writerow(_safe_csv_row(row))
    return output


def _alphamissense_by_protein_change(
    result: AnalysisResult,
) -> dict[str, tuple[EvidenceRecord, ...]]:
    """Role: Compute alphamissense by protein change for the export workflow."""
    grouped: dict[str, list[EvidenceRecord]] = {}
    for record in result.evidence_records:
        if record.source != "alphamissense" or not record.protein_change:
            continue
        grouped.setdefault(record.protein_change, []).append(record)
    return {change: tuple(records) for change, records in grouped.items()}


def write_therapy_json(result: TherapyResult, path: str | Path) -> Path:
    """Role: Serialize and write therapy json output."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(
            result.to_dict(), ensure_ascii=False, indent=2, allow_nan=False
        )
        + "\n",
        encoding="utf-8",
    )
    return output


def write_therapy_csv(result: TherapyResult, path: str | Path) -> Path:
    """Role: Serialize and write therapy csv output."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "layer",
        "source",
        "source_record_id",
        "gene",
        "association_scope",
        "protein_change",
        "molecular_profile",
        "therapy_name",
        "therapy_components",
        "therapy_type",
        "chembl_ids",
        "canonical_smiles",
        "therapy_id",
        "therapy_status",
        "effect",
        "disease_ids",
        "disease_names",
        "evidence_type",
        "evidence_level",
        "evidence_direction",
        "evidence_rating",
        "jurisdiction",
        "trial_id",
        "trial_phase",
        "trial_status",
        "trial_role",
        "citation_id",
        "description",
        "source_url",
        "source_version",
        "retrieved_at",
        "license",
        "details",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in result.associations:
            row = {field: getattr(record, field) for field in fieldnames[:-1]}
            for field in (
                "therapy_components",
                "chembl_ids",
                "disease_ids",
                "disease_names",
            ):
                row[field] = " | ".join(getattr(record, field))
            row["details"] = json.dumps(
                dict(record.details), ensure_ascii=False, sort_keys=True
            )
            writer.writerow(_safe_csv_row(row))
    return output


def _safe_csv_row(row: dict[str, object]) -> dict[str, object]:
    """Role: Sanitize csv row for safe serialization."""
    return {key: _safe_csv_value(value) for key, value in row.items()}


def _safe_csv_value(value: object) -> object:
    """Role: Sanitize csv value for safe serialization."""
    if isinstance(value, str):
        trimmed = value.lstrip(" \t\r\n")
        if value.startswith(("\t", "\r", "\n")) or trimmed.startswith(
            ("=", "+", "-", "@")
        ):
            return "'" + value
    return value
