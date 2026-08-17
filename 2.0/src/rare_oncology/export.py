"""Result exporters used by the command-line interface."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .models import AnalysisResult


def write_json(result: AnalysisResult, path: str | Path) -> Path:
    """Role: Serialize and write json output."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result.to_dict(), ensure_ascii=False, indent=2) + "\n",
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
        "orphanet_codes",
        "clinvar_accessions",
        "germline_classifications",
        "tcga_projects",
        "cancer_patient_count",
        "cancer_event_count",
        "genomic_match_count",
        "genomic_confirmations",
        "match_level",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in result.shared_mutations:
            writer.writerow(
                {
                    "gene": item.gene,
                    "protein_change": item.protein_change,
                    "diseases": " | ".join(item.diseases),
                    "orphanet_codes": " | ".join(item.orphanet_codes),
                    "clinvar_accessions": " | ".join(item.clinvar_accessions),
                    "germline_classifications": " | ".join(
                        item.germline_classifications
                    ),
                    "tcga_projects": " | ".join(item.tcga_projects),
                    "cancer_patient_count": item.cancer_patient_count,
                    "cancer_event_count": item.cancer_event_count,
                    "genomic_match_count": item.genomic_match_count,
                    "genomic_confirmations": " | ".join(item.genomic_confirmations),
                    "match_level": item.match_level,
                }
            )
    return output
