"""Public API for the rare disease/oncology mutation cross-checker."""

from .analysis import find_shared_mutations
from .clients import GdcClient, NcbiClinVarClient
from .models import (
    AnalysisResult,
    CancerObservation,
    MutationPlotRecord,
    RareVariant,
    SharedMutation,
)
from .normalization import normalize_protein_change
from .plots import generate_gene_plots
from .protein_annotations import (
    ProteinAnnotation,
    ProteinFeature,
    UniProtProteinClient,
)

__all__ = [
    "AnalysisResult",
    "CancerObservation",
    "GdcClient",
    "MutationPlotRecord",
    "NcbiClinVarClient",
    "ProteinAnnotation",
    "ProteinFeature",
    "RareVariant",
    "SharedMutation",
    "UniProtProteinClient",
    "find_shared_mutations",
    "generate_gene_plots",
    "normalize_protein_change",
]
