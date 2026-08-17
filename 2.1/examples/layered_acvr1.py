"""ACVR1 example with optional evidence layers."""

from rare_oncology import (
    GnomadEvidenceProvider,
    OrphadataXmlEvidenceProvider,
    PanelAppEvidenceProvider,
    find_shared_mutations,
)


result = find_shared_mutations(
    "ACVR1",
    evidence_providers=[
        GnomadEvidenceProvider(max_variants=25),
        PanelAppEvidenceProvider(),
        OrphadataXmlEvidenceProvider("../data/Orphanet/en_product6.xml"),
    ],
)

print(result.gene, result.shared_mutation_count, result.evidence_record_count)
for status in result.source_statuses:
    print(status.source, status.status, status.record_count, status.evidence_layer)
