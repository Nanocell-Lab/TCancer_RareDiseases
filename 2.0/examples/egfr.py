"""Live EGFR example. Run from the 2.0 directory."""

from rare_oncology import find_shared_mutations, generate_gene_plots
from rare_oncology.export import write_csv, write_json


result = find_shared_mutations("EGFR")
write_json(result, "results/egfr_report.json")
write_csv(result, "results/egfr_shared.csv")
generate_gene_plots(result, "results/plots")
print(result.to_dict())
