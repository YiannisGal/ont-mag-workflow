"""Phylum relative abundance from a SingleM taxonomic profile.

Filled coverage of a phylum = coverage assigned to the phylum itself plus all
its descendants. Relative abundance = filled phylum coverage / sum of filled
phylum coverage (reads not assigned to any phylum are excluded).
"""
import collections

import pandas as pd

smk = snakemake  # noqa: F821

prof = pd.read_csv(smk.input[0], sep="\t")
filled = collections.Counter()
for tax, cov in zip(prof["taxonomy"], prof["coverage"]):
    levels = [x.strip() for x in str(tax).split(";")]
    phy = [x for x in levels if x.startswith("p__")]
    if phy:
        domain = next((x for x in levels if x.startswith("d__")), "")
        filled[(domain, phy[0])] += float(cov)

total = sum(filled.values()) or 1.0
rows = [{"domain": d, "phylum": p, "filled_coverage": c, "relative_abundance_pct": 100 * c / total}
        for (d, p), c in filled.items()]
pd.DataFrame(rows).sort_values("filled_coverage", ascending=False).to_csv(smk.output[0], sep="\t", index=False)
