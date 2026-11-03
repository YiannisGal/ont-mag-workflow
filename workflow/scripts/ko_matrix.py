"""MAG x KO count matrix from KofamScan significant hits (curated thresholds)."""
import os

import pandas as pd

smk = snakemake  # noqa: F821

rows = []
for path in smk.input:
    mag = os.path.basename(path).replace(".significant.tsv", "")
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[0] == "*":
            rows.append((mag, f[1], f[2]))  # gene, KO

df = pd.DataFrame(rows, columns=["mag", "gene", "KO"]).drop_duplicates()
if df.empty:
    pd.DataFrame().to_csv(smk.output[0], sep="\t")
else:
    m = df.groupby(["KO", "mag"]).size().unstack(fill_value=0)
    m.to_csv(smk.output[0], sep="\t")
