"""Competitive best-hit assignment of reads to a labelled protein panel.

Each read is assigned to the gene label of its single best hit (highest
bitscore). Sibling families must be in the panel: if, e.g., mtmB is searched
without the other methyltransferases, reads from siblings will be mis-assigned.
"""
import pandas as pd

smk = snakemake  # noqa: F821

labels = pd.read_csv(smk.input.labels, sep="\t", header=None, names=["protein", "gene"], dtype=str)
lab = dict(zip(labels["protein"], labels["gene"]))

cols = ["read", "protein", "pident", "length", "evalue", "bitscore"]
try:
    h = pd.read_csv(smk.input.hits, sep="\t", header=None, names=cols)
except pd.errors.EmptyDataError:
    h = pd.DataFrame(columns=cols)
h = h[h["bitscore"] >= smk.params.min_bits]
best = h.sort_values("bitscore", ascending=False).drop_duplicates("read")
best["gene"] = best["protein"].map(lab).fillna("unlabelled")

out = best.groupby("gene").agg(reads=("read", "size"), median_pident=("pident", "median")).reset_index()
out.sort_values("reads", ascending=False).to_csv(smk.output[0], sep="\t", index=False)
best.to_csv(str(smk.output[0]).replace(".counts.tsv", ".best_hits.tsv"), sep="\t", index=False)
