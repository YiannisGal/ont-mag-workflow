"""Per-genome detection call for one read set.

For every primary alignment covering >= min_aln % of the read, compute
BLAST-like identity = 1 - NM / (M + I + D) and build a per-genome histogram
(1 % bins). The mode of that histogram is compared with breadth at >= 90 %
identity (CoverM):

  same-species population : mode >= mode_thr  and breadth90 >= same_breadth
  related population      : breadth90 >= related_breadth (lower mode)
  not detected            : otherwise
"""
import collections
import os

import pandas as pd
import pysam

smk = snakemake  # noqa: F821

min_aln = smk.params.min_aln / 100.0
hist = collections.defaultdict(collections.Counter)

with pysam.AlignmentFile(smk.input.bam, "rb") as bam:
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        qlen = a.infer_read_length() or a.query_length
        if not qlen or a.query_alignment_length / qlen < min_aln:
            continue
        cols = sum(n for op, n in a.cigartuples if op in (0, 1, 2, 7, 8))
        if cols == 0 or not a.has_tag("NM"):
            continue
        ident = 1 - a.get_tag("NM") / cols
        genome = a.reference_name.split("~", 1)[0]
        hist[genome][int(ident * 100)] += 1


def read_coverm(path, suffix):
    d = pd.read_csv(path, sep="\t")
    d.columns = ["genome"] + [c.split(" ", 1)[1].strip() + suffix if " " in c else c + suffix
                              for c in d.columns[1:]]
    return d


df = read_coverm(smk.input.base, "")
for ident, path in zip(smk.params.identities, smk.input.ids):
    df = df.merge(read_coverm(path, f" (id>={ident})"), on="genome", how="left")
df = df[df["genome"] != "unmapped"]

df["n_reads_identity_hist"] = df["genome"].map(lambda g: sum(hist[g].values()))
df["identity_mode"] = df["genome"].map(lambda g: hist[g].most_common(1)[0][0] if hist[g] else float("nan"))

b90_col = [c for c in df.columns if c.startswith("Covered Fraction") and "(id>=90)" in c]
if not b90_col:
    raise SystemExit("detection.identities must include 90 for the detection call")
b90 = df[b90_col[0]]


def call(mode, b):
    if b >= smk.params.same_breadth and mode >= smk.params.mode_thr:
        return "same-species population"
    if b >= smk.params.related_breadth:
        return "related population"
    return "not detected"


df["detection_call"] = [call(m if m == m else 0, b) for m, b in zip(df["identity_mode"], b90)]
df.to_csv(smk.output[0], sep="\t", index=False)
