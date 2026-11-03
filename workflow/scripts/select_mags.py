"""Combine CheckM2, GUNC, barrnap and tRNAscan-SE results, classify bins by
MIMAG, select MAGs and give them stable published IDs.

MIMAG (Bowers et al. 2017):
  high-quality   : completeness > 90, contamination < 5, 23S+16S+5S present, >= 18 tRNA types
  medium-quality : completeness >= 50, contamination < 10
  low-quality    : everything else
"""
import glob
import os
import re
import shutil

import pandas as pd

smk = snakemake  # noqa: F821  (injected by Snakemake)

bins_dir = smk.input.bins
out_dir = smk.output.mags
os.makedirs(out_dir, exist_ok=True)

# ---- CheckM2 ---------------------------------------------------------
ck = pd.read_csv(smk.input.checkm2, sep="\t").rename(columns={"Name": "bin"})
keep = [c for c in ["bin", "Completeness", "Contamination", "Genome_Size", "Contig_N50",
                    "GC_Content", "Coding_Density", "Average_Gene_Length",
                    "Total_Coding_Sequences"] if c in ck.columns]
df = ck[keep].copy()

# ---- GUNC ------------------------------------------------------------
gunc_files = glob.glob(os.path.join(smk.input.gunc, "*maxCSS_level.tsv"))
if gunc_files:
    g = pd.read_csv(gunc_files[0], sep="\t")
    g = g.rename(columns={"genome": "bin"})[["bin", "clade_separation_score", "pass.GUNC"]]
    g["pass.GUNC"] = g["pass.GUNC"].astype(str).str.lower().isin(["true", "1"])
    df = df.merge(g, on="bin", how="left")
else:
    df["clade_separation_score"] = float("nan")
    df["pass.GUNC"] = pd.NA

# ---- rRNA / tRNA -----------------------------------------------------
STANDARD_AA = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
               "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}


def rrna(b):
    found = set()
    for k in ("bac", "arc"):
        p = os.path.join(smk.input.rrna_trna, f"{b}.{k}.gff")
        if not os.path.exists(p):
            continue
        for line in open(p):
            m = re.search(r"Name=(5S|16S|23S)_rRNA", line)
            if m:
                found.add(m.group(1))
    return found


def trna_types(b):
    p = os.path.join(smk.input.rrna_trna, f"{b}.trna.txt")
    types = set()
    if not os.path.exists(p):
        return 0
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or not f[1].strip().isdigit():
            continue  # header lines
        note = f[-1].lower() if len(f) > 9 else ""
        if "pseudo" in note:
            continue
        aa = re.sub(r"\d+$", "", f[4].strip())
        aa = {"iMet": "Met", "fMet": "Met", "Ile2": "Ile"}.get(aa, aa)
        if aa in STANDARD_AA:
            types.add(aa)
    return len(types)


rr = df["bin"].apply(rrna)
for s in ("5S", "16S", "23S"):
    df[f"has_{s}"] = rr.apply(lambda x, s=s: s in x)
df["n_trna_types"] = df["bin"].apply(trna_types)

# ---- MIMAG -----------------------------------------------------------


def mimag(r):
    if (r.Completeness > 90 and r.Contamination < 5 and r.has_5S and r.has_16S
            and r.has_23S and r.n_trna_types >= 18):
        return "high-quality"
    if r.Completeness >= 50 and r.Contamination < 10:
        return "medium-quality"
    return "low-quality"


df["mimag"] = df.apply(mimag, axis=1)
df["sample"] = df["bin"].str.split("__").str[0]
df["winning_binner"] = df["bin"].str.split("__").str[1].str.split("_").str[0]
df["quality_score"] = df["Completeness"] - 5 * df["Contamination"]

gunc_ok = df["pass.GUNC"].fillna(True).astype(bool) if smk.params.gunc_required else True
df["selected"] = ((df["Completeness"] >= smk.params.min_comp)
                  & (df["Contamination"] <= smk.params.max_cont)
                  & gunc_ok)

# ---- naming: <PREFIX>_<SAMPLE>_MAG_<NNN>, best first ------------------
rank = {"high-quality": 0, "medium-quality": 1, "low-quality": 2}
df["_r"] = df["mimag"].map(rank)
df["mag_id"] = ""
for sample, sub in df[df["selected"]].groupby("sample"):
    sub = sub.sort_values(["_r", "quality_score"], ascending=[True, False])
    for i, idx in enumerate(sub.index, 1):
        df.loc[idx, "mag_id"] = f"{smk.params.prefix}_{sample}_MAG_{i:03d}"

for _, r in df[df["selected"]].iterrows():
    shutil.copy(os.path.join(bins_dir, f"{r['bin']}.fa"), os.path.join(out_dir, f"{r['mag_id']}.fa"))

df = df.drop(columns="_r").rename(columns={"bin": "working_name"})
cols = ["mag_id", "working_name", "sample", "winning_binner", "selected", "mimag"] + \
       [c for c in df.columns if c not in {"mag_id", "working_name", "sample", "winning_binner", "selected", "mimag"}]
df = df[cols].sort_values(["selected", "sample", "mag_id"], ascending=[False, True, True])
df.to_csv(smk.output.summary, sep="\t", index=False)
df[df["selected"]][["mag_id", "working_name", "sample", "winning_binner", "Completeness",
                    "Contamination", "pass.GUNC", "mimag"]].to_csv(smk.output.crosswalk, sep="\t", index=False)
print(f"{int(df['selected'].sum())} MAGs selected out of {len(df)} bins")
