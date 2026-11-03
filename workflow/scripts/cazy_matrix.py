"""MAG x CAZy family matrix from dbCAN overview.txt.

A protein counts as a CAZyme when >= 2 of the 3 dbCAN tools agree
(#ofTools >= 2). Families are taken from the HMMER column (e.g. GH5_2(10-300)
-> GH5); each protein counts once per family.
"""
import os
import re

import pandas as pd

smk = snakemake  # noqa: F821

rows = []
for path in smk.input:
    mag = os.path.basename(os.path.dirname(path))
    o = pd.read_csv(path, sep="\t", dtype=str)
    o.columns = [c.strip() for c in o.columns]
    ntools = next(c for c in o.columns if c.lower().replace(" ", "") in ("#oftools", "#oftools"))
    hmm = next(c for c in o.columns if c.upper().startswith("HMMER"))
    gid = o.columns[0]
    for _, r in o.iterrows():
        try:
            if int(r[ntools]) < 2:
                continue
        except (TypeError, ValueError):
            continue
        fams = set()
        for tok in str(r[hmm]).split("+"):
            m = re.match(r"([A-Z]+\d+)", tok.strip())
            if m:
                fams.add(m.group(1))
        for fam in fams:
            rows.append((mag, r[gid], fam))

df = pd.DataFrame(rows, columns=["mag", "protein", "family"])
if df.empty:
    pd.DataFrame().to_csv(smk.output[0], sep="\t")
else:
    df.groupby(["family", "mag"]).size().unstack(fill_value=0).to_csv(smk.output[0], sep="\t")
