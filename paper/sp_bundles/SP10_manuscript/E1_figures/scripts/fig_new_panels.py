#!/usr/bin/env python3
"""SP10 E1: the three figures with no thesis script (steps 27, 18/25/28, 33/33b).

    S_rrna_corrected   27_rrna_corrected/{SED,LIQ}_corrected_profile.tsv
    S_viruses_defence  18_viruses/provirus_host_links.tsv, 28_crispr/defence_systems.tsv,
                       28_crispr/spacers.fasta (array = header field 2)
    S_bgc              33_antismash/mibig_hits.tsv (filled = n_mibig_hits > 0)

Values are read and plotted; nothing is recomputed.  python3 fig_new_panels.py <WS2> <outdir>
"""
import csv, os, sys
from collections import Counter, defaultdict
import etoliko_figstyle as S
import matplotlib.pyplot as plt
import numpy as np

WS2, OUT = sys.argv[1], sys.argv[2]
S.apply()

def rows(rel):
    with open(os.path.join(WS2, rel)) as f:
        return list(csv.DictReader(f, delimiter="\t"))

cw = {r["mag_id"]: r for r in rows("07_mags/MAG_CROSSWALK.tsv")}
master = {r["mag_id"]: r for r in rows("47_master/MASTER.tsv")}
MAGS = sorted(cw)
dom = {m: master[m]["domain"] for m in MAGS}
def short(m): return m.replace("ETO_", "")

# ------------------------------------------------------------- S_rrna_corrected
fig, axes = plt.subplots(1, 2, figsize=(S.FULL_W_IN, 2.6))
for ax, lib, letter in zip(axes, ("SED", "LIQ"), "ab"):
    d = rows("27_rrna_corrected/%s_corrected_profile.tsv" % lib)
    d.sort(key=lambda r: float(r["raw_pct_of_subset"]))
    y = np.arange(len(d)); h = 0.38
    raw = [float(r["raw_pct_of_subset"]) for r in d]; cor = [float(r["corrected_pct_of_subset"]) for r in d]
    ax.barh(y + h / 2, raw, h, color="#bdbdbd", label="raw")
    ax.barh(y - h / 2, cor, h, color=S.DOMAIN["Bacteria"], label="rRNA-copy corrected")
    for i, r in enumerate(d):
        if master[r["source_mag"]]["domain"] == "Archaea":
            ax.patches[len(d) + i].set_color(S.DOMAIN["Archaea"])
    ax.set_yticks(y, ["%s (%s)" % (r["genus"], r["operons"]) for r in d], fontstyle="italic")
    ax.set_xlabel("% of genus-assigned SSU reads in subset")
    S.despine(ax); S.recessive_grid(ax)
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=9, ha="right")
    ax.text(0.98, 0.02, "sediment" if lib == "SED" else "enrichment cultures", transform=ax.transAxes, ha="right", fontsize=7, color=S.INK_2)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color="#bdbdbd", label="raw"), Patch(color=S.DOMAIN["Bacteria"], label="corrected, Bacteria"),
                    Patch(color=S.DOMAIN["Archaea"], label="corrected, Archaea")], loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.08))
fig.tight_layout(rect=(0, 0.06, 1, 1))
S.save(fig, "S_rrna_corrected", OUT)

# ------------------------------------------------------------- S_viruses_defence
prov = Counter(r["mag_id"] for r in rows("18_viruses/provirus_host_links.tsv"))
defn = {r["mag_id"]: int(r["total"]) for r in rows("28_crispr/defence_systems.tsv")}
arrays = defaultdict(set)
for ln in open(os.path.join(WS2, "28_crispr/spacers.fasta")):
    if ln.startswith(">"):
        p = ln[1:].strip().split("|"); arrays[p[0]].add(p[1])
order = sorted(MAGS, key=lambda m: (dom[m] != "Archaea", -(prov[m] + defn.get(m, 0) + len(arrays[m]))))
fig, ax = plt.subplots(figsize=(S.FULL_W_IN, 2.9))
x = np.arange(len(order))
a = [prov[m] for m in order]; b = [defn.get(m, 0) for m in order]; c = [len(arrays[m]) for m in order]
ax.bar(x, a, 0.75, color="#6a3d9a", label="proviruses (geNomad)")
ax.bar(x, b, 0.75, bottom=a, color="#1b7f5f", label="defence systems")
ax.bar(x, c, 0.75, bottom=np.add(a, b), color="#e6ab02", label="CRISPR arrays")
ax.set_xticks(x, [short(m) for m in order], rotation=90)
for t, m in zip(ax.get_xticklabels(), order):
    t.set_color(S.DOMAIN[dom[m]]); t.set_fontweight("bold")
ax.set_ylabel("count per MAG"); S.despine(ax); S.recessive_grid(ax, "y")
ax.legend(ncol=3, loc="upper right")
ax.set_xlim(-0.6, len(order) - 0.4); ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
S.save(fig, "S_viruses_defence", OUT)

# ------------------------------------------------------------- S_bgc
def bclass(p):
    cls = set()
    for t in p.split(";"):
        if "terpene" in t: cls.add("terpene")
        elif t in ("NRPS", "NRPS-like") or "PKS" in t: cls.add("NRPS/PKS")
        elif t in ("ranthipeptide", "RiPP-like", "thioamitides", "triceptide", "sactipeptide", "proteusin",
                   "azole-containing-RiPP", "RRE-containing", "lanthipeptide", "lassopeptide", "thiopeptide") or "peptide" in t:
            cls.add("RiPP")
        else: cls.add("other")
    return cls.pop() if len(cls) == 1 else "hybrid"
CL = ["terpene", "RiPP", "NRPS/PKS", "hybrid", "other"]
COL = dict(zip(CL, ["#1b7f5f", "#6a3d9a", "#b2182b", "#e6ab02", "#7d7c78"]))
cnt = defaultdict(Counter)
for r in rows("33_antismash/mibig_hits.tsv"):
    hit = r["n_mibig_hits"] not in ("", "0", "NA")
    cnt[r["mag_id"]][(bclass(r["products"]), hit)] += 1
order = sorted(MAGS, key=lambda m: (dom[m] != "Archaea", -sum(cnt[m].values())))
fig, ax = plt.subplots(figsize=(S.FULL_W_IN, 2.9))
x = np.arange(len(order)); base = np.zeros(len(order))
for c in CL:
    for hit in (True, False):
        v = np.array([cnt[m][(c, hit)] for m in order])
        ax.bar(x, v, 0.75, bottom=base, color=COL[c] if hit else "white", edgecolor=COL[c], linewidth=0.8)
        base += v
ax.set_xticks(x, [short(m) for m in order], rotation=90)
for t, m in zip(ax.get_xticklabels(), order):
    t.set_color(S.DOMAIN[dom[m]]); t.set_fontweight("bold")
ax.set_ylabel("antiSMASH regions per MAG"); S.despine(ax); S.recessive_grid(ax, "y")
h = [Patch(facecolor=COL[c], edgecolor=COL[c], label=c) for c in CL] + \
    [Patch(facecolor="#7d7c78", label="filled: MIBiG hit"), Patch(facecolor="white", edgecolor="#7d7c78", label="open: no MIBiG hit")]
ax.legend(handles=h, ncol=4, loc="upper right")
ax.set_xlim(-0.6, len(order) - 0.4); ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
S.save(fig, "S_bgc", OUT)
print("regions:", sum(sum(c.values()) for c in cnt.values()), "with MIBiG hit:", sum(v for c in cnt.values() for (k, h_), v in c.items() if h_))
