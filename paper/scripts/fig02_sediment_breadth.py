#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig02_sediment_breadth.py -- το σύνολο in situ, ορισμένο από ένα κατώφλι.

    python3 fig02_sediment_breadth.py [MASTER.tsv] [outdir]

FORM. Magnitude across 29 named entities with long labels -> horizontal ranked bars.
A vertical chart would rotate 29 genus names to 90 degrees; a dot plot would lose the
"how much of the genome was covered" reading a bar length gives directly.

ENCODED ONCE, EACH:
    bar length  sediment breadth
    bar colour  DOMAIN only -- not rank, not set membership
    position    membership of the in-situ ten IS the dashed rule; above it is in.
    one callout MTP4, because the chapter argues about an organism below the cut.

Genus names are drawn ITALIC, as scientific names must be, by placing the labels by
hand rather than as tick labels -- mathtext would choke on GCA-2746605 and on the
underscore in Methanobacterium_C, and would not use the Greek face.

The threshold is the datum, not decoration: the ten were chosen as the top ten by
sediment breadth (step 44), the cut is 0,8936, and the gap beneath it (0,8936 -> 0,6912)
is why the boundary is defensible. The figure exists to show that gap.
"""
import sys, csv
import etoliko_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory as blend

SRC = sys.argv[1] if len(sys.argv) > 1 else "MASTER.tsv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures"
CALLOUT = "ETO_LIQ_MAG_003"          # MTP4

EL = {
    "title":   "Κάλυψη των 29 MAGs στο ίζημα του Αιτωλικού",
    "xlabel":  "Εύρος κάλυψης στο ίζημα",
    "thresh":  "σύνολο in situ\nτα %d πρώτα, όριο %s",
    "callout": "MTP4\nμεθυλοτροφικό\nμεθανογόνο",
    "caption": ("Το εύρος κάλυψης στο ίζημα είναι το κλάσμα κάθε MAG που καλύπτεται από "
                "αναγνώσεις ιζήματος (26_crossmap/cross_abundance.tsv). Η διακεκομμένη "
                "γραμμή σημειώνει το όριο που ορίζει το σύνολο in situ — τα δέκα υψηλότερα, "
                "έως %s. Κανένα γονιδίωμα εκτός του συνόλου δεν το φτάνει, και το κενό από "
                "κάτω (%s έως %s) είναι ο λόγος που το όριο δεν είναι αυθαίρετο. Το χρώμα "
                "δηλώνει επικράτεια· η ένταξη στο σύνολο φαίνεται από τη θέση, όχι από το χρώμα."),
}

# ---------------------------------------------------------------- data
rows = []
with open(SRC, encoding="utf-8") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        try:
            b = float(r.get("sed_breadth") or "nan")
        except ValueError:
            continue
        if b != b:
            continue
        rows.append({
            "mag":    r.get("mag_id") or r.get("working_name", "?"),
            "genus":  (lambda g: "—" if g in ("", "-", "NA", "None") else g)((r.get("genus") or "").strip()),
            "domain": "Archaea" if str(r.get("domain", "")).lower().startswith("arch") else "Bacteria",
            "ten":    str(r.get("in_situ_ten", "")).lower().startswith("y"),
            "b":      b,
        })
if not rows:
    sys.exit("no usable rows in %s" % SRC)
rows.sort(key=lambda r: r["b"])

n_ten = sum(1 for r in rows if r["ten"])
cut = min((r["b"] for r in rows if r["ten"]), default=0.0)
next_below = max((r["b"] for r in rows if not r["ten"]), default=0.0)
print("   %d genomes, %d in the in-situ set, cut at %s" % (len(rows), n_ten, S.num(cut)))

# ---------------------------------------------------------------- figure
S.apply()
H = 0.155 * len(rows) + 2.05
fig = plt.figure(figsize=(7.6, H))
bottom_in, top_in = 1.16, 0.34
ax = fig.add_axes([0.360, bottom_in / H, 0.530, 1 - (bottom_in + top_in) / H])
S.recessive_grid(ax, axis="x")
S.despine(ax, keep=())

for i, r in enumerate(rows):
    ax.barh(i, r["b"], height=0.62, color=S.DOMAIN[r["domain"]], zorder=3, linewidth=0)
    ax.text(r["b"] + 0.012, i, S.num(r["b"]), va="center", ha="left",
            fontsize=S.PT_VALUE, color=S.INK_MUTED, zorder=4)

# ---- labels drawn by hand: MAG id roman, genus italic
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([])
ax.tick_params(axis="y", length=0)
tr = blend(ax.transAxes, ax.transData)
for i, r in enumerate(rows):
    strong = r["ten"]
    # two columns: MAG id flush left, genus flush right against the axis, so the
    # long ones (Methanobacterium_C) grow leftwards into empty space instead of
    # into the neighbouring column -- which is what the first Greek render did.
    ax.text(-0.62, i, r["mag"].replace("ETO_", ""), transform=tr,
            ha="left", va="center", fontsize=S.PT_TICK,
            color=S.INK if strong else S.INK_2)
    ax.text(-0.02, i, r["genus"], transform=tr,
            ha="right", va="center", fontsize=S.PT_TICK, style="italic",
            color=S.INK if strong else S.INK_2)

MARGIN_X, XMAX = 1.09, 1.52

if n_ten:
    ax.axhline(len(rows) - n_ten - 0.5, color=S.INK, lw=1.0, ls=(0, (4, 2)), zorder=5)
    ax.text(MARGIN_X, len(rows) - n_ten - 0.5, EL["thresh"] % (n_ten, S.num(cut)),
            ha="left", va="center", fontsize=S.PT_VALUE, color=S.INK,
            linespacing=1.4, zorder=6)

for i, r in enumerate(rows):
    if r["mag"] == CALLOUT:
        ax.annotate(EL["callout"], xy=(0.76, i), xytext=(MARGIN_X, i),
                    fontsize=S.PT_VALUE, color=S.HILITE, va="center", ha="left",
                    linespacing=1.4,
                    arrowprops=dict(arrowstyle="-", color=S.HILITE, lw=0.8,
                                    shrinkA=2, shrinkB=0),
                    zorder=7)
        break

ax.set_xlim(0, XMAX)
ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xticklabels([S.num(v, 1) for v in (0, 0.2, 0.4, 0.6, 0.8, 1.0)])
ax.set_ylim(-0.8, len(rows) - 0.2)
ax.set_xlabel(EL["xlabel"])
ax.set_title(EL["title"], loc="left", pad=8)

ax.legend(handles=[Patch(facecolor=S.DOMAIN["Bacteria"], label=S.DOMAIN_EL["Bacteria"]),
                   Patch(facecolor=S.DOMAIN["Archaea"],  label=S.DOMAIN_EL["Archaea"])],
          loc="lower right", bbox_to_anchor=(0.62, 0.015), ncols=2,
          handlelength=1.0, handleheight=1.0, borderpad=0.2, columnspacing=1.2)

S.caption(fig, EL["caption"] % (S.num(cut), S.num(cut), S.num(next_below)))
S.save(fig, "fig02_sediment_breadth", OUT)
