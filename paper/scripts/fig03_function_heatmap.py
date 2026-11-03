#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig03_function_heatmap.py -- τι μπορεί να κάνει κάθε γονιδίωμα.

    python3 fig03_function_heatmap.py [MASTER.tsv] [outdir]
    python3 fig03_function_heatmap.py [MASTER.tsv] --values   # inspect first

FORM. One quantity measured the same way across 29 genomes and a dozen functions
is a matrix, and a matrix reads as a heatmap. The alternative -- twelve small bar
charts -- would make the reader compare across panels to answer the only question
this figure exists for: which genomes can do which things, and do the archaea
differ from the bacteria.

THE ONE PLACE THIS CHAPTER LETS COLOUR CARRY A VALUE. Everywhere else colour is
domain and nothing else. A heatmap cannot work that way, so: the CELLS carry the
quantity on a single sequential ramp with a key, every cell also prints its own
value, and DOMAIN moves to a strip down the left edge. Nothing is identified by
colour alone, which is the actual rule.

READ THE VALUES BEFORE BELIEVING THE PICTURE. The module columns in MASTER.tsv
may be "3/5", a count, or a yes/no, and the three normalise differently. The
script detects the format per column, PRINTS what it decided, and refuses to
guess silently. Run with --values first if anything looks off.
"""
import sys, csv, os, re
import etoliko_figstyle as S
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch, Rectangle

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SRC = args[0] if args else "MASTER.tsv"
OUT = args[1] if len(args) > 1 else "figures"
SHOW_VALUES = "--values" in sys.argv
FIG_W_MM = 160.0

# column -> label on the panel. Edit here to add, drop or reorder.
PANEL = [
    ("mcr",                "mcr"),
    ("methylamine_MT",     "μεθυλαμίνη"),
    ("methanol_MT",        "μεθανόλη"),
    ("fwd_fmd",            "fwd/fmd"),
    ("hydrogenase",        "υδρογονάση"),
    ("fdh",                "fdh"),
    ("ack_pta",            "ack/pta"),
    ("dsr_apr_sat",        "dsr/apr/sat"),
    ("nif",                "nif"),
    ("B12_biosynth",       "B$_{12}$ σύνθεση"),
    ("B12_dependent_wide", "B$_{12}$ εξάρτηση"),
    ("B12_transport",      "B$_{12}$ μεταφορά"),
]

EL = {
    "title": "Λειτουργικό δυναμικό των 29 MAGs",
    "cbar":  "αριθμός γονιδίων, ως προς το μέγιστο της στήλης",
    "caption": ("Κάθε κελί δείχνει τον αριθμό των γονιδίων του αντίστοιχου συνόλου σε κάθε "
                "γονιδίωμα· ο αριθμός τυπώνεται μέσα στο κελί και "
                "το χρώμα τον αποδίδει κανονικοποιημένο ως προς το μέγιστο της κάθε στήλης, "
                "ώστε οι στήλες να είναι συγκρίσιμες μεταξύ τους — δεν είναι πληρότητα "
                "μονοπατιού. Τα άδεια κελιά σημαίνουν μηδέν γονίδια. Η λωρίδα αριστερά "
                "δηλώνει επικράτεια· τα "
                "γονιδιώματα είναι ταξινομημένα κατά επικράτεια και, μέσα σε κάθε ομάδα, "
                "κατά φθίνον εύρος κάλυψης στο ίζημα, ώστε η σειρά να είναι η ίδια με το "
                "σχήμα 2. Οι απουσίες είναι ερμηνεύσιμες μόνο σε γονιδιώματα επαρκούς "
                "πληρότητας· δείτε το σχήμα 1."),
}

# ---------------------------------------------------------------- read
if not os.path.exists(SRC):
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, os.path.basename(SRC)),
                 os.path.join(here, "..", os.path.basename(SRC)),
                 os.path.join(here, "..", "47_master", os.path.basename(SRC))):
        if os.path.exists(cand):
            SRC = cand; print("   found %s" % os.path.normpath(cand)); break
    else:
        sys.exit("fig03: cannot find %s -- pass the path as the first argument" % SRC)

with open(SRC, encoding="utf-8") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
if not rows:
    sys.exit("fig03: %s is empty" % SRC)
have = set(rows[0].keys())

missing = [c for c, _ in PANEL if c not in have]
if missing:
    print("   !! not in %s, dropped from the panel: %s" % (SRC, ", ".join(missing)))
PANEL = [(c, l) for c, l in PANEL if c in have]
if not PANEL:
    sys.exit("fig03: none of the panel columns exist. Columns present:\n   %s"
             % "\n   ".join(sorted(have)))

if SHOW_VALUES:
    print("distinct values per panel column:")
    for c, _ in PANEL:
        vals = sorted({(r.get(c) or "").strip() for r in rows})
        print("   %-20s %s" % (c, ", ".join(v or "<empty>" for v in vals[:14])))
    sys.exit(0)

# ---------------------------------------------------------------- parse
# CLASSIFY THE COLUMN BEFORE THE CELL. The first version judged each cell on its
# own, so in a count column like mcr (values 0, 3, 5) every "0" looked boolean
# and every "3" looked numeric -- and since most genomes carry no mcr, "boolean"
# won the vote and nine columns were mis-read. A column is one kind of thing;
# decide what it is from all of its values, then read the cells that way.
FRAC = re.compile(r"^\s*(\d+(?:[.,]\d+)?)\s*/\s*(\d+(?:[.,]\d+)?)\s*$")
YES = {"y", "yes", "true", "t", "+", "present"}
NO = {"n", "no", "false", "f", "-", "absent"}


def raw(v):
    v = (v or "").strip()
    return None if v == "" or v.lower() in ("na", "nan", "none", ".") else v


def as_num(v):
    try: return float(v.replace(",", "."))
    except ValueError: return None


grid, kinds = [], {}
for c, _ in PANEL:
    vals = [raw(r.get(c)) for r in rows]
    nb = [v for v in vals if v is not None]
    nums = [as_num(v) for v in nb]

    if nb and all(FRAC.match(v) for v in nb):
        kind = "frac"
    elif nb and all(n is not None for n in nums):
        kind = "bool01" if set(nums) <= {0.0, 1.0} else "count"
    elif nb and all(v.lower() in YES | NO for v in nb):
        kind = "booltext"
    else:
        kind = "text"

    mx = max([n for n in nums if n is not None] or [1.0])
    col = []
    for v in vals:
        if v is None:
            col.append((None, "")); continue
        if kind == "frac":
            a, b = (float(x.replace(",", ".")) for x in FRAC.match(v).groups())
            col.append(((a / b if b else None), v.replace(" ", "")))
        elif kind == "count":
            n = as_num(v) or 0.0
            # 0 is drawn as an empty cell and carries no number: 29 x 12 zeros is
            # noise, and absence reads better as emptiness than as a printed nought
            col.append((n / mx if mx else 0.0, ("%g" % n) if n > 0 else ""))
        elif kind == "bool01":
            n = as_num(v) or 0.0
            col.append((n, "+" if n else ""))
        elif kind == "booltext":
            y = v.lower() in YES
            col.append((1.0 if y else 0.0, "+" if y else ""))
        else:
            col.append((None, v[:5]))
    kinds[c] = kind + ("" if kind != "count" else "  (max %g)" % mx)
    grid.append(col)

print("   %d genomes x %d functions" % (len(rows), len(PANEL)))
print("   column formats detected:")
for c, _ in PANEL:
    print("      %-20s %s" % (c, kinds[c]))

# ---------------------------------------------------------------- order
def f(r, key, default=0.0):
    try: return float((r.get(key) or "").replace(",", "."))
    except (TypeError, ValueError): return default

def is_arc(r): return str(r.get("domain", "")).lower().startswith("arch")

order = sorted(range(len(rows)),
               key=lambda i: (0 if is_arc(rows[i]) else 1,
                              -f(rows[i], "sed_breadth", -1),
                              -f(rows[i], "completeness")))

# ---------------------------------------------------------------- figure
S.apply()
W = FIG_W_MM / 25.4
RH = 0.150                                  # inches per genome row
H = RH * len(rows) + 2.15
fig = plt.figure(figsize=(W, H))
left, right = 0.345, 0.935
bot_in, top_in = 1.34, 1.00
ax = fig.add_axes([left, bot_in / H, right - left, 1 - (bot_in + top_in) / H])

ramp = mcolors.LinearSegmentedColormap.from_list(
    "eto", ["#eef2f5", "#9cc0d8", S.DOMAIN["Bacteria"], "#0c3f61"])

n, m = len(rows), len(PANEL)
for j in range(m):
    for rank, i in enumerate(order):
        val, txt = grid[j][i]
        y = n - 1 - rank
        if val is None:
            ax.add_patch(Rectangle((j, y), 1, 1, facecolor=S.SURFACE,
                                   edgecolor=S.GRID, linewidth=0.5, zorder=2))
            continue
        v = min(max(val, 0.0), 1.0)
        ax.add_patch(Rectangle((j, y), 1, 1, facecolor=ramp(v),
                               edgecolor="white", linewidth=0.6, zorder=2))
        if txt:
            ax.text(j + 0.5, y + 0.5, txt.replace(".", S.DECIMAL),
                    ha="center", va="center", fontsize=6.2, zorder=3,
                    color="white" if v > 0.62 else S.INK)

ax.set_xlim(0, m); ax.set_ylim(0, n)
ax.set_xticks([]); ax.set_yticks([])
for s in ("top", "right", "left", "bottom"):
    ax.spines[s].set_visible(False)

# column labels, angled so twelve of them fit at 160 mm
for j, (_, lab) in enumerate(PANEL):
    ax.text(j + 0.5, n + 0.25, lab, rotation=38, rotation_mode="anchor",
            ha="left", va="bottom", fontsize=S.PT_TICK, color=S.INK_2)

# row labels: MAG id roman, genus italic -- fig02's convention
for rank, i in enumerate(order):
    r = rows[i]
    y = n - 1 - rank + 0.5
    mag = (r.get("mag_id") or r.get("working_name") or "?").replace("ETO_", "")
    gen = (r.get("genus") or "").strip() or "—"
    ax.text(-6.20, y, mag, ha="left", va="center", fontsize=S.PT_TICK, color=S.INK_2,
            transform=ax.transData)
    ax.text(-0.45, y, gen, ha="right", va="center", fontsize=S.PT_TICK,
            style="italic", color=S.INK_2, transform=ax.transData)
    # domain strip
    ax.add_patch(Rectangle((-0.28, y - 0.5), 0.20, 1.0,
                           facecolor=S.DOMAIN["Archaea" if is_arc(r) else "Bacteria"],
                           edgecolor="none", zorder=3, clip_on=False))

ax.set_title(EL["title"], loc="left", pad=34)

# key: the ramp, and the domain strip it replaces
cax = fig.add_axes([left, (bot_in - 0.20) / H, 0.24, 0.050 / H])
grad = [[i / 255.0 for i in range(256)]]
cax.imshow(grad, aspect="auto", cmap=ramp)
cax.set_xticks([0, 255]); cax.set_xticklabels(["0", "1"], fontsize=S.PT_VALUE)
cax.set_yticks([])
for s in cax.spines.values():
    s.set_visible(False)
# label ABOVE the bar: below it, it lands in the caption
cax.set_title(EL["cbar"], fontsize=S.PT_VALUE, loc="left", pad=4, color=S.INK_2,
              fontweight="normal")

fig.legend(handles=[Patch(facecolor=S.DOMAIN["Bacteria"], label=S.DOMAIN_EL["Bacteria"]),
                    Patch(facecolor=S.DOMAIN["Archaea"], label=S.DOMAIN_EL["Archaea"])],
           loc="lower right", bbox_to_anchor=(right, (bot_in - 0.26) / H),
           ncols=2, handlelength=1.0, handleheight=1.0, borderpad=0.2, columnspacing=1.2)

S.caption(fig, EL["caption"])

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))
if bb.width * 25.4 > FIG_W_MM + 1:
    print("   !! wider than the text block; drop a column or shorten the row labels.")

S.save(fig, "fig03_function_heatmap", OUT)
