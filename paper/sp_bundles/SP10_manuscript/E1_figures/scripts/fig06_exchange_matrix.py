#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig06_exchange_matrix.py -- ποιος τρέφει ποιον, στην ανακατασκευή.

    python3 fig06_exchange_matrix.py [67_network/fig06_matrix.tsv] [outdir] [MASTER.tsv]

WHY A MATRIX AND NOT A NODE-LINK GRAPH. Eleven genomes and thirty-two directed
edges is exactly where a graph stops paying its way: at that density the reader
spends the figure untangling crossings, and the two facts that matter here are
both about WHOLE ROWS AND COLUMNS -- one genome receives from all ten others, two
receive from nobody. An empty column states that in one glance. In a node-link
layout it is a node you have to trace every edge away from to be sure of.

An adjacency matrix is a network diagram. It is the same graph, drawn so that
absence is as visible as presence -- and absence is the result here.

ENCODED ONCE, EACH:
    row          the donor
    column       the receiver
    cell colour  summed SMETANA score for that ordered pair
    cell number  the same value, because colour alone should never carry a number
    grey cell    no predicted transfer -- a real zero, not missing data
    hatched      the diagonal: a genome cannot feed itself, undefined not zero
    margin bars  totals, so the hubs are readable without summing by eye
    domain strip archaea against bacteria, the one thing colour carries elsewhere

ORDER. Both axes are sorted by total score RECEIVED, descending, and the two axes
carry the same order, so the matrix is square and comparable across the diagonal.
The consequence is that the genomes receiving nothing fall together at the end --
that is the point, not a coincidence of sorting.

WHAT THIS PANEL IS A CLAIM ABOUT. Capacity, not activity. Every value is a
SMETANA score from gap-filled draft models on a defined medium, which is a model
of Etoliko and not Etoliko (CHAPTER_STRUCTURE §7). The caption says so.
"""
import sys, csv, os, math, collections
import etoliko_figstyle as S
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Patch

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SRC = args[0] if args else "67_network/fig06_matrix.tsv"
OUT = args[1] if len(args) > 1 else "figures"
MASTER = args[2] if len(args) > 2 else "47_master/MASTER.tsv"
FIG_W_MM = 170.0

EL = {
    "title": "",
    "xlabel": "receiver",
    "ylabel": "donor",
    "cbar": "summed SMETANA score per pair",
    "mdon": "total\ndonated",
    "mrec": "total received",
    "none": "no predicted transfer",
    "diag": "not defined",
    "lowc": "genome < 90% complete (°)",
    "caption": (
        "Γραμμές = δότες, στήλες = δέκτες· κάθε κελί δίνει το άθροισμα των βαθμών "
        "SMETANA για όλες τις ενώσεις που μεταφέρονται από το ένα γονιδίωμα στο "
        "άλλο, σε ζευγαρωτές προσομοιώσεις πάνω στο θρεπτικό μέσο του ιζήματος "
        "Τα γκρίζα κελιά είναι πραγματικά μηδενικά — καμία προβλεπόμενη "
        "μεταφορά — ενώ η διαγώνιος δεν ορίζεται. Οι μπάρες στο περιθώριο δίνουν "
        "τα σύνολα. Και οι δύο άξονες έχουν την ίδια σειρά, κατά φθίνον σύνολο "
        "λήψης. Πρόκειται για ΔΥΝΑΤΟΤΗΤΑ και όχι για δραστηριότητα: οι βαθμοί "
        "προκύπτουν από μοντέλα με συμπλήρωση κενών, πάνω σε ένα ορισμένο "
        "θρεπτικό μέσο που είναι μοντέλο του Αιτωλικού και όχι το Αιτωλικό. Δίπλα "
        "σε κάθε γονιδίωμα δίνεται η πληρότητά του· ένα γονιδίωμα κάτω από 90 % "
        "(°) παράγει μοντέλο του οποίου οι ανάγκες καλύφθηκαν από συμπλήρωση "
        "κενών, οπότε μια κενή γραμμή ή στήλη εκεί δεν ερμηνεύεται."),
}


def find(p, extra=()):
    if os.path.exists(p):
        return p
    here = os.path.dirname(os.path.abspath(__file__))
    b = os.path.basename(p)
    for c in (os.path.join(here, b), os.path.join(here, "..", b),
              os.path.join(here, "..", os.path.dirname(p), b)) + tuple(
              os.path.join(here, "..", d, b) for d in extra):
        if os.path.exists(c):
            print("   found %s" % os.path.normpath(c))
            return c
    return None


SRC = find(SRC, ("67_network",)) or sys.exit(
    "fig06: cannot find the matrix table -- pass it as the first argument")
MASTER = find(MASTER, ("47_master",))


def f(v, d=0.0):
    try: return float(str(v).replace(",", ".").strip())
    except (TypeError, ValueError): return d


edges = []
with open(SRC, encoding="utf-8") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    need = ("donor_mag", "receiver_mag", "sum_smetana")
    miss = [c for c in need if c not in (rdr.fieldnames or [])]
    if miss:
        sys.exit("fig06: %s lacks column(s): %s\n       has: %s"
                 % (SRC, ", ".join(miss), ", ".join(rdr.fieldnames or [])))
    for r in rdr:
        v = f(r["sum_smetana"])
        if v <= 0:
            continue
        edges.append((r["donor_mag"], r["receiver_mag"], v,
                      int(f(r.get("n_compounds"), 0)),
                      r.get("donor_domain", "B"), r.get("receiver_domain", "B")))
if not edges:
    sys.exit("fig06: no non-zero edges in %s -- there is no network to draw, and\n"
             "       that is a sentence for the text, not a panel" % SRC)

# THE NODE SET IS THE EDGE TABLE'S, AND ONLY THE EDGE TABLE'S. The first version
# built it from whatever MASTER.tsv happened to contain and drew a 29x29 grid with
# eighteen empty rows -- MASTER holds the whole catalogue, the SMETANA run held
# eleven. MASTER is a LOOKUP here (genus, domain), never a source of membership.
nodes = sorted({d for d, r, v, n, dd, rd in edges} | {r for d, r, v, n, dd, rd in edges})
dom, genus, compl = {}, {}, {}
for d, r, v, n, dd, rd in edges:
    dom.setdefault(d, "Archaea" if str(dd).upper().startswith("A") else "Bacteria")
    dom.setdefault(r, "Archaea" if str(rd).upper().startswith("A") else "Bacteria")
if MASTER:
    with open(MASTER, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r.get("mag_id") or "").replace("ETO_", "")
            if k in nodes:
                genus[k] = (r.get("genus") or "").strip() or "—"
                try: compl[k] = float(r.get("completeness", ""))
                except (TypeError, ValueError): pass
                if r.get("domain"):
                    dom[k] = ("Archaea" if str(r["domain"]).lower().startswith("arch")
                              else "Bacteria")

S_ = collections.defaultdict(float)
N_ = collections.Counter()
for d, r, v, n, _, _ in edges:
    S_[(d, r)] += v
    N_[(d, r)] += n

rec = {g: sum(v for (a, b), v in S_.items() if b == g) for g in nodes}
don = {g: sum(v for (a, b), v in S_.items() if a == g) for g in nodes}
order = sorted(nodes, key=lambda g: (-rec[g], -don[g], g))
n = len(order)

print("   %d genome(s), %d directed edge(s)" % (n, len(S_)))
print("   receives nothing: %s"
      % (", ".join(g for g in order if rec[g] <= 0) or "none"))
print("   donates nothing : %s"
      % (", ".join(g for g in order if don[g] <= 0) or "none"))
full = [g for g in order if sum(1 for a in order if a != g and S_.get((a, g), 0) > 0) == n - 1]
print("   receives from every other genome: %s" % (", ".join(full) or "none"))
LOWC = 90.0
low = sorted(g for g in nodes if compl.get(g, 100.0) < LOWC)
if low:
    print("   genome(s) below %.0f%% complete, marked on the panel: %s"
          % (LOWC, ", ".join("%s (%.1f%%)" % (g, compl[g]) for g in low)))
    print("   an empty row or column for one of those is not interpretable: a")
    print("   half-complete genome makes a model whose needs were supplied by")
    print("   gap-filling, not by annotation")

VMAX = max(S_.values())
ramp = mcolors.LinearSegmentedColormap.from_list(
    "eto", ["#eef2f5", "#9cc0d8", S.DOMAIN["Bacteria"], "#0c3f61"])

# ---------------------------------------------------------------- layout
S.apply()
W = FIG_W_MM / 25.4
CELL = 0.295                      # inches per cell; 11 cells = 3.25 in = 82 mm
MARG = 0.62                       # inches for the marginal bar strips
LEFT_IN = 2.18                    # row labels: MAG id, genus, completeness
BOT_IN = 1.62                     # receive-strip + axis label + key + caption room

# the cell size is a preference, not a promise: whatever n is, the grid plus its
# label column and margin strip has to fit the text block
FIT = (W - LEFT_IN - MARG - 0.22) / n
if FIT < CELL:
    print("   %d genomes will not fit at %.3f in per cell; using %.3f"
          % (n, CELL, FIT))
    CELL = FIT

# TOP_IN HAS TO HOLD THE ROTATED COLUMN LABELS, AND THEIR HEIGHT DEPENDS ON THE
# LONGEST LABEL, NOT ON A NUMBER I PICKED. A fixed 1,30 in was enough for eleven
# and put the title straight through the column labels the moment they grew.
LAB_ROT = 42.0
_lab_in = (max(len(g) for g in nodes) * 0.54 * S.PT_TICK / 72.0
           * math.sin(math.radians(LAB_ROT)))
TOP_IN = _lab_in + 0.62           # labels + the title's own line + air

GRID_W = CELL * n
H = TOP_IN + GRID_W + MARG + BOT_IN
fig = plt.figure(figsize=(W, H))

ax = fig.add_axes([LEFT_IN / W, (BOT_IN + MARG) / H, GRID_W / W, GRID_W / H])
axd = fig.add_axes([(LEFT_IN + GRID_W + 0.10) / W, (BOT_IN + MARG) / H,
                    MARG / W, GRID_W / H])          # donated, per row
axr = fig.add_axes([LEFT_IN / W, BOT_IN / H,
                    GRID_W / W, (MARG - 0.12) / H])  # received, per column

for a in (ax, axd, axr):
    for s in a.spines.values():
        s.set_visible(False)

# ---------------------------------------------------------------- the cells
for i, a in enumerate(order):            # row = donor, drawn top-down
    y = n - 1 - i
    for j, b in enumerate(order):        # column = receiver
        if a == b:
            ax.add_patch(Rectangle((j, y), 1, 1, facecolor="#f2f1ee",
                                   edgecolor="white", linewidth=0.6, hatch="///",
                                   zorder=2))
            continue
        v = S_.get((a, b), 0.0)
        if v <= 0:
            ax.add_patch(Rectangle((j, y), 1, 1, facecolor=S.GRID,
                                   edgecolor="white", linewidth=0.6, zorder=2))
            continue
        t = v / VMAX
        ax.add_patch(Rectangle((j, y), 1, 1, facecolor=ramp(t),
                               edgecolor="white", linewidth=0.6, zorder=2))
        ax.text(j + 0.5, y + 0.5, S.num(v, 1), ha="center", va="center",
                fontsize=6.2, zorder=3, color="white" if t > 0.55 else S.INK)

ax.set_xlim(0, n); ax.set_ylim(0, n)
ax.set_xticks([]); ax.set_yticks([])

# ---- row labels: MAG id roman, genus italic, then the domain strip
for i, g in enumerate(order):
    y = n - 1 - i + 0.5
    # the label column is LEFT_IN inches wide; in cell units that is LEFT_IN/CELL,
    # so the id starts there rather than at a number picked by eye
    # COMPLETENESS IS ON THE PANEL BECAUSE AN EMPTY COLUMN MEANS TWO DIFFERENT
    # THINGS. MTP4's zero comes from a 99,8 % genome; MAG_025's from a 52,5 % one,
    # whose model needed gap-filling to grow at all. Without the number beside
    # them the reader has no way to weigh one against the other.
    c = compl.get(g)
    mark = " °" if (c is not None and c < LOWC) else ""
    ax.text(-(LEFT_IN / CELL) + 0.10, y, g + mark, ha="left", va="center",
            fontsize=S.PT_TICK, color=S.INK_2)
    ax.text(-(LEFT_IN / CELL) + 3.05, y, genus.get(g, "—"), ha="left",
            va="center", fontsize=S.PT_TICK, style="italic", color=S.INK_2)
    ax.text(-0.40, y, ("%s%%" % S.num(c, 1)) if c is not None else "—",
            ha="right", va="center", fontsize=6.2,
            color=S.HILITE if mark else S.INK_MUTED)
    ax.add_patch(Rectangle((-0.24, y - 0.5), 0.18, 1.0, facecolor=S.DOMAIN[dom[g]],
                           edgecolor="none", zorder=3, clip_on=False))

# ---- column labels, angled; the domain strip repeated so the axes match
collab = []
for j, g in enumerate(order):
    collab.append(ax.text(j + 0.52, n + 0.32,
                          g + (" °" if compl.get(g, 100.0) < LOWC else ""),
                          rotation=LAB_ROT,
                          rotation_mode="anchor", ha="left", va="bottom",
                          fontsize=S.PT_TICK, color=S.INK_2))
    ax.add_patch(Rectangle((j, n + 0.06), 1.0, 0.18, facecolor=S.DOMAIN[dom[g]],
                           edgecolor="none", zorder=3, clip_on=False))

# ---------------------------------------------------------------- margins
DMAX = max(max(don.values()), 1e-9)
for i, g in enumerate(order):
    axd.barh(n - 1 - i + 0.5, don[g], height=0.62, color=S.DOMAIN[dom[g]],
             linewidth=0, zorder=3)
    axd.text(don[g] + DMAX * 0.04, n - 1 - i + 0.5, S.num(don[g], 1), va="center",
             ha="left", fontsize=6.2, color=S.INK_MUTED, zorder=4)
axd.set_ylim(0, n); axd.set_xlim(0, DMAX * 1.42)
axd.set_xticks([]); axd.set_yticks([])
axd.set_title(EL["mdon"], fontsize=S.PT_VALUE, color=S.INK_2, loc="left", pad=4)

RMAX = max(max(rec.values()), 1e-9)
for j, g in enumerate(order):
    axr.bar(j + 0.5, rec[g], width=0.62, color=S.DOMAIN[dom[g]], linewidth=0, zorder=3)
    if rec[g] > 0:
        # the axis is inverted so the bars hang from the matrix: the value goes
        # BELOW the bar end, which after inversion means va="top"
        axr.text(j + 0.5, rec[g] + RMAX * 0.06, S.num(rec[g], 1), ha="center",
                 va="top", fontsize=6.2, color=S.INK_MUTED, zorder=4)
axr.set_xlim(0, n); axr.set_ylim(0, RMAX * 1.30)
axr.set_xticks([]); axr.set_yticks([])
axr.invert_yaxis()
axr.set_ylabel(EL["mrec"], fontsize=S.PT_VALUE, color=S.INK_2, rotation=0,
               ha="right", va="center", labelpad=8)

# THE TITLE IS PLACED ABOVE THE COLUMN LABELS BY MEASUREMENT. A title pad is a
# guess about how tall the rotated labels are, and the guess was wrong as soon as
# the label set changed -- the title printed through them.
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_top = max(t.get_window_extent(_r).y1 for t in collab) / fig.dpi / H
fig.text(LEFT_IN / W, min(_top + 0.012, 0.995), EL["title"], ha="left", va="bottom",
         fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)

# ---------------------------------------------------------------- keys
# THE COLOUR KEY GOES HARD LEFT, THE DOMAIN KEY HARD RIGHT. Both sat in the
# middle of the same line, so the colour bar's caption ran straight into the
# Bacteria swatch and the two keys read as one confused row. Two keys on one line
# need the page's full width between them, not a polite gap.
cax = fig.add_axes([0.14 / W, (BOT_IN - 0.32) / H, 1.05 / W, 0.055 / H])
cax.imshow([[i / 255.0 for i in range(256)]], aspect="auto", cmap=ramp)
cax.set_xticks([0, 255])
cax.set_xticklabels(["0", S.num(VMAX, 1)], fontsize=S.PT_VALUE)
cax.set_yticks([])
for s in cax.spines.values():
    s.set_visible(False)
cax.set_title(EL["cbar"], fontsize=S.PT_VALUE, loc="left", pad=4, color=S.INK_2)

fig.legend(handles=[Patch(facecolor=S.DOMAIN["Bacteria"], label=S.DOMAIN_EL["Bacteria"]),
                    Patch(facecolor=S.DOMAIN["Archaea"], label=S.DOMAIN_EL["Archaea"]),
                    Patch(facecolor=S.GRID, label=EL["none"]),
                    Patch(facecolor="#f2f1ee", hatch="///", label=EL["diag"])],
           loc="lower right",
           bbox_to_anchor=((LEFT_IN + GRID_W + MARG + 0.10) / W, (BOT_IN - 0.44) / H),
           ncols=2, handlelength=1.1, handleheight=1.1, borderpad=0.3,
           columnspacing=1.3, frameon=False)


# A CHARACTER COUNT IS NOT A WIDTH -- see fig04. Wrap to a measured width so every
# caption in the set has the same physical measure.
def caption_fit(fig, text, lo=70, hi=150):
    lim = fig.get_figwidth() * (1.0 - 0.012)
    best = lo
    for w in range(lo, hi + 1, 2):
        S.caption(fig, text, width=w)
        art = fig.texts[-1]
        fits = art.get_window_extent(fig.canvas.get_renderer()).x1 / fig.dpi <= lim
        art.remove()
        if not fits:
            break
        best = w
    S.caption(fig, text, width=best)
    return best


fig.canvas.draw()
print("   caption wrapped at %d characters"
      % caption_fit(fig, EL["caption"]))

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))
if bb.width * 25.4 > FIG_W_MM + 1:
    print("   !! wider than the text block -- shorten the row labels or CELL")

S.save(fig, "S_exchange_matrix", OUT)
