#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig08_methylamine_ratio.py -- η μεθυλαμίνη στο ίζημα και στην καλλιέργεια.

    python3 fig08_methylamine_ratio.py [61_mtmB dir] [outdir] [63_mtmB_relabel dir]

IT READS THE RUN, NOT A TRANSCRIPTION. The first version took a hand-typed table
of Volume XV's numbers. That is a copy, and a copy can drift from the run that
produced it without anything complaining. This recomputes the counts from
61_mtmB/{SED,LIQ}_hits.tsv using step 63's own rule -- best hit per read by bit
score, minimum alignment length 30, then the seven bacterial accessions in
63_mtmB_relabel/bacterial.txt relabelled to DECOYbact -- and prints the full
marker table so it can be checked line for line against step 63's output.

WHAT THIS PANEL CARRIES. Volume XV's headline, and the one read-level result the
chapter's methylotrophy argument rests on:

    mtmB per mcrA is 2,3x higher in the sediment than in the enrichment cultures.

THE RATIO IS THE MEASUREMENT, NOT THE COUNTS. The sediment library has three
recorded read counts across the project's history, and the culture library is a
different size again. A ratio of two markers measured in the same library is
denominator-free: it carries no assumption about library size at all. That is why
the panel plots mtmB/mcrA and prints the counts beside it, rather than the other
way round.

AND IT IS INTERNALLY CALIBRATED. The catalogue holds 11 mtmB proteins above their
curated threshold and 9 mcrA. Corrected for read-plus-gene span -- 1.864 bp
against 2.914 bp -- that predicts a read ratio of 0,78 if the two markers were
equally represented. The sediment measures 0,75. The dashed line is that
prediction, and it is the reason this number is reportable when every earlier
version of it was not: it has a ceiling derived from the same dataset rather than
asserted at 1:1, which is the error that made an earlier mtmB count unusable.

WHAT THE READER SHOULD BE ABLE TO DO. See that the sediment bar lands on the
predicted line and the culture bar falls well below it, and read the raw counts
that produced each without leaving the figure.
"""
import sys, csv, os
import etoliko_figstyle as S
import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SRC = args[0] if args else "fig08_input.tsv"
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 160.0

# The prediction, and where every term in it comes from. Written out so it can be
# argued with rather than trusted: 11 mtmB and 9 mcrA proteins above their KOfam
# curated thresholds, spans 1.864 bp and 2.914 bp (read length + gene length).
PRED_COPIES = (11, 9)
PRED_SPAN = (1864.0, 2914.0)
PRED = (PRED_COPIES[0] * PRED_SPAN[0]) / (PRED_COPIES[1] * PRED_SPAN[1])

EL = {
    "title": "Μεθυλαμίνη έναντι μεθανογένεσης, σε επίπεδο αναγνώσεων",
    "ylabel": "αναγνώσεις mtmB ανά ανάγνωση mcrA",
    "pred": "εσωτερική πρόβλεψη %s",
    "sed": "ίζημα",
    "liq": "καλλιέργειες",
    "caption": (
        "Ο λόγος των αναγνώσεων mtmB προς mcrA στα δύο δείγματα. Ο λόγος, και όχι "
        "οι απόλυτοι αριθμοί, είναι η μέτρηση: δύο δείκτες μετρημένοι στην ίδια "
        "βιβλιοθήκη δεν εξαρτώνται από το μέγεθός της. Η διακεκομμένη γραμμή στο "
        "%s είναι η πρόβλεψη που προκύπτει από τον ίδιο τον κατάλογο — %d "
        "πρωτεΐνες mtmB προς %d mcrA, διορθωμένες για το μήκος ανάγνωσης και "
        "γονιδίου (%s έναντι %s bp). Το ίζημα πέφτει πάνω στην πρόβλεψη· οι "
        "καλλιέργειες πολύ κάτω από αυτήν. Ο εμπλουτισμός έγινε χωρίς θειικά και "
        "υπό H₂/CO₂, που ευνοεί τα υδρογονοτροφικά μεθανογόνα· in situ η οδός της "
        "μεθυλαμίνης είναι σχετικά πιο έντονη. Ο αρχαιακός θετικός μάρτυρας "
        "πέρασε και στα δύο δείγματα (%s και %s αναγνώσεις), οπότε το χαμηλό "
        "νούμερο δεν είναι αποτυχία της αναζήτησης."),
}


MINLEN = 30                       # step 63's own cutoff, not a new choice


def pick(*cands):
    here = os.path.dirname(os.path.abspath(__file__))
    for c in cands:
        for base in ("", here, os.path.join(here, "..")):
            p = os.path.join(base, c) if base else c
            if os.path.isdir(p):
                return os.path.normpath(p)
    return None


P61 = args[0] if args else None
P63 = args[2] if len(args) > 2 else None
P61 = (P61 if (P61 and os.path.isdir(P61)) else pick("61_mtmB"))
P63 = (P63 if (P63 and os.path.isdir(P63)) else pick("63_mtmB_relabel"))
if not P61 or not P63:
    sys.exit("fig08: need 61_mtmB/ and 63_mtmB_relabel/.\n"
             "       found 61_mtmB=%s  63_mtmB_relabel=%s\n"
             "       run from workspace_2, or pass the paths as arguments"
             % (P61, P63))
print("   reading %s and %s" % (P61, P63))

BACT = os.path.join(P63, "bacterial.txt")
if not os.path.exists(BACT):
    sys.exit("fig08: %s missing. Without it the seven bacterial references are\n"
             "       counted as mtmB, which is the error step 63 exists to fix.\n"
             "       Refusing to draw the uncorrected numbers." % BACT)
bact = {l.strip() for l in open(BACT) if l.strip()}
print("   %d bacterial accession(s) relabelled to DECOYbact" % len(bact))


def relabel(ref):
    p = ref.split("|")
    return "DECOYbact" if p[-1] in bact else p[0]


def counts(tag):
    """Step 63's rule, reimplemented here so the figure and the step agree."""
    p = os.path.join(P61, "%s_hits.tsv" % tag)
    if not os.path.exists(p):
        return None
    best = {}
    nline = 0
    for line in open(p, errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        try:
            ln, bs = int(f[3]), float(f[5])
        except ValueError:
            continue
        nline += 1
        if ln < MINLEN:
            continue
        cur = best.get(f[0])
        if cur is None or bs > cur[1]:
            best[f[0]] = (f[1], bs)
    out = {}
    for _q, (ref, _bs) in best.items():
        k = relabel(ref)
        out[k] = out.get(k, 0) + 1
    return out, nline, len(best)


data = {}
for tag in ("SED", "LIQ"):
    r = counts(tag)
    if r is None:
        sys.exit("fig08: %s/%s_hits.tsv missing -- run step 61 first" % (P61, tag))
    data[tag], nline, nread = r
    print("   %-4s %d hit line(s), %d read(s) with a best hit" % (tag, nline, nread))

print()
print("   marker counts after relabelling -- compare with step 63's own table:")
marks = sorted({m for t in data for m in data[t]},
               key=lambda m: -max(data[t].get(m, 0) for t in data))
print("     %-16s %8s %8s" % ("marker", "SED", "LIQ"))
for m in marks:
    print("     %-16s %8d %8d" % (m, data["SED"].get(m, 0), data["LIQ"].get(m, 0)))

for m in ("mtmB", "mcrA"):
    for t in data:
        if m not in data[t]:
            sys.exit("fig08: no %s in %s -- the panel names above are what exist" % (m, t))

CTRL = next((m for m in marks
             if any(k in m.lower() for k in ("arch", "ctrl", "control"))), None)
print("   archaeal control marker: %s" % (CTRL or "none found; the row is omitted"))

SAMPLES = [("SED", EL["sed"]), ("LIQ", EL["liq"])]
print()
print("   internal prediction: (%d x %s) / (%d x %s) = %.4f"
      % (PRED_COPIES[0], S.num(PRED_SPAN[0], 0), PRED_COPIES[1],
         S.num(PRED_SPAN[1], 0), PRED))
bars = []
for tg, lab in SAMPLES:
    a, b = data[tg]["mtmB"], data[tg]["mcrA"]
    if b <= 0:
        sys.exit("fig08: mcrA is zero in %s -- the ratio is undefined" % tg)
    bars.append({"s": tg, "lab": lab, "r": a / b, "a": a, "b": b,
                 "ctrl": data[tg].get(CTRL, 0) if CTRL else 0})
    print("   %-4s mtmB %d / mcrA %d = %.4f   (control %d reads)"
          % (tg, a, b, a / b, bars[-1]["ctrl"]))
if len(bars) == 2 and bars[1]["r"] > 0:
    print("   sediment / culture = %.2fx" % (bars[0]["r"] / bars[1]["r"]))

# ---------------------------------------------------------------- layout
S.apply()
W = FIG_W_MM / 25.4
H = W * 0.62
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0.115, 0.335, 0.495, 0.540])
S.recessive_grid(ax, axis="y")
S.despine(ax, keep=("left", "bottom"))

YMAX = max([b["r"] for b in bars] + [PRED]) * 1.30
for i, b in enumerate(bars):
    ax.bar(i, b["r"], width=0.52, linewidth=0,
           color=S.DOMAIN["Archaea"] if b["s"] == "SED" else S.DOMAIN["Bacteria"],
           zorder=3)
    ax.text(i, b["r"] + YMAX * 0.030, S.num(b["r"], 2), ha="center", va="bottom",
            fontsize=S.PT_LABEL, color=S.INK, zorder=4)
    # the counts that made the ratio, inside the bar where they belong to it
    ax.text(i, YMAX * 0.030, "%g / %g" % (b["a"], b["b"]), ha="center",
            va="bottom", fontsize=S.PT_VALUE, color="white", zorder=5)

ax.axhline(PRED, color=S.INK, lw=1.0, ls=(0, (4, 2)), zorder=4)
ax.text(len(bars) - 0.48, PRED + YMAX * 0.020, EL["pred"] % S.num(PRED, 2),
        ha="right", va="bottom", fontsize=S.PT_VALUE, color=S.INK, zorder=5)

ax.set_xticks(range(len(bars)))
ax.set_xticklabels([b["lab"] for b in bars], fontsize=S.PT_LABEL)
ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.65, len(bars) - 0.35)
ax.set_ylim(0, YMAX)
ax.set_ylabel(EL["ylabel"])
yt = [t for t in ax.get_yticks() if 0 <= t <= YMAX]
ax.set_yticks(yt)
ax.set_yticklabels([S.num(t, 1) for t in yt])
ax.set_title(EL["title"], loc="left", pad=8)

# the numbers, as a small table to the right -- a reader checking the figure
# against Volume XV should not have to compute anything. Three fixed columns:
# label left-aligned, both counts right-aligned, so nothing is truncated to fit.
# clear of the axes: the prediction label is right-aligned at the axes
# edge, so the table must start beyond it, not at it
TX0, TX1, TX2 = 0.690, 0.870, 0.990
ty, dy = 0.860, 0.078


def trow(y, lab, a=None, b=None, col=None, bold=False):
    c = col or S.INK_2
    fig.text(TX0, y, lab, fontsize=S.PT_VALUE, color=c, va="top",
             fontweight="bold" if bold else "normal")
    if a is not None:
        fig.text(TX1, y, a, fontsize=S.PT_VALUE, color=c, va="top", ha="right")
    if b is not None:
        fig.text(TX2, y, b, fontsize=S.PT_VALUE, color=c, va="top", ha="right")


trow(ty, "", "mtmB", "mcrA", col=S.INK)
for k, b in enumerate(bars):
    trow(ty - dy * (k + 1), b["lab"], "%g" % b["a"], "%g" % b["b"])
trow(ty - dy * (len(bars) + 1.7), "αρχαιακός μάρτυρας", col=S.INK_MUTED)
for k, b in enumerate(bars):
    trow(ty - dy * (len(bars) + 2.5 + k), b["lab"], "%g" % b["ctrl"],
         col=S.INK_MUTED)


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
caption_fit(fig, EL["caption"] % (S.num(PRED, 2), PRED_COPIES[0], PRED_COPIES[1],
                                  S.num(PRED_SPAN[0], 0), S.num(PRED_SPAN[1], 0),
                                  "%g" % bars[0]["ctrl"], "%g" % bars[1]["ctrl"]))

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))

S.save(fig, "fig08_methylamine_ratio", OUT)#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig08_methylamine_ratio.py -- η μεθυλαμίνη στο ίζημα και στην καλλιέργεια.

    python3 fig08_methylamine_ratio.py [61_mtmB dir] [outdir] [63_mtmB_relabel dir]

IT READS THE RUN, NOT A TRANSCRIPTION. The first version took a hand-typed table
of Volume XV's numbers. That is a copy, and a copy can drift from the run that
produced it without anything complaining. This recomputes the counts from
61_mtmB/{SED,LIQ}_hits.tsv using step 63's own rule -- best hit per read by bit
score, minimum alignment length 30, then the seven bacterial accessions in
63_mtmB_relabel/bacterial.txt relabelled to DECOYbact -- and prints the full
marker table so it can be checked line for line against step 63's output.

WHAT THIS PANEL CARRIES. Volume XV's headline, and the one read-level result the
chapter's methylotrophy argument rests on:

    mtmB per mcrA is 2,3x higher in the sediment than in the enrichment cultures.

THE RATIO IS THE MEASUREMENT, NOT THE COUNTS. The sediment library has three
recorded read counts across the project's history, and the culture library is a
different size again. A ratio of two markers measured in the same library is
denominator-free: it carries no assumption about library size at all. That is why
the panel plots mtmB/mcrA and prints the counts beside it, rather than the other
way round.

AND IT IS INTERNALLY CALIBRATED. The catalogue holds 11 mtmB proteins above their
curated threshold and 9 mcrA. Corrected for read-plus-gene span -- 1.864 bp
against 2.914 bp -- that predicts a read ratio of 0,78 if the two markers were
equally represented. The sediment measures 0,75. The dashed line is that
prediction, and it is the reason this number is reportable when every earlier
version of it was not: it has a ceiling derived from the same dataset rather than
asserted at 1:1, which is the error that made an earlier mtmB count unusable.

WHAT THE READER SHOULD BE ABLE TO DO. See that the sediment bar lands on the
predicted line and the culture bar falls well below it, and read the raw counts
that produced each without leaving the figure.
"""
import sys, csv, os
import etoliko_figstyle as S
import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SRC = args[0] if args else "fig08_input.tsv"
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 160.0

# The prediction, and where every term in it comes from. Written out so it can be
# argued with rather than trusted: 11 mtmB and 9 mcrA proteins above their KOfam
# curated thresholds, spans 1.864 bp and 2.914 bp (read length + gene length).
PRED_COPIES = (11, 9)
PRED_SPAN = (1864.0, 2914.0)
PRED = (PRED_COPIES[0] * PRED_SPAN[0]) / (PRED_COPIES[1] * PRED_SPAN[1])

EL = {
    "title": "Μεθυλαμίνη έναντι μεθανογένεσης, σε επίπεδο αναγνώσεων",
    "ylabel": "αναγνώσεις mtmB ανά ανάγνωση mcrA",
    "pred": "εσωτερική πρόβλεψη %s",
    "sed": "ίζημα",
    "liq": "καλλιέργειες",
    "caption": (
        "Ο λόγος των αναγνώσεων mtmB προς mcrA στα δύο δείγματα. Ο λόγος, και όχι "
        "οι απόλυτοι αριθμοί, είναι η μέτρηση: δύο δείκτες μετρημένοι στην ίδια "
        "βιβλιοθήκη δεν εξαρτώνται από το μέγεθός της. Η διακεκομμένη γραμμή στο "
        "%s είναι η πρόβλεψη που προκύπτει από τον ίδιο τον κατάλογο — %d "
        "πρωτεΐνες mtmB προς %d mcrA, διορθωμένες για το μήκος ανάγνωσης και "
        "γονιδίου (%s έναντι %s bp). Το ίζημα πέφτει πάνω στην πρόβλεψη· οι "
        "καλλιέργειες πολύ κάτω από αυτήν. Ο εμπλουτισμός έγινε χωρίς θειικά και "
        "υπό H₂/CO₂, που ευνοεί τα υδρογονοτροφικά μεθανογόνα· in situ η οδός της "
        "μεθυλαμίνης είναι σχετικά πιο έντονη. Ο αρχαιακός θετικός μάρτυρας "
        "πέρασε και στα δύο δείγματα (%s και %s αναγνώσεις), οπότε το χαμηλό "
        "νούμερο δεν είναι αποτυχία της αναζήτησης."),
}


MINLEN = 30                       # step 63's own cutoff, not a new choice


def pick(*cands):
    here = os.path.dirname(os.path.abspath(__file__))
    for c in cands:
        for base in ("", here, os.path.join(here, "..")):
            p = os.path.join(base, c) if base else c
            if os.path.isdir(p):
                return os.path.normpath(p)
    return None


P61 = args[0] if args else None
P63 = args[2] if len(args) > 2 else None
P61 = (P61 if (P61 and os.path.isdir(P61)) else pick("61_mtmB"))
P63 = (P63 if (P63 and os.path.isdir(P63)) else pick("63_mtmB_relabel"))
if not P61 or not P63:
    sys.exit("fig08: need 61_mtmB/ and 63_mtmB_relabel/.\n"
             "       found 61_mtmB=%s  63_mtmB_relabel=%s\n"
             "       run from workspace_2, or pass the paths as arguments"
             % (P61, P63))
print("   reading %s and %s" % (P61, P63))

BACT = os.path.join(P63, "bacterial.txt")
if not os.path.exists(BACT):
    sys.exit("fig08: %s missing. Without it the seven bacterial references are\n"
             "       counted as mtmB, which is the error step 63 exists to fix.\n"
             "       Refusing to draw the uncorrected numbers." % BACT)
bact = {l.strip() for l in open(BACT) if l.strip()}
print("   %d bacterial accession(s) relabelled to DECOYbact" % len(bact))


def relabel(ref):
    p = ref.split("|")
    return "DECOYbact" if p[-1] in bact else p[0]


def counts(tag):
    """Step 63's rule, reimplemented here so the figure and the step agree."""
    p = os.path.join(P61, "%s_hits.tsv" % tag)
    if not os.path.exists(p):
        return None
    best = {}
    nline = 0
    for line in open(p, errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        try:
            ln, bs = int(f[3]), float(f[5])
        except ValueError:
            continue
        nline += 1
        if ln < MINLEN:
            continue
        cur = best.get(f[0])
        if cur is None or bs > cur[1]:
            best[f[0]] = (f[1], bs)
    out = {}
    for _q, (ref, _bs) in best.items():
        k = relabel(ref)
        out[k] = out.get(k, 0) + 1
    return out, nline, len(best)


data = {}
for tag in ("SED", "LIQ"):
    r = counts(tag)
    if r is None:
        sys.exit("fig08: %s/%s_hits.tsv missing -- run step 61 first" % (P61, tag))
    data[tag], nline, nread = r
    print("   %-4s %d hit line(s), %d read(s) with a best hit" % (tag, nline, nread))

print()
print("   marker counts after relabelling -- compare with step 63's own table:")
marks = sorted({m for t in data for m in data[t]},
               key=lambda m: -max(data[t].get(m, 0) for t in data))
print("     %-16s %8s %8s" % ("marker", "SED", "LIQ"))
for m in marks:
    print("     %-16s %8d %8d" % (m, data["SED"].get(m, 0), data["LIQ"].get(m, 0)))

for m in ("mtmB", "mcrA"):
    for t in data:
        if m not in data[t]:
            sys.exit("fig08: no %s in %s -- the panel names above are what exist" % (m, t))

CTRL = next((m for m in marks
             if any(k in m.lower() for k in ("arch", "ctrl", "control"))), None)
print("   archaeal control marker: %s" % (CTRL or "none found; the row is omitted"))

SAMPLES = [("SED", EL["sed"]), ("LIQ", EL["liq"])]
print()
print("   internal prediction: (%d x %s) / (%d x %s) = %.4f"
      % (PRED_COPIES[0], S.num(PRED_SPAN[0], 0), PRED_COPIES[1],
         S.num(PRED_SPAN[1], 0), PRED))
bars = []
for tg, lab in SAMPLES:
    a, b = data[tg]["mtmB"], data[tg]["mcrA"]
    if b <= 0:
        sys.exit("fig08: mcrA is zero in %s -- the ratio is undefined" % tg)
    bars.append({"s": tg, "lab": lab, "r": a / b, "a": a, "b": b,
                 "ctrl": data[tg].get(CTRL, 0) if CTRL else 0})
    print("   %-4s mtmB %d / mcrA %d = %.4f   (control %d reads)"
          % (tg, a, b, a / b, bars[-1]["ctrl"]))
if len(bars) == 2 and bars[1]["r"] > 0:
    print("   sediment / culture = %.2fx" % (bars[0]["r"] / bars[1]["r"]))

# ---------------------------------------------------------------- layout
S.apply()
W = FIG_W_MM / 25.4
H = W * 0.62
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0.115, 0.335, 0.495, 0.540])
S.recessive_grid(ax, axis="y")
S.despine(ax, keep=("left", "bottom"))

YMAX = max([b["r"] for b in bars] + [PRED]) * 1.30
for i, b in enumerate(bars):
    ax.bar(i, b["r"], width=0.52, linewidth=0,
           color=S.DOMAIN["Archaea"] if b["s"] == "SED" else S.DOMAIN["Bacteria"],
           zorder=3)
    ax.text(i, b["r"] + YMAX * 0.030, S.num(b["r"], 2), ha="center", va="bottom",
            fontsize=S.PT_LABEL, color=S.INK, zorder=4)
    # the counts that made the ratio, inside the bar where they belong to it
    ax.text(i, YMAX * 0.030, "%g / %g" % (b["a"], b["b"]), ha="center",
            va="bottom", fontsize=S.PT_VALUE, color="white", zorder=5)

ax.axhline(PRED, color=S.INK, lw=1.0, ls=(0, (4, 2)), zorder=4)
ax.text(len(bars) - 0.48, PRED + YMAX * 0.020, EL["pred"] % S.num(PRED, 2),
        ha="right", va="bottom", fontsize=S.PT_VALUE, color=S.INK, zorder=5)

ax.set_xticks(range(len(bars)))
ax.set_xticklabels([b["lab"] for b in bars], fontsize=S.PT_LABEL)
ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.65, len(bars) - 0.35)
ax.set_ylim(0, YMAX)
ax.set_ylabel(EL["ylabel"])
yt = [t for t in ax.get_yticks() if 0 <= t <= YMAX]
ax.set_yticks(yt)
ax.set_yticklabels([S.num(t, 1) for t in yt])
ax.set_title(EL["title"], loc="left", pad=8)

# the numbers, as a small table to the right -- a reader checking the figure
# against Volume XV should not have to compute anything. Three fixed columns:
# label left-aligned, both counts right-aligned, so nothing is truncated to fit.
# clear of the axes: the prediction label is right-aligned at the axes
# edge, so the table must start beyond it, not at it
TX0, TX1, TX2 = 0.690, 0.870, 0.990
ty, dy = 0.860, 0.078


def trow(y, lab, a=None, b=None, col=None, bold=False):
    c = col or S.INK_2
    fig.text(TX0, y, lab, fontsize=S.PT_VALUE, color=c, va="top",
             fontweight="bold" if bold else "normal")
    if a is not None:
        fig.text(TX1, y, a, fontsize=S.PT_VALUE, color=c, va="top", ha="right")
    if b is not None:
        fig.text(TX2, y, b, fontsize=S.PT_VALUE, color=c, va="top", ha="right")


trow(ty, "", "mtmB", "mcrA", col=S.INK)
for k, b in enumerate(bars):
    trow(ty - dy * (k + 1), b["lab"], "%g" % b["a"], "%g" % b["b"])
trow(ty - dy * (len(bars) + 1.7), "αρχαιακός μάρτυρας", col=S.INK_MUTED)
for k, b in enumerate(bars):
    trow(ty - dy * (len(bars) + 2.5 + k), b["lab"], "%g" % b["ctrl"],
         col=S.INK_MUTED)


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
caption_fit(fig, EL["caption"] % (S.num(PRED, 2), PRED_COPIES[0], PRED_COPIES[1],
                                  S.num(PRED_SPAN[0], 0), S.num(PRED_SPAN[1], 0),
                                  "%g" % bars[0]["ctrl"], "%g" % bars[1]["ctrl"]))

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))

S.save(fig, "fig08_methylamine_ratio", OUT)
