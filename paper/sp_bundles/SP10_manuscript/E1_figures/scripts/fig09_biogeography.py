#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig09_biogeography.py -- 29 γονιδιώματα που δεν υπάρχουν πουθενά αλλού.

    python3 fig09_biogeography.py [14_comparative dir] [outdir] [MASTER.tsv]

THE CLAIM. Part 1's strongest single result: not one of the 29 genomes recruits
meaningfully from any public dataset, while every one of them recruits from
Etoliko's own reads. The number that carries it is the public maximum -- across
every MAG and every public dataset, the best breadth anywhere else is 0,008.

FORM. One row per genome, two dots joined by a hairline: breadth in Etoliko's own
reads, and the best breadth in ANY public dataset. A paired plot, not a scatter,
because the public values collapse onto zero and a scatter would pile 29 points
into one corner. Here they form a column against the axis, and the length of each
connector is the result.

LINEAR AXIS, 0 TO 1, AND NO LOG. Breadth is a fraction of a genome covered; zero
is a real value and the distance from 0,008 to 1,0 is the whole point. A log axis
would turn the strongest result in the chapter into a gentle slope.

IT READS THE RUN. 14_comparative/recruitment_matrix.tsv, the same file step 37
read: every column ending in _breadth, split into ETOLIKO_* and everything else.
Nothing is transcribed and no threshold is invented here -- the 30% line is the
one step 37 tested against.
"""
import sys, csv, os, glob
import etoliko_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

args = [a for a in sys.argv[1:] if not a.startswith("--")]
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 170.0
THR = 0.30                         # step 37's test, not a new choice

EL = {
    "title": "",
    "xlabel": "breadth of coverage (fraction of genome)",
    "eto": "Aitoliko reads",
    "pub": "best public dataset",
    "thr": "detection threshold %s",
    "caption": (
        "Το εύρος κάλυψης είναι το κλάσμα του γονιδιώματος που καλύπτεται από "
        "τουλάχιστον μία ανάγνωση, όταν οι αναγνώσεις ενός δείγματος "
        "χαρτογραφηθούν πάνω του. Κάθε γραμμή είναι ένα MAG: ο γεμάτος κύκλος "
        "δίνει το εύρος από τις αναγνώσεις του Αιτωλικού — την καλύτερη από τις "
        "δύο βιβλιοθήκες, ιζήματος και καλλιέργειας — και ο ανοιχτός το καλύτερο "
        "εύρος σε οποιοδήποτε από %d δημόσια μεταγονιδιωματικά σύνολα δεδομένων. "
        "ΟΙ ΑΝΟΙΧΤΟΙ ΚΥΚΛΟΙ ΔΕΝ ΕΧΟΥΝ ΤΟΠΟΘΕΤΗΘΕΙ ΣΤΟ ΜΗΔΕΝ: σχεδιάζονται στις "
        "πραγματικές τους τιμές, που κυμαίνονται από %s έως %s, δηλαδή κάτω από "
        "το ένα εκατοστό του γονιδιώματος, και σε αυτή την κλίμακα δεν "
        "ξεχωρίζουν από τον άξονα. Κάλυψη αυτού του μεγέθους προκύπτει από "
        "λίγες συντηρημένες περιοχές που χαρτογραφούνται διασταυρούμενα από "
        "συγγενικούς οργανισμούς· είναι το υπόβαθρο της μεθόδου και όχι ένδειξη "
        "παρουσίας. Η διακεκομμένη γραμμή στο %s είναι το κριτήριο παρουσίας: "
        "γονιδίωμα που καλύπτεται σε λιγότερο από το ένα τρίτο του μήκους του "
        "δεν θεωρείται ότι ανιχνεύεται. Και τα %d γονιδιώματα το ξεπερνούν στο "
        "Αιτωλικό, με χαμηλότερη τιμή %s, και κανένα δεν το πλησιάζει σε δημόσιο "
        "σύνολο — διαφορά %s φορές ανάμεσα στη χαμηλότερη τιμή εδώ και την "
        "υψηλότερη οπουδήποτε αλλού."),
}


def pick(*cands):
    here = os.path.dirname(os.path.abspath(__file__))
    for c in cands:
        for base in ("", here, os.path.join(here, "..")):
            p = os.path.join(base, c) if base else c
            if os.path.exists(p):
                return os.path.normpath(p)
    return None


CDIR = args[0] if args and os.path.isdir(args[0]) else pick("14_comparative")
if not CDIR:
    sys.exit("fig09: cannot find 14_comparative/. Pass it as the first argument.")
SRC = os.path.join(CDIR, "recruitment_matrix.tsv")
if not os.path.exists(SRC):
    sys.exit("fig09: %s missing -- step 14 has not produced a matrix.\n"
             "       what is in %s:\n         %s"
             % (SRC, CDIR, "\n         ".join(sorted(os.listdir(CDIR))[:20])))
MASTER = args[2] if len(args) > 2 else pick("47_master/MASTER.tsv", "MASTER.tsv")
print("   reading %s" % SRC)

rows = list(csv.DictReader(open(SRC, errors="replace"), delimiter="\t"))
if not rows:
    sys.exit("fig09: %s is empty" % SRC)
cols = [c for c in rows[0] if c]
bcols = [c for c in cols if c.endswith("_breadth")]
ETO = [c for c in bcols if c.startswith("ETOLIKO_")]
PUB = [c for c in bcols if not c.startswith("ETOLIKO_")]
print("   %d MAG(s); %d breadth column(s) = %d Etoliko + %d public"
      % (len(rows), len(bcols), len(ETO), len(PUB)))
if not ETO or not PUB:
    sys.exit("fig09: need both ETOLIKO_* and public *_breadth columns.\n"
             "       columns present:\n         %s" % "\n         ".join(cols))

IDC = next((c for c in cols if c.lower() in
            ("mag_id", "mag", "genome", "bin", "working_name", "id")), cols[0])
print("   identifier column: %s" % IDC)

genus = {}
if MASTER and os.path.exists(MASTER):
    for r in csv.DictReader(open(MASTER, errors="replace"), delimiter="\t"):
        g = (r.get("genus") or "").strip() or "—"
        for k in (r.get("mag_id"), r.get("working_name")):
            if k:
                genus[k] = g
                genus[k.replace("ETO_", "")] = g


def val(r, c):
    try: return float(r[c])
    except (TypeError, ValueError): return float("nan")


def mx(r, cs):
    v = [x for x in (val(r, c) for c in cs) if x == x]
    return max(v) if v else float("nan")


dat = []
for r in rows:
    k = (r.get(IDC) or "?").strip()
    e, p = mx(r, ETO), mx(r, PUB)
    if e != e and p != p:
        continue
    bestpub = max(((val(r, c), c) for c in PUB if val(r, c) == val(r, c)),
                  default=(float("nan"), ""))
    dat.append({"id": k.replace("ETO_", ""), "genus": genus.get(k, genus.get(
        k.replace("ETO_", ""), "—")), "eto": 0.0 if e != e else e,
        "pub": 0.0 if p != p else p, "site": bestpub[1].replace("_breadth", "")})
if not dat:
    sys.exit("fig09: no usable rows")

dat.sort(key=lambda d: -d["eto"])
pubmax = max(d["pub"] for d in dat)
worst = max(dat, key=lambda d: d["pub"])
below = [d["id"] for d in dat if d["eto"] < THR]
above_pub = [d["id"] for d in dat if d["pub"] >= THR]
print("   highest public breadth anywhere: %.4f  (%s, %s)"
      % (pubmax, worst["id"], worst["site"] or "?"))
print("   Etoliko breadth range: %.4f .. %.4f" % (dat[-1]["eto"], dat[0]["eto"]))
print("   genome(s) below the %.0f%% test in Etoliko: %s"
      % (THR * 100, ", ".join(below) or "none"))
print("   genome(s) reaching %.0f%% in a public dataset: %s"
      % (THR * 100, ", ".join(above_pub) or "none"))
etomin = min(d["eto"] for d in dat)
pubmin = min(d["pub"] for d in dat)
FOLD = (etomin / pubmax) if pubmax > 0 else float("inf")
print("   lowest Etoliko breadth %.4f against best public %.4f -- %.0fx"
      % (etomin, pubmax, FOLD))

# ---------------------------------------------------------------- layout
S.apply()
W = FIG_W_MM / 25.4
ROW = 0.150
TOP_IN, BOT_IN = 0.52, 2.12
H = TOP_IN + ROW * len(dat) + BOT_IN
fig = plt.figure(figsize=(W, H))
LEFT = 1.96
ax = fig.add_axes([LEFT / W, BOT_IN / H, (W - LEFT - 0.30) / W,
                   ROW * len(dat) / H])
S.recessive_grid(ax, axis="x")
S.despine(ax, keep=("bottom",))

for i, d in enumerate(dat):
    y = len(dat) - 1 - i
    ax.plot([d["pub"], d["eto"]], [y, y], color=S.RULE, lw=0.8, zorder=2)
    ax.scatter(d["eto"], y, s=26, color=S.DOMAIN["Bacteria"], linewidths=0, zorder=4)
    ax.scatter(d["pub"], y, s=26, facecolors="none", edgecolors=S.INK_MUTED,
               linewidths=0.9, zorder=5)
    # fig02's convention: MAG id roman, genus italic, both in the label column
    ax.text(-(LEFT / (W - LEFT - 0.30)) + 0.015, y, d["id"], ha="left",
            va="center", fontsize=S.PT_TICK, color=S.INK_2,
            transform=ax.get_yaxis_transform(), clip_on=False)
    ax.text(-0.018, y, d["genus"], ha="right", va="center", fontsize=S.PT_TICK,
            style="italic", color=S.INK_2,
            transform=ax.get_yaxis_transform(), clip_on=False)

ax.axvline(THR, color=S.INK, lw=1.0, ls=(0, (4, 2)), zorder=3)
ax.text(THR + 0.012, len(dat) - 0.4, EL["thr"] % S.num(THR, 2), ha="left",
        va="top", fontsize=S.PT_VALUE, color=S.INK, zorder=6)

ax.set_xlim(-0.015, 1.03)
ax.set_ylim(-0.9, len(dat) - 0.1)
ax.set_yticks([])
ax.set_xlabel(EL["xlabel"])
xt = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
ax.set_xticks(xt)
ax.set_xticklabels([S.num(t, 1) for t in xt])

fig.text(0.012, 1 - 0.24 / H, EL["title"], ha="left", va="top",
         fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)

# THE LEGEND GOES OUTSIDE THE AXES. Every filled dot sits between 0,91 and 1,00,
# so the right-hand corners are occupied on all 29 rows, and the left-hand ones
# hold the open dots -- there is no empty corner to put a key in. Below the axis
# label there is, and it costs 4 mm of height.
fig.legend(handles=[
    Line2D([], [], marker="o", linestyle="none", markersize=5.2,
           markerfacecolor=S.DOMAIN["Bacteria"], markeredgewidth=0, label=EL["eto"]),
    Line2D([], [], marker="o", linestyle="none", markersize=5.2,
           markerfacecolor="none", markeredgecolor=S.INK_MUTED,
           markeredgewidth=0.9, label=EL["pub"])],
    loc="lower left", bbox_to_anchor=(LEFT / W, (BOT_IN - 0.62) / H),
    ncols=2, handlelength=1.0, borderpad=0.3, columnspacing=2.0, frameon=False)


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
caption_fit(fig, EL["caption"] % (
    len(PUB), S.num(pubmin, 4), S.num(pubmax, 4), S.num(THR, 2),
    len(dat), S.num(etomin, 3), "%.0f" % FOLD))

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))

S.save(fig, "S_biogeography", OUT)
