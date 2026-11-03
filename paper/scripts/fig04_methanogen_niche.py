#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig04_methanogen_niche.py -- εξοπλισμός έναντι ανιχνευσιμότητας.

    python3 fig04_methanogen_niche.py [MASTER.tsv] [outdir]

THE ARGUMENT THIS PANEL HAS TO CARRY. The lagoon is sulfate-replete, so sulfate
reducers should outcompete methanogens for H2 and acetate, and the way out is the
non-competitive substrates -- methylamines and methanol, which sulfate reducers
cannot use. The catalogue contains methanogens equipped for that route. It also
shows that the best-equipped ones are NOT the ones the sediment reads can see.
Both halves have to be visible in one figure or the argument reads as special
pleading.

FORM. Capability on x, in-situ detectability on y, one dot per archaeal genome,
every dot named. Eight points is small enough that direct labels beat any legend,
and the threshold that defines the in-situ set is drawn because the whole point is
which genomes fall on which side of it.

ENCODED ONCE, EACH:
    x            methylamine + methanol methyltransferase genes
    y            sediment breadth -- can the undisturbed sediment see it at all
    marker area  hydrogenotrophic capacity (fwd/fmd + hydrogenase + fdh)
    the rule     the in-situ cut, taken from the data, not chosen here
    colour       archaea throughout: this panel is one domain, so colour carries
                 nothing and is deliberately uniform

WHAT THE READER SHOULD BE ABLE TO DO. Count the dots above the rule, and see that
the ones furthest right are not among them.
"""
import sys, csv, os
import etoliko_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SRC = args[0] if args else "MASTER.tsv"
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 160.0

METHYL = ("methylamine_MT", "methanol_MT")
HYDRO = ("fwd_fmd", "hydrogenase", "fdh")

EL = {
    "title":  "Μεθανογόνα: εξοπλισμός για μη ανταγωνιστικά υποστρώματα\nέναντι ανιχνευσιμότητας στο ίζημα",
    "xlabel": "γονίδια μεθυλοτρανσφερασών μεθυλαμίνης και μεθανόλης",
    "ylabel": "εύρος κάλυψης στο ίζημα",
    "cut":    "όριο του συνόλου in situ (%s)",
    "size":   "υδρογονοτροφικά\nγονίδια",
    "caption": ("Τα οκτώ αρχαία του καταλόγου. Ο οριζόντιος άξονας μετρά τα γονίδια των "
                "μεθυλοτρανσφερασών μεθυλαμίνης και μεθανόλης — το μη ανταγωνιστικό "
                "υπόστρωμα, που οι θειικοαναγωγείς δεν μπορούν να χρησιμοποιήσουν. Ο "
                "κατακόρυφος άξονας μετρά το εύρος κάλυψης από αναγνώσεις ιζήματος, "
                "δηλαδή αν το αδιατάρακτο ίζημα μπορεί καν να δει το γονιδίωμα. Το "
                "εμβαδόν του κύκλου δίνει την υδρογονοτροφική ικανότητα (fwd/fmd, "
                "υδρογονάση, fdh). Η διακεκομμένη γραμμή είναι το όριο που ορίζει το "
                "σύνολο in situ (σχήμα 2). Τα γονιδιώματα με τον πληρέστερο "
                "μεθυλοτροφικό εξοπλισμό βρίσκονται κάτω από τη γραμμή: ο εξοπλισμός "
                "υπάρχει στη λιμνοθάλασσα, αλλά οι φορείς του είναι σπάνιοι in situ."),
}

if not os.path.exists(SRC):
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, os.path.basename(SRC)),
                 os.path.join(here, "..", os.path.basename(SRC)),
                 os.path.join(here, "..", "47_master", os.path.basename(SRC))):
        if os.path.exists(cand):
            SRC = cand; print("   found %s" % os.path.normpath(cand)); break
    else:
        sys.exit("fig04: cannot find %s -- pass the path as the first argument" % SRC)


def fnum(v, d=0.0):
    try: return float(str(v).replace(",", ".").replace("%", "").strip())
    except (TypeError, ValueError): return d


with open(SRC, encoding="utf-8") as fh:
    allrows = list(csv.DictReader(fh, delimiter="\t"))
if not allrows:
    sys.exit("fig04: %s is empty" % SRC)

need = [c for c in METHYL + HYDRO + ("sed_breadth",) if c not in allrows[0]]
if need:
    sys.exit("fig04: missing column(s): %s\n       columns present:\n         %s"
             % (", ".join(need), "\n         ".join(sorted(allrows[0]))))

rows, cut_pool = [], []
for r in allrows:
    b = fnum(r.get("sed_breadth"), -1.0)
    if str(r.get("in_situ_ten", "")).lower().startswith(("y", "t", "1")) and b >= 0:
        cut_pool.append(b)
    if not str(r.get("domain", "")).lower().startswith("arch"):
        continue
    rows.append({
        "mag": (r.get("mag_id") or r.get("working_name") or "?").replace("ETO_", ""),
        "genus": (r.get("genus") or "").strip() or "—",
        "x": sum(fnum(r.get(c)) for c in METHYL),
        "y": max(b, 0.0),
        "h": sum(fnum(r.get(c)) for c in HYDRO),
    })

if not rows:
    sys.exit("fig04: no archaeal rows found -- check the domain column")
CUT = min(cut_pool) if cut_pool else None
print("   %d archaea; methylotrophic genes %g–%g; breadth %.4f–%.4f"
      % (len(rows), min(r["x"] for r in rows), max(r["x"] for r in rows),
         min(r["y"] for r in rows), max(r["y"] for r in rows)))
if CUT is None:
    print("   !! no in_situ_ten column -- the threshold line is omitted")
else:
    above = [r["mag"] for r in rows if r["y"] >= CUT]
    print("   in-situ cut %.4f; archaea above it: %s"
          % (CUT, ", ".join(above) if above else "none"))
    top = sorted(rows, key=lambda r: -r["x"])[:3]
    print("   best-equipped three: %s"
          % ", ".join("%s (x=%g, breadth %.4f)" % (r["mag"], r["x"], r["y"]) for r in top))

hmax = max([r["h"] for r in rows] or [1.0])
hmin = min([r["h"] for r in rows] or [0.0])


def area(h):
    f = 0.0 if hmax == hmin else (h - hmin) / (hmax - hmin)
    return 40.0 + f * 180.0


# A CHARACTER COUNT IS NOT A WIDTH. S.caption wraps at 112 characters, which is a
# different number of millimetres in every caption: Greek at 7,5 pt runs wider per
# character than the Latin the default was set against, and a caption heavy in
# long words wraps early and leaves the block half empty. Wrapping to a MEASURED
# width gives every caption in the set the same physical measure, which is what
# the 112 was standing in for.
def caption_fit(fig, text, lo=70, hi=150):
    lim = fig.get_figwidth() * (1.0 - 0.012)     # mirror the left inset
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


S.apply()
W = FIG_W_MM / 25.4
fig = plt.figure(figsize=(W, W * 0.78))

# THE CAPTION IS SET FIRST AND THE AXES ARE FITTED ABOVE IT. A fixed axes
# rectangle assumes a caption of a fixed number of lines; this one runs to five
# or eight depending on how it wraps, so a rectangle tuned for one leaves a band
# of white the height of three lines under the other. Measure, then place.
fig.canvas.draw()
print("   caption wrapped at %d characters to fill the text block"
      % caption_fit(fig, EL["caption"]))
fig.canvas.draw()
_ctop = (max(t.get_window_extent(fig.canvas.get_renderer()).y1 for t in fig.texts)
         / fig.dpi / fig.get_figheight())
BOT = min(max(_ctop + 0.105, 0.22), 0.45)      # + room for the axis label
ax = fig.add_axes([0.105, BOT, 0.860, 0.875 - BOT])
S.recessive_grid(ax, axis="both")
S.despine(ax, keep=("left", "bottom"))

xmax = max(r["x"] for r in rows)
XHI = xmax + 1.0          # provisional; the size key claims its own column below
ax.set_xlim(-0.6, XHI)
ax.set_ylim(-0.04, 1.06)

for r in rows:
    ax.scatter(r["x"], r["y"], s=area(r["h"]), color=S.DOMAIN["Archaea"],
               alpha=0.9, linewidths=0, zorder=5)

# THE THRESHOLD CAPTION GOES BELOW ITS LINE. Above it, it crowds whichever genome
# sits nearest the cut -- and the genome nearest the cut is by definition the one
# the line is about, so the collision is guaranteed, not bad luck.
cut_txt = None
if CUT is not None:
    # drawn to the edge of the DATA, not the edge of the axes: past the last
    # genome the rule has nothing to separate, and it would run under the key
    ax.plot([-0.6, xmax + 0.6], [CUT, CUT], color=S.INK, lw=1.0,
            ls=(0, (4, 2)), zorder=3, clip_on=False)
    cut_txt = ax.text(-0.4, CUT - 0.025, EL["cut"] % S.num(CUT), ha="left",
                      va="top", fontsize=S.PT_VALUE, color=S.INK, zorder=4)

# THE SIZE KEY IS DRAWN BEFORE THE LABELS ARE PLACED, ON PURPOSE. It is made of
# circles in the same colour as the data, so without a frame it reads as three
# more genomes -- it gets the page colour behind it and a hairline round it. And
# because it is a real artist by the time the labels are solved, it can be handed
# to the solver as an obstacle instead of being hoped around.
keys = sorted({hmin, (hmin + hmax) / 2, hmax})
leg = ax.legend(handles=[Line2D([], [], marker="o", linestyle="none",
                                markersize=area(k) ** 0.5,
                                markerfacecolor=S.DOMAIN["Archaea"],
                                markeredgewidth=0, alpha=0.9, label=S.num(k, 0))
                         for k in keys],
                loc="upper right", bbox_to_anchor=(0.998, 0.998), title=EL["size"],
                title_fontsize=S.PT_VALUE, labelspacing=1.1, handletextpad=0.8,
                borderpad=0.5, frameon=True, facecolor=S.SURFACE,
                edgecolor=S.RULE, framealpha=1.0)
leg.get_frame().set_linewidth(0.6)
leg.set_zorder(6)

# AND IT GETS ITS OWN COLUMN, MEASURED. Reserving "xmax + 2,4" was a guess, and
# the guess was short: the key overlapped the right-most genome. The axes
# rectangle is fixed, so the key's width in pixels does not depend on the x
# limit -- which makes the limit that puts the key entirely to the right of the
# data solvable in one step instead of tuned by eye.
fig.canvas.draw()
_R = fig.canvas.get_renderer()
AW = ax.get_window_extent(_R).width
LW = leg.get_window_extent(_R).width + 0.10 * fig.dpi     # key + 2,5 mm of air
if LW < AW - 20:
    XHI = (xmax + 1.0) * AW / (AW - LW) - 0.6
    ax.set_xlim(-0.6, XHI)
    print("   size key reserved %.1f x-units; x limit set to %.2f" % (XHI - xmax, XHI))

# DIRECT LABELS, DE-COLLIDED IN DISPLAY SPACE. Eight genomes, several sharing an
# x value of 1, so a fixed offset stacks them. Two changes over fig01's solver:
# every label box is MEASURED rather than estimated from a character count (the
# estimate is what let 017 and 007 print on top of each other -- two-line labels
# in a mixed Greek/Latin face are not 0,52 em per character), and rectangular
# obstacles -- the threshold caption, the size key -- are pushed out of along
# whichever axis they overlap least, so a box in the corner does not shove a
# label down through the data.
fig.canvas.draw()
REND = fig.canvas.get_renderer()
P2PX = fig.dpi / 72.0
LIFT = 13.0 * P2PX


def measured(s):
    t = fig.text(0.5, 0.5, s, fontsize=S.PT_VALUE, ha="center", va="center")
    bb = t.get_window_extent(REND)
    t.remove()
    return bb.width, bb.height


lab = []
for r in rows:
    ax_, ay = ax.transData.transform((r["x"], r["y"]))
    txt = "%s\n%s" % (r["mag"], r["genus"])
    w, h = measured(txt)
    lab.append({"r": r, "t": txt, "ax": ax_, "ay": ay, "x": ax_, "y": ay + LIFT,
                "w": w + 3.0 * P2PX, "h": h + 2.0 * P2PX})

circles = [(ax.transData.transform((r["x"], r["y"]))[0],
            ax.transData.transform((r["x"], r["y"]))[1],
            (area(r["h"]) / 3.14159) ** 0.5 * P2PX + 4.0 * P2PX) for r in rows]
boxes = []
for art in ((cut_txt, leg) if cut_txt is not None else (leg,)):
    bb = art.get_window_extent(REND)
    boxes.append((bb.x0 - 2.0, bb.y0 - 2.0, bb.x1 + 2.0, bb.y1 + 2.0))

x0, y0 = ax.transData.transform((-0.6, -0.04))
x1, y1 = ax.transData.transform((XHI, 1.06))


def clamp():
    for L in lab:
        L["x"] = min(max(L["x"], x0 + L["w"] / 2), x1 - L["w"] / 2)
        L["y"] = min(max(L["y"], y0 + L["h"] / 2), y1 - L["h"] / 2)


def separate(spring, avoid=True):
    moved = False
    for i in range(len(lab)):
        for j in range(i + 1, len(lab)):
            a, b = lab[i], lab[j]
            ox = (a["w"] + b["w"]) / 2 - abs(a["x"] - b["x"])
            oy = (a["h"] + b["h"]) / 2 - abs(a["y"] - b["y"])
            if ox <= 0 or oy <= 0:
                continue
            # SEPARATE ALONG WHICHEVER AXIS IS CHEAPER. Pushing only in y is what
            # kept 017 and 007 welded together: they overlapped by 13 px
            # sideways and 40 px vertically, so the vertical fix had to win
            # against the spring and against 017's own marker, and it lost.
            if ox <= oy:
                d = (ox / 2 + 0.5) * (1 if a["x"] <= b["x"] else -1)
                a["x"] -= d; b["x"] += d
            else:
                d = (oy / 2 + 0.5) * (1 if a["y"] <= b["y"] else -1)
                a["y"] -= d; b["y"] += d
            moved = True
    for L in (lab if avoid else ()):
        for ox_, oy_, rad in circles:
            dy_ = L["y"] - oy_
            if abs(L["x"] - ox_) < rad + L["w"] / 2 and abs(dy_) < rad + L["h"] / 2:
                push = rad + L["h"] / 2 - abs(dy_) + 0.5
                L["y"] += push if dy_ >= 0 else -push
                moved = True
        for bx0, by0, bx1, by1 in boxes:
            lx0, lx1 = L["x"] - L["w"] / 2, L["x"] + L["w"] / 2
            ly0, ly1 = L["y"] - L["h"] / 2, L["y"] + L["h"] / 2
            if lx1 <= bx0 or lx0 >= bx1 or ly1 <= by0 or ly0 >= by1:
                continue
            left, right = lx1 - bx0, bx1 - lx0        # cost of leaving each way
            down, up = ly1 - by0, by1 - ly0
            step = min(left, right, down, up) + 0.5
            if step == left + 0.5:    L["x"] -= step
            elif step == right + 0.5: L["x"] += step
            elif step == down + 0.5:  L["y"] -= step
            else:                     L["y"] += step
            moved = True
    if spring:
        for L in lab:
            L["x"] += (L["ax"] - L["x"]) * 0.03
            L["y"] += (L["ay"] + LIFT - L["y"]) * 0.03
    clamp()
    return moved


for _ in range(900):
    if not separate(True): break
for _ in range(600):
    if not separate(False): break
# LAST RESORT: DROP THE MARKERS FROM THE CONSTRAINT SET. Labels pushed apart and
# markers pushing them back can deadlock at a few pixels of overlap, and of the
# two failures, two labels touching each other is the one a reader cannot parse;
# a label grazing a circle is merely untidy. So the final pass keeps only the
# label-label rule -- and only if the full rule set has not already converged.
for _ in range(400):
    if not separate(False, avoid=False): break

rest = 0
for i in range(len(lab)):
    for j in range(i + 1, len(lab)):
        a, b = lab[i], lab[j]
        if ((a["w"] + b["w"]) / 2 > abs(a["x"] - b["x"])
                and (a["h"] + b["h"]) / 2 > abs(a["y"] - b["y"])):
            rest += 1
print("   label pairs still overlapping after the solver: %d" % rest)

for L in lab:
    dx = (L["x"] - L["ax"]) / P2PX
    dy = (L["y"] - L["ay"]) / P2PX
    far = (dx * dx + dy * dy) ** 0.5 > 12.0
    ax.annotate(L["t"], (L["r"]["x"], L["r"]["y"]), textcoords="offset points",
                xytext=(dx, dy), ha="center", va="center", fontsize=S.PT_VALUE,
                color=S.INK, zorder=8,
                arrowprops=dict(arrowstyle="-", lw=0.55, shrinkA=1, shrinkB=3,
                                color=S.INK_MUTED, alpha=0.75) if far else None)

ax.set_xlabel(EL["xlabel"]); ax.set_ylabel(EL["ylabel"])
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticklabels([S.num(v, 1) for v in (0, 0.2, 0.4, 0.6, 0.8, 1.0)])
ax.set_xticks([t for t in ax.get_xticks() if 0 <= t <= xmax + 0.4])
ax.set_xticklabels([S.num(t, 0) for t in ax.get_xticks()])
ax.set_title(EL["title"], loc="left", pad=8)

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))
if bb.width * 25.4 > FIG_W_MM + 1:
    print("   !! wider than the text block -- shorten the labels or the size key.")

S.save(fig, "fig04_methanogen_niche", OUT)
