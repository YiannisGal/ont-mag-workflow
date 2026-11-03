#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig01_mag_quality.py -- η ποιότητα των 29 MAGs, κατά MIMAG.

    python3 fig01_mag_quality.py [MASTER.tsv] [outdir]

FORM. Two quantities that jointly define a standard (completeness, contamination)
-> position on both axes, with the standard's own boundaries drawn as the
reference. A bar chart of completeness alone would hide contamination, which is
the axis that actually demotes a genome; two bar charts would make the reader do
the join by eye.

ENCODED ONCE, EACH:
    x            completeness
    y            contamination
    marker area  genome size          (area, not radius -- radius would triple
                 the apparent difference)
    colour       DOMAIN only
    shaded band  the MIMAG tier, which is the datum this figure is about
    direct label the high-quality genomes and anything near a boundary

WHY THE BOUNDARIES ARE DRAWN AND NOT ASSUMED. MIMAG high-quality is >=90%
complete and <5% contaminated; medium is >=50% and <10%. Those two lines are
what the 8/21/0 split in the text means, and a reader should be able to count
the dots above the line and get 8.

WIDTH. Authored at FIG_W_MM so the fonts land at their stated point size in the
thesis. The script MEASURES the rendered width after layout and prints it: a
figure exported wider than the text block gets scaled down on insert and every
font in it shrinks with it, which is the one figure defect that cannot be fixed
by editing the caption.
"""
import sys, csv, os
import etoliko_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

SRC = sys.argv[1] if len(sys.argv) > 1 else "MASTER.tsv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures"
FIG_W_MM = 170.0
SHOW_COLUMNS = "--columns" in sys.argv      # print the header and stop

# WHICH GENOMES GET NAMED ON THE PANEL.
# Default: none. Twelve of the twenty-nine sit inside the box, in one corner, and
# labelling all twelve produces a block of text with leader lines crossing the
# data -- three renders confirmed it. Position already carries the message, and
# the twelve are named in the accompanying table where they can be read.
# Put the two or three the text actually argues about here, e.g.
#     HILITE_MAGS = ("LIQ_MAG_003", "LIQ_MAG_011", "LIQ_MAG_022")
# Names are as they appear on the panel: the "ETO_" prefix is stripped.
HILITE_MAGS = ()

EL = {
    "title":  "",
    "xlabel": "Completeness (%)",
    "ylabel": "Contamination (%)",
    # NOT "υψηλή ποιότητα": MIMAG high-quality also requires 23S, 16S and 5S rRNA
    # genes plus >= 18 tRNAs. This box is the completeness/contamination half only,
    # and it must not be labelled as if it were the whole standard.
    "high":   "≥ 90% complete, < 5% contaminated",
    "med":    "medium quality",
    "size":   "genome size",
    "caption": ("Κάθε κύκλος είναι ένα MAG· η θέση δίνει την πληρότητα και την επιμόλυνση "
                "κατά CheckM2 (06_binqc), το εμβαδόν του κύκλου το μέγεθος του γονιδιώματος "
                "και το χρώμα την επικράτεια. Οι γραμμές είναι τα ίδια τα όρια του προτύπου "
                "MIMAG: υψηλή ποιότητα ≥ 90 %% πληρότητα και < 5 %% επιμόλυνση, μέτρια "
                "ποιότητα ≥ 50 %% και < 10 %%. Η κατάταξη φαίνεται από τη θέση ως προς τις "
                "γραμμές, όχι από το χρώμα. %s γονιδιώματα πληρούν τα δύο αυτά κατώφλια "
                "και %s τα κατώφλια της μέτριας ποιότητας. Το πρότυπο MIMAG απαιτεί "
                "επιπλέον, για την υψηλή ποιότητα, την παρουσία των γονιδίων rRNA 23S, "
                "16S και 5S και τουλάχιστον 18 tRNA· ο αριθμός των γονιδιωμάτων υψηλής "
                "ποιότητας είναι επομένως μικρότερος και δίνεται στο κείμενο."),
}

# ---------------------------------------------------------------- columns
# Never assume a column name. Resolve from aliases, and if a required one is
# missing print every column the file HAS -- a figure script that dies saying
# "KeyError: completeness" tells you nothing you can act on.
ALIASES = {
    "mag":    ("mag_id", "working_name", "MAG", "genome", "bin_id", "name"),
    "genus":  ("genus", "g", "gtdb_genus"),
    "domain": ("domain", "Domain", "gtdb_domain", "kingdom"),
    "comp":   ("completeness", "Completeness", "checkm2_completeness",
               "completeness_checkm2", "comp", "Completeness_General"),
    "cont":   ("contamination", "Contamination", "checkm2_contamination",
               "contamination_checkm2", "contam"),
    "size":   ("bp", "genome_size", "size_bp", "total_length", "assembly_size",
               "length", "size", "Genome_Size", "total_bp"),
    # the authoritative tier, if the table carries one. MIMAG high-quality is NOT
    # just >=90/<5 -- it also needs 23S + 16S + 5S rRNA and >= 18 tRNAs -- so a
    # recomputed tier will over-count. Prefer the column; verify it against the
    # recomputation and show the difference rather than hiding it.
    "mimag":  ("mimag", "MIMAG", "quality", "mimag_tier", "tier"),
    "rrna16": ("has_16S", "has_16s", "16S", "ssu", "rrna_16s"),
}
OPTIONAL = {"genus", "size", "mimag", "rrna16"}


def resolve(fieldnames):
    got, missing = {}, []
    low = {f.lower(): f for f in fieldnames}
    for key, cands in ALIASES.items():
        hit = next((low[c.lower()] for c in cands if c.lower() in low), None)
        if hit: got[key] = hit
        elif key not in OPTIONAL: missing.append(key)
        else: got[key] = None
    return got, missing


def fnum(v):
    try:
        x = float(str(v).replace(",", ".").replace("%", "").strip())
        return x
    except (TypeError, ValueError):
        return None


# LOCATE THE TABLE. Run from the workspace root, "MASTER.tsv" is not in the cwd
# and the traceback says only FileNotFoundError, which does not tell you where to
# look. Search the obvious places and say what was tried.
if not os.path.exists(SRC):
    here = os.path.dirname(os.path.abspath(__file__))
    tried = [SRC]
    for cand in (os.path.join(here, os.path.basename(SRC)),
                 os.path.join(here, "..", os.path.basename(SRC)),
                 os.path.join(here, "..", "12_results", os.path.basename(SRC)),
                 os.path.join(here, "..", "47_master", os.path.basename(SRC)),
                 os.path.join(here, "..", "12_results", "MAG_CATALOGUE.tsv")):
        tried.append(os.path.normpath(cand))
        if os.path.exists(cand):
            SRC = cand
            print("   found %s" % os.path.normpath(cand))
            break
    else:
        sys.exit("fig01: cannot find the genome table. Tried:\n         %s\n"
                 "       Pass the path:  python3 %s <table.tsv> <outdir>"
                 % ("\n         ".join(tried), os.path.basename(__file__)))

with open(SRC, encoding="utf-8") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    if SHOW_COLUMNS:
        print("columns in %s:" % SRC)
        for c in (rdr.fieldnames or []):
            print("   %s" % c)
        sys.exit(0)
    cols, missing = resolve(rdr.fieldnames or [])
    if missing:
        sys.exit("fig01: %s has no column for: %s\n"
                 "       columns present:\n         %s\n"
                 "       add the right name to ALIASES at the top of this script."
                 % (SRC, ", ".join(missing),
                    "\n         ".join(rdr.fieldnames or ["<none>"])))
    rows = []
    for r in rdr:
        c, k = fnum(r[cols["comp"]]), fnum(r[cols["cont"]])
        if c is None or k is None:
            continue
        sz = fnum(r[cols["size"]]) if cols["size"] else None
        g = (r.get(cols["genus"]) or "").strip() if cols["genus"] else ""
        tier = (r.get(cols["mimag"]) or "").strip().lower() if cols["mimag"] else ""
        r16 = (r.get(cols["rrna16"]) or "").strip().lower() if cols["rrna16"] else ""
        rows.append({
            "mag": (r[cols["mag"]] or "?").replace("ETO_", ""),
            "genus": "—" if g in ("", "-", "NA", "None") else g,
            "domain": "Archaea" if str(r[cols["domain"]]).lower().startswith("arch")
                      else "Bacteria",
            "comp": c, "cont": k, "size": sz,
            "tier": tier,
            "has16": r16 in ("y", "yes", "true", "1", "t", "+"),
        })

if not rows:
    sys.exit("fig01: no usable rows in %s" % SRC)

# completeness is sometimes stored as a fraction. Detect rather than assume.
if max(r["comp"] for r in rows) <= 1.0:
    for r in rows:
        r["comp"] *= 100.0
        r["cont"] *= 100.0

n_high = sum(1 for r in rows if r["comp"] >= 90 and r["cont"] < 5)
n_med = sum(1 for r in rows if not (r["comp"] >= 90 and r["cont"] < 5)
            and r["comp"] >= 50 and r["cont"] < 10)
n_low = len(rows) - n_high - n_med
print("   %d genomes: %d high, %d medium, %d low  (completeness/contamination only)"
      % (len(rows), n_high, n_med, n_low))

# the table's own tier, which is the one to trust
tiers = sorted({r["tier"] for r in rows if r["tier"]})
n_mimag = sum(1 for r in rows if r["tier"].startswith("h"))
HAVE_TIER = bool(tiers)
if HAVE_TIER:
    print("   MIMAG column values found: %s" % ", ".join(tiers))
    print("   MIMAG high-quality by that column: %d" % n_mimag)
    if n_mimag != n_high:
        miss = [r["mag"] for r in rows
                if r["comp"] >= 90 and r["cont"] < 5 and not r["tier"].startswith("h")]
        print("   -> %d genome(s) clear 90/5 but are NOT MIMAG high: %s"
              % (len(miss), ", ".join(miss)))
        no16 = [m for m in miss
                if not next(r["has16"] for r in rows if r["mag"] == m)]
        if no16:
            print("      of those, %d lack 16S: %s" % (len(no16), ", ".join(no16)))
        print("      That difference IS the rRNA/tRNA criterion. Both numbers are")
        print("      correct; the panel shows both.")
if n_low:
    print("   !! %d genome(s) fall outside MIMAG medium." % n_low)
print("   NOTE: these are the completeness/contamination thresholds ONLY.")
print("         MIMAG high-quality additionally requires 23S + 16S + 5S rRNA and")
print("         >= 18 tRNAs, so the high-quality count in the text is expected to")
print("         be LOWER than %d. If the text says 8, the four-genome difference is"
      % n_high)
print("         the rRNA/tRNA criterion and should be stated, not reconciled away.")

sizes = [r["size"] for r in rows if r["size"]]
smin, smax = (min(sizes), max(sizes)) if sizes else (1.0, 1.0)


def area(r):
    if not r["size"] or smax == smin:
        return 46.0
    frac = (r["size"] - smin) / (smax - smin)
    return 24.0 + frac * 150.0          # AREA scales linearly with genome size


# ---------------------------------------------------------------- figure
S.apply()
W = FIG_W_MM / 25.4
fig = plt.figure(figsize=(W, W * 0.80))
ax = fig.add_axes([0.098, 0.345, 0.870, 0.555])
S.recessive_grid(ax, axis="both")
S.despine(ax, keep=("left", "bottom"))

xlo = max(45.0, min(r["comp"] for r in rows) - 4)
yhi = max(6.0, max(r["cont"] for r in rows) * 1.35)

ax.add_patch(Rectangle((90, 0), 100 - 90, 5, facecolor=S.DOMAIN["Bacteria"],
                       alpha=0.055, zorder=0, linewidth=0))
ax.axvline(90, color=S.RULE, lw=0.8, ls=(0, (4, 2)), zorder=1)
ax.axhline(5, color=S.RULE, lw=0.8, ls=(0, (4, 2)), zorder=1)
if yhi > 10:
    ax.axhline(10, color=S.RULE, lw=0.8, ls=(0, (1, 2)), zorder=1)

for r in rows:
    full = HAVE_TIER and r["tier"].startswith("h")
    ax.scatter(r["comp"], r["cont"], s=area(r), color=S.DOMAIN[r["domain"]],
               alpha=0.85, zorder=5 if full else 4,
               linewidths=1.1 if full else 0,
               edgecolors=S.INK if full else "none")

# DIRECT LABELS, DE-COLLIDED.
# The high-quality genomes cluster in one corner by definition, so labels placed
# at a fixed offset land on top of each other -- the first render stacked eight
# of them into an unreadable block. Offsets are solved in DISPLAY space (points
# are what overlap, not data units) by pushing overlapping labels apart, and any
# label pushed far from its dot gets a leader line so the pairing stays readable.
# Label the box and nothing else. The first run on real data labelled eighteen
# of twenty-nine points -- every genome within 1,5 points of a boundary as well --
# and the near-boundary ones carry no argument. A genome at 89,5 % completeness
# is not a story; the twelve inside the box are.
lab = []
for r in rows:
    inbox = r["comp"] >= 90 and r["cont"] < 5
    called_out = r["mag"] in HILITE_MAGS
    if called_out:
        lab.append({"r": r, "edge": False, "inbox": inbox})
if not HILITE_MAGS:
    print("   no genomes named on the panel (HILITE_MAGS is empty). Add the ones")
    print("   the text argues about to HILITE_MAGS at the top of this script.")
elif len(lab) > 6:
    print("   !! %d labels in one corner will crowd. Keep HILITE_MAGS to a few."
          % len(lab))

if lab:
    # THE UNITS BIT, WHICH I GOT WRONG THE FIRST TIME. ax.transData gives DISPLAY
    # PIXELS at fig.dpi; font sizes are POINTS. A separation of "11" in pixels at
    # dpi 100 is 7.9 pt -- less than the line height of 7.5 pt text -- so the
    # labels were pushed apart by less than their own height and still overlapped.
    # Everything below is computed in points and converted once.
    fig.canvas.draw()
    P2PX = fig.dpi / 72.0
    LIFT = 9.0 * P2PX
    for L in lab:
        ax_, ay = ax.transData.transform((L["r"]["comp"], L["r"]["cont"]))
        L["ax"], L["ay"] = ax_, ay
        L["w"] = len(L["r"]["mag"]) * 0.54 * S.PT_VALUE * P2PX
        L["h"] = 1.35 * S.PT_VALUE * P2PX
        L["x"], L["y"] = ax_, ay + LIFT

    # THE DOTS ARE OBSTACLES TOO. Separating labels from each other is only half
    # the job -- in the first fix they stopped overlapping one another and sat on
    # top of the markers instead, which is just as unreadable. Every point is a
    # fixed circle the labels must stay outside of.
    obs = []
    for r in rows:
        ox_, oy_ = ax.transData.transform((r["comp"], r["cont"]))
        rad = (area(r) / 3.14159) ** 0.5 * P2PX + 1.5 * P2PX
        obs.append((ox_, oy_, rad))

    for _ in range(800):
        moved = False
        for i in range(len(lab)):
            for j in range(i + 1, len(lab)):
                a, b = lab[i], lab[j]
                ox = (a["w"] + b["w"]) / 2 - abs(a["x"] - b["x"])
                oy = (a["h"] + b["h"]) / 2 - abs(a["y"] - b["y"])
                if ox <= 0 or oy <= 0:
                    continue
                # separate along whichever axis needs the smaller move
                if oy <= ox * 0.30:
                    d = (oy / 2 + 0.4) * (1 if a["y"] <= b["y"] else -1)
                    a["y"] -= d; b["y"] += d
                else:
                    d = (ox / 2 + 0.4) * (1 if a["x"] <= b["x"] else -1)
                    a["x"] -= d; b["x"] += d
                moved = True
        for L in lab:                      # clear of every marker
            for ox_, oy_, rad in obs:
                dx_ = L["x"] - ox_
                dy_ = L["y"] - oy_
                need_y = rad + L["h"] / 2
                need_x = rad + L["w"] / 2
                if abs(dx_) < need_x and abs(dy_) < need_y:
                    push = need_y - abs(dy_) + 0.4
                    L["y"] += push if dy_ >= 0 else -push
                    moved = True
        for L in lab:                      # weak spring home, so nothing drifts off
            L["x"] += (L["ax"] - L["x"]) * 0.03
            L["y"] += (L["ay"] + LIFT - L["y"]) * 0.03
        if not moved:
            break

    # keep labels inside the axes DURING the solve, not after it. Clamping at the
    # end can shove two already-separated labels back on top of each other, which
    # is exactly what left one pair superimposed.
    x0, y0 = ax.transData.transform((xlo, -yhi * 0.045))
    x1, y1 = ax.transData.transform((101.5, yhi))

    def clamp():
        for L in lab:
            L["x"] = min(max(L["x"], x0 + L["w"] / 2), x1 - L["w"] / 2)
            L["y"] = min(max(L["y"], y0 + L["h"] / 2), y1 - L["h"] / 2)

    clamp()

    # FINISHING PASS, NO SPRING. The spring keeps labels near their dots, but it
    # also pulls them back into overlaps that the same iteration had just fixed --
    # which is why two pairs were still touching after the first loop. Run the
    # separation alone at the end so the last thing that happens is un-overlapping.
    for _ in range(1500):
        moved = False
        for i in range(len(lab)):
            for j in range(i + 1, len(lab)):
                a, b = lab[i], lab[j]
                ox = (a["w"] + b["w"]) / 2 - abs(a["x"] - b["x"])
                oy = (a["h"] + b["h"]) / 2 - abs(a["y"] - b["y"])
                if ox <= 0 or oy <= 0:
                    continue
                d = (oy / 2 + 0.5) * (1 if a["y"] <= b["y"] else -1)
                a["y"] -= d; b["y"] += d
                moved = True
        for L in lab:
            for ox_, oy_, rad in obs:
                dy_ = L["y"] - oy_
                if abs(L["x"] - ox_) < rad + L["w"] / 2 and abs(dy_) < rad + L["h"] / 2:
                    push = rad + L["h"] / 2 - abs(dy_) + 0.5
                    L["y"] += push if dy_ >= 0 else -push
                    moved = True
        clamp()
        if not moved:
            break

    for L in lab:
        dx = (L["x"] - L["ax"]) / P2PX
        dy = (L["y"] - L["ay"]) / P2PX
        far = (dx * dx + dy * dy) ** 0.5 > 11.0
        ax.annotate(L["r"]["mag"], (L["r"]["comp"], L["r"]["cont"]),
                    textcoords="offset points", xytext=(dx, dy),
                    ha="center", va="center", fontsize=S.PT_VALUE,
                    color=S.HILITE if L["edge"] else S.INK, zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.55, shrinkA=1, shrinkB=3,
                                    color=S.INK_MUTED, alpha=0.75) if far else None)

ax.set_xlim(xlo, 101.5)
ax.set_ylim(-yhi * 0.045, yhi)
ax.set_xlabel(EL["xlabel"])
ax.set_ylabel(EL["ylabel"])
# set_ticks BEFORE set_ticklabels -- matplotlib warns otherwise because a
# locator can hand back a different tick count than the labels you just built.
xt = [t for t in ax.get_xticks() if xlo <= t <= 101.5]
yt = [t for t in ax.get_yticks() if -yhi * 0.045 <= t <= yhi]
ax.set_xticks(xt); ax.set_xticklabels([S.num(t, 0) for t in xt])
ax.set_yticks(yt); ax.set_yticklabels([S.num(t, 0) for t in yt])
ax.set_title(EL["title"], loc="left", pad=8)

# the tier note sits LEFT of the 90 line: to its right is the shaded band and,
# beyond the axes, the legend -- the first render put all three on top of each other
ax.text(101.0, yhi * 0.965, EL["high"], ha="right", va="top",
        fontsize=S.PT_VALUE, color=S.INK_2, linespacing=1.4, zorder=6)

# legend: domain, then a size key built from the real minimum and maximum
dom = [Patch(facecolor=S.DOMAIN["Bacteria"], label=S.DOMAIN_EL["Bacteria"]),
       Patch(facecolor=S.DOMAIN["Archaea"], label=S.DOMAIN_EL["Archaea"])]
if HAVE_TIER:
    dom.append(Line2D([], [], marker="o", linestyle="none", markersize=6,
                      markerfacecolor=S.INK_MUTED, markeredgecolor=S.INK,
                      markeredgewidth=1.1, label="MIMAG high-quality"))
leg1 = ax.legend(handles=dom, loc="upper left", bbox_to_anchor=(0.012, 0.995),
                 handlelength=1.0, handleheight=1.0, borderpad=0.2)
ax.add_artist(leg1)
if sizes:
    keys = [smin, (smin + smax) / 2, smax]
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=(24.0 + f * 150.0) ** 0.5,
                      markerfacecolor=S.INK_MUTED, markeredgewidth=0,
                      label="%s Mb" % S.num(k / 1e6, 1))
               for k, f in zip(keys, (0.0, 0.5, 1.0))]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.012, 0.995 - (0.075 * len(dom))),
              title=EL["size"], title_fontsize=S.PT_VALUE, labelspacing=1.05,
              handletextpad=1.0, borderpad=0.3)
else:
    print("   !! no genome-size column found, so every marker is the same size and")
    print("      the size key is omitted. Add the column name to ALIASES['size'].")

S.caption(fig, EL["caption"] % (S.num(n_high, 0), S.num(n_med, 0))
          + ("" if not HAVE_TIER else
             " Τα γονιδιώματα με περίγραμμα (%s) πληρούν το πρότυπο MIMAG υψηλής "
             "ποιότητας στο σύνολό του." % S.num(n_mimag, 0)))

# ---- measure what was actually drawn, in mm, before saving
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))
if bb.width * 25.4 > FIG_W_MM + 1:
    print("   !! wider than the text block. Word will scale it down on insert and")
    print("      every font in it shrinks with it. Reduce FIG_W_MM or the legend.")

S.save(fig, "S_MAG_quality", OUT)
