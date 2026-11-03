#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig05_sed_vs_liq.py -- τι διάλεξε η καλλιέργεια.

    python3 fig05_sed_vs_liq.py [MASTER.tsv] [outdir]

FORM. Two abundances for the same 29 genomes, and the question is not how big
either is but how they DIFFER -- so both go on axes and the 1:1 line does the
work. A genome above the line was enriched by the medium; one below it is
relatively commoner in the undisturbed sediment than in the bottles. A paired bar
chart would show the same numbers and hide the comparison.

LOG AXES, BECAUSE THE DATA SPANS ORDERS OF MAGNITUDE. On a linear axis the three
dominant culture organisms would compress everything else onto the origin. The
consequence is that ZERO CANNOT BE PLOTTED: a genome absent from one sample is
drawn on the margin at the detection floor and marked open, not dropped, because
"not detected" is a value and dropping it would quietly shrink the catalogue.

ENCODED ONCE, EACH:
    x, y         relative abundance in sediment and in cultures
    the diagonal equality -- the only reference the figure needs
    colour       DOMAIN only
    open marker  below the floor in one sample, i.e. not detected there
    ring         in the in-situ ten, if MASTER.tsv says so
"""
import sys, csv, os, math
import etoliko_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SRC = args[0] if args else "MASTER.tsv"
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 160.0

# genomes the text argues about; everything else stays unlabelled
HILITE_MAGS = ("LIQ_MAG_003", "LIQ_MAG_011", "LIQ_MAG_022")

EL = {
    "title":  "Σχετική αφθονία: ίζημα έναντι καλλιέργειας",
    "xlabel": "σχετική αφθονία στο ίζημα (%)",
    "ylabel": "σχετική αφθονία στην καλλιέργεια (%)",
    "eq":     "ίση αφθονία",
    "enr":    "εμπλουτισμός\nστην καλλιέργεια",
    "dep":    "σχετικά συχνότερα\nστο ίζημα",
    "nd":     "μη ανιχνεύσιμο σε ένα δείγμα",
    "ten":    "σύνολο in situ",
    # The caption is assembled from these pieces at the end, because two of them
    # describe things that are only on the panel FOR SOME DATA. The first render
    # carried a fixed caption explaining open circles on the margin while the
    # legend, correctly, drew no such entry -- the figure had none. A caption
    # that describes a marker the reader cannot find is the same class of error
    # as a legend entry with nothing behind it.
    "cap_base": ("Κάθε κύκλος είναι ένα MAG· οι άξονες δίνουν τη σχετική αφθονία στα δύο "
                 "δείγματα, σε λογαριθμική κλίμακα. Η διαγώνιος είναι η ισότητα: πάνω "
                 "από αυτήν το γονιδίωμα είναι σχετικά αφθονότερο στην καλλιέργεια απ' ό,τι "
                 "στο ίζημα, κάτω από αυτήν συμβαίνει το αντίστροφο. Το χρώμα δηλώνει "
                 "επικράτεια και ο μαύρος δακτύλιος σημειώνει τα γονιδιώματα του συνόλου "
                 "in situ."),
    # drawn only when at least one genome falls below the detection floor
    "cap_nd":   ("Οι ανοιχτοί κύκλοι στο περιθώριο δεν ανιχνεύθηκαν στο αντίστοιχο δείγμα "
                 "και σχεδιάζονται στο όριο ανίχνευσης· δεν είναι μηδέν, είναι «κάτω από "
                 "το όριο»."),
    # drawn only when the lower triangle is empty, so that the empty half of the
    # panel reads as the result and not as wasted space
    "cap_allup": ("Και τα %d γονιδιώματα βρίσκονται πάνω από τη διαγώνιο. Το κενό κάτω "
                  "ημιεπίπεδο δεν είναι αδιάθετος χώρος: είναι το αποτέλεσμα, ότι δηλαδή "
                  "κανένα γονιδίωμα του καταλόγου δεν είναι σχετικά συχνότερο στο "
                  "αδιατάρακτο ίζημα απ' ό,τι στις καλλιέργειες. Η σχετική αφθονία είναι "
                  "κλάσμα του συνόλου, οπότε αυτό σημαίνει ανακατανομή του πληθυσμού υπέρ "
                  "των οργανισμών που το θρεπτικό μέσο ευνοεί, και όχι ότι κάθε οργανισμός "
                  "πολλαπλασιάστηκε."),
    "cap_mixed": ("Πάνω από τη διαγώνιο βρίσκονται %d γονιδιώματα και κάτω από αυτήν %d."),
}

if not os.path.exists(SRC):
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.path.join(here, os.path.basename(SRC)),
                 os.path.join(here, "..", os.path.basename(SRC)),
                 os.path.join(here, "..", "47_master", os.path.basename(SRC))):
        if os.path.exists(cand):
            SRC = cand; print("   found %s" % os.path.normpath(cand)); break
    else:
        sys.exit("fig05: cannot find %s -- pass the path as the first argument" % SRC)

SEDC = ("sed_relabund", "sed_rel_abund", "sed_abundance", "sed_abund")
LIQC = ("liq_relabund", "liq_rel_abund", "liq_abundance", "liq_abund")


def pick(names, have):
    low = {h.lower(): h for h in have}
    return next((low[n] for n in names if n in low), None)


def fnum(v):
    try: return float(str(v).replace(",", ".").replace("%", "").strip())
    except (TypeError, ValueError): return None


with open(SRC, encoding="utf-8") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    cols = rdr.fieldnames or []
    cs, cl = pick(SEDC, cols), pick(LIQC, cols)
    if not cs or not cl:
        sys.exit("fig05: need a sediment and a culture abundance column.\n"
                 "       columns present:\n         %s" % "\n         ".join(cols))
    rows = []
    for r in rdr:
        a, b = fnum(r.get(cs)), fnum(r.get(cl))
        if a is None and b is None:
            continue
        rows.append({
            "mag": (r.get("mag_id") or r.get("working_name") or "?").replace("ETO_", ""),
            "genus": (r.get("genus") or "").strip() or "—",
            "domain": "Archaea" if str(r.get("domain", "")).lower().startswith("arch")
                      else "Bacteria",
            "ten": str(r.get("in_situ_ten", "")).lower().startswith(("y", "t", "1")),
            "sed": a or 0.0, "liq": b or 0.0,
        })
if not rows:
    sys.exit("fig05: no usable rows in %s" % SRC)

pos = [v for r in rows for v in (r["sed"], r["liq"]) if v > 0]
if not pos:
    sys.exit("fig05: every abundance is zero")
FLOOR = min(pos) / 3.0           # margin position for "not detected"
hi = max(pos) * 2.2
print("   %d genomes; detected range %.4g – %.4g %%; floor drawn at %.4g"
      % (len(rows), min(pos), max(pos), FLOOR))
nd = sum(1 for r in rows if r["sed"] <= 0 or r["liq"] <= 0)
if nd:
    print("   %d genome(s) not detected in one sample, drawn open on the margin" % nd)

S.apply()
W = FIG_W_MM / 25.4
fig = plt.figure(figsize=(W, W * 0.82))
ax = fig.add_axes([0.115, 0.255, 0.855, 0.665])
ax.set_xscale("log"); ax.set_yscale("log")
S.recessive_grid(ax, axis="both")
S.despine(ax, keep=("left", "bottom"))

# with nothing on the margin there is no need to reserve room for it, and the
# axes must stay equal or the 45-degree diagonal stops meaning equality
lo = (FLOOR / 1.6) if nd else (min(pos) / 2.0)
ax.plot([lo, hi], [lo, hi], color=S.RULE, lw=0.9, ls=(0, (4, 2)), zorder=1)
ax.text(hi, hi, EL["eq"], ha="right", va="bottom", fontsize=S.PT_VALUE,
        color=S.INK_MUTED, rotation=45, rotation_mode="anchor", zorder=1)
# The "enriched" / "commoner in sediment" corner labels are deliberately NOT
# drawn. Where the two triangles are empty depends entirely on the data, and in
# the first render one landed on the not-detected column and the other on the
# legend. The diagonal is labelled and the caption says which side is which;
# a label that sometimes sits on the data is worse than no label.

for r in rows:
    x = r["sed"] if r["sed"] > 0 else FLOOR
    y = r["liq"] if r["liq"] > 0 else FLOOR
    off = (r["sed"] <= 0 or r["liq"] <= 0)
    ax.scatter(x, y, s=52, zorder=5 if r["ten"] else 4,
               facecolors="none" if off else S.DOMAIN[r["domain"]],
               edgecolors=S.DOMAIN[r["domain"]] if off else
                          (S.INK if r["ten"] else "none"),
               linewidths=1.1 if (off or r["ten"]) else 0, alpha=0.9)

ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel(EL["xlabel"]); ax.set_ylabel(EL["ylabel"])
ax.set_title(EL["title"], loc="left", pad=8)

h = [Patch(facecolor=S.DOMAIN["Bacteria"], label=S.DOMAIN_EL["Bacteria"]),
     Patch(facecolor=S.DOMAIN["Archaea"], label=S.DOMAIN_EL["Archaea"])]
if any(r["ten"] for r in rows):
    h.append(Line2D([], [], marker="o", linestyle="none", markersize=6,
                    markerfacecolor=S.INK_MUTED, markeredgecolor=S.INK,
                    markeredgewidth=1.1, label=EL["ten"]))
if nd:
    h.append(Line2D([], [], marker="o", linestyle="none", markersize=6,
                    markerfacecolor="none", markeredgecolor=S.INK_MUTED,
                    markeredgewidth=1.1, label=EL["nd"]))
ax.legend(handles=h, loc="lower right", bbox_to_anchor=(0.995, 0.012),
          handlelength=1.0, handleheight=1.0, borderpad=0.3)

# --- direct labels, measured and de-collided --------------------------------
# A fixed offset above the marker is only safe on a sparse panel. Here the
# genomes the text argues about sit inside the crowd, and in the first render
# LIQ_MAG_011's name landed on a neighbouring point. Each label is now tried in
# eight directions at three distances; the position whose RENDERED box overlaps
# the fewest markers, no other label and no part of the legend is kept, and a
# hairline leader is drawn whenever the winner had to move far from its marker.
lab = [r for r in rows if r["mag"] in HILITE_MAGS]
if not HILITE_MAGS:
    print("   no genomes named on the panel; add them to HILITE_MAGS")

if lab:
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    PPX = fig.dpi / 72.0
    MRAD = (0.5 * math.sqrt(52.0) + 1.4) * PPX      # marker radius, display px

    marker_px = []
    for r in rows:
        mx = r["sed"] if r["sed"] > 0 else FLOOR
        my = r["liq"] if r["liq"] > 0 else FLOOR
        marker_px.append(ax.transData.transform((mx, my)))

    keepout = [ax.get_legend().get_window_extent(rend)]
    axbb = ax.get_window_extent(rend)

    CAND = [(d * math.cos(math.radians(a)), d * math.sin(math.radians(a)))
            for d in (10.0, 16.0, 23.0) for a in range(0, 360, 45)]

    def score(box, dist):
        """lower is better"""
        pen = 0.0
        for px, py in marker_px:                     # marker centre inside the box
            if (box.x0 - MRAD <= px <= box.x1 + MRAD and
                    box.y0 - MRAD <= py <= box.y1 + MRAD):
                pen += 100.0
        for other in keepout:
            if box.overlaps(other):
                pen += 400.0
        if not (axbb.x0 <= box.x0 and box.x1 <= axbb.x1 and
                axbb.y0 <= box.y0 and box.y1 <= axbb.y1):
            pen += 800.0                             # never let a name leave the axes
        return pen + dist                            # ties go to the closest offset

    for r in lab:
        x = r["sed"] if r["sed"] > 0 else FLOOR
        y = r["liq"] if r["liq"] > 0 else FLOOR
        best = None
        for dx, dy in CAND:
            t = ax.annotate(r["mag"], (x, y), textcoords="offset points",
                            xytext=(dx, dy),
                            ha="left" if dx > 1 else ("right" if dx < -1 else "center"),
                            va="bottom" if dy > 1 else ("top" if dy < -1 else "center"),
                            fontsize=S.PT_VALUE, color=S.INK, zorder=7)
            box = t.get_window_extent(rend)
            s = score(box, math.hypot(dx, dy))
            t.remove()
            if best is None or s < best[0]:
                best = (s, dx, dy, box)
        s, dx, dy, box = best
        # A leader is drawn ALWAYS, not only when the offset is large. The label
        # is pushed off its marker precisely because other markers are near, so
        # the moment it moves at all there is more than one point it could be
        # naming; the hairline is what removes the ambiguity.
        ax.annotate(r["mag"], (x, y), textcoords="offset points", xytext=(dx, dy),
                    ha="left" if dx > 1 else ("right" if dx < -1 else "center"),
                    va="bottom" if dy > 1 else ("top" if dy < -1 else "center"),
                    fontsize=S.PT_VALUE, color=S.INK, zorder=7,
                    arrowprops=dict(arrowstyle="-", color=S.INK_MUTED, lw=0.6,
                                    shrinkA=1.0, shrinkB=4.5,
                                    connectionstyle="arc3,rad=0"))
        keepout.append(box)
        print("   label %-14s offset (%+5.1f, %+5.1f) pt, penalty %.0f"
              % (r["mag"], dx, dy, s))
        if s >= 100.0:
            print("   ! %s could not be placed clear of every marker; "
                  "check the render" % r["mag"])

above = sum(1 for r in rows if r["liq"] > r["sed"] > 0)
below = sum(1 for r in rows if 0 < r["liq"] < r["sed"])
print("   above the diagonal (enriched in culture): %d;  below: %d" % (above, below))
if below == 0:
    print("   every genome is relatively commoner in the cultures. That empty lower")
    print("   triangle is the result, not wasted space -- say so in the text.")

cap = [EL["cap_base"]]
if nd:
    cap.append(EL["cap_nd"])
if below == 0 and above == len(rows):
    cap.append(EL["cap_allup"] % above)
else:
    cap.append(EL["cap_mixed"] % (above, below))
S.caption(fig, " ".join(cap))

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))

S.save(fig, "fig05_sed_vs_liq", OUT)
