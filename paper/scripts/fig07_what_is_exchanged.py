#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig07_what_is_exchanged.py -- τι ανταλλάσσεται στην ανακατασκευή.

    python3 fig07_what_is_exchanged.py [67_network/compound_summary.tsv] [outdir]

THE ARGUMENT THIS PANEL EXISTS TO SETTLE. The chapter promises syntrophy: a
fermenting community handing hydrogen, formate and acetate to methanogens, with
sulfate reducers competing for the same currency. The reconstruction says
something else. Amino acids carry 61,9% of the total exchange; H2, formate,
acetate, methanol and the C1 compounds together carry 2,2%. HYDROGEN ALONE
CARRIES 0,2%.

That is not a weaker result. A sulfate-replete lagoon in which the model predicts
no H2 syntrophy says the same thing as every other line in this thesis: the
sulfate reducers take the hydrogen. But the text has to say THIS, and a reader
will not believe it from a sentence. They will believe a bar 15,55 long beside a
bar 0,30 long.

FORM. Ranked horizontal bars, linear axis. Linear, not log, precisely because the
contrast is the finding -- a log axis would flatter hydrogen into visibility. The
four substrates the chapter argues about are pulled out below a rule and drawn to
the same scale even though they rank 23rd, 42nd, 47th and 55th, because "not in
the top twelve" is a fact about my cutoff and "0,30 against 15,55" is a fact
about the data.

ENCODED ONCE, EACH:
    bar length   summed SMETANA score across every pair that moved it
    colour       three classes only -- amino acid, energetic/C1, everything else.
                 Colour carries the argument, and the argument is two-sided
    the rule     separates the ranking from the named substrates below it
"""
import sys, csv, os
import etoliko_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SRC = args[0] if args else "67_network/compound_summary.tsv"
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 160.0
TOPN = 12

# The four the chapter argues about, in the order it argues about them.
ARGUED = ["cpd11640", "cpd00047", "cpd00029", "cpd00116"]

# 🔴 THE CLASSIFIER IN STEP 67 MISSED FOUR, AND THE ONLY REASON I KNOW IS THAT IT
# PRINTED ALL 64 RATHER THAN JUST THE CLASS TOTALS. PAN is pantothenate and mql7
# is menaquinol-7 -- both vitamins, both filed under "other" because the rule
# matched on names. Formaldehyde and CO are C1 compounds. Corrected here, by
# accession, in the open.
RECLASS = {
    "cpd00644": "vitamin/cofactor",   # PAN = pantothenate
    "cpd11451": "vitamin/cofactor",   # mql7 = menaquinol-7
    "cpd00055": "energy/C1",          # formaldehyde
    "cpd00204": "energy/C1",          # carbon monoxide
}

# three groups, because the claim is two-sided and everything else is context
GROUP = {"amino acid": "aa", "energy/C1": "c1"}
COL = {"aa": S.DOMAIN["Bacteria"], "c1": S.DOMAIN["Archaea"], "other": S.INK_MUTED}

EL = {
    "title": "Τι ανταλλάσσεται: αμινοξέα, όχι υδρογόνο",
    "xlabel": "άθροισμα βαθμών SMETANA σε όλα τα ζεύγη",
    "aa": "αμινοξέα",
    "c1": "ενεργειακά και C1",
    "other": "άλλα (νουκλεοσίδια, οργανικά οξέα, βιταμίνες)",
    "argued": "τα υποστρώματα που συζητά το κεφάλαιο",
    "caption": (
        "Οι δώδεκα ενώσεις με το μεγαλύτερο άθροισμα βαθμών SMETANA, και από κάτω "
        "— χωριστά, στην ίδια κλίμακα — τα τέσσερα υποστρώματα γύρω από τα οποία "
        "χτίζεται το κεφάλαιο. Τα αμινοξέα μεταφέρουν %s%% του συνόλου· το "
        "υδρογόνο, το μυρμηκικό, το οξικό και η μεθανόλη μαζί %s%%, και το "
        "υδρογόνο μόνο του %s%%. Ο άξονας είναι γραμμικός επίτηδες: σε "
        "λογαριθμική κλίμακα η αντίθεση, που είναι το ίδιο το αποτέλεσμα, θα "
        "εξαφανιζόταν. Πρόκειται για δυνατότητα και όχι για δραστηριότητα — "
        "μοντέλα με συμπλήρωση κενών πάνω σε ορισμένο θρεπτικό μέσο (σχήμα 6)."),
}

# Greek names for what is actually moving. Anything not here prints as it comes
# out of the SBML and is reported, so the gap is visible rather than silent.
GR = {
    "L-Cysteine": "L-κυστεΐνη", "L-Methionine": "L-μεθειονίνη",
    "L-Arginine": "L-αργινίνη", "L-Aspartate": "L-ασπαρτικό",
    "L-Valine": "L-βαλίνη", "L-Tryptophan": "L-τρυπτοφάνη",
    "L-Glutamate": "L-γλουταμικό", "L-Alanine": "L-αλανίνη",
    "L-Glutamine": "L-γλουταμίνη", "L-Threonine": "L-θρεονίνη",
    "L-Serine": "L-σερίνη", "L-Leucine": "L-λευκίνη",
    "L-Isoleucine": "L-ισολευκίνη", "L-Proline": "L-προλίνη",
    "L-Lysine": "L-λυσίνη", "L-Histidine": "L-ιστιδίνη",
    "L-Tyrosine": "L-τυροσίνη", "L-Phenylalanine": "L-φαινυλαλανίνη",
    "L-Asparagine": "L-ασπαραγίνη", "Glycine": "γλυκίνη",
    "Ornithine": "ορνιθίνη", "D-Serine": "D-σερίνη",
    "PPi": "πυροφωσφορικό", "Formaldehyde": "φορμαλδεΰδη",
    "indol": "ινδόλιο", "Ethanol": "αιθανόλη", "Pyruvate": "πυροσταφυλικό",
    "Uridine": "ουριδίνη", "D-Ribose": "D-ριβόζη", "D-Lactate": "D-γαλακτικό",
    "L-Lactate": "L-γαλακτικό", "Citrate": "κιτρικό", "PAN": "παντοθενικό",
    "Acetaldehyde": "ακεταλδεΰδη", "L-Malate": "L-μηλικό", "NAD": "NAD",
    "D-Fructose": "D-φρουκτόζη", "CoA": "συνένζυμο Α", "Formate": "μυρμηκικό",
    "Acetate": "οξικό", "H2": "υδρογόνο (H₂)", "Methanol": "μεθανόλη",
    "Succinate": "ηλεκτρικό", "Fumarate": "φουμαρικό", "CO": "CO",
    "Glycerol": "γλυκερόλη", "Putrescine": "πουτρεσκίνη",
}


def find(p, extra=()):
    if os.path.exists(p):
        return p
    here = os.path.dirname(os.path.abspath(__file__))
    b = os.path.basename(p)
    for c in [os.path.join(here, b), os.path.join(here, "..", b),
              os.path.join(here, "..", os.path.dirname(p), b)] + \
             [os.path.join(here, "..", d, b) for d in extra]:
        if os.path.exists(c):
            print("   found %s" % os.path.normpath(c))
            return c
    return None


SRC = find(SRC, ("67_network",)) or sys.exit(
    "fig07: cannot find compound_summary.tsv -- pass it as the first argument")


def f(v, d=0.0):
    try: return float(str(v).replace(",", ".").strip())
    except (TypeError, ValueError): return d


rows = []
with open(SRC, encoding="utf-8") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    need = ("cpd", "name", "class", "sum_smetana")
    miss = [c for c in need if c not in (rdr.fieldnames or [])]
    if miss:
        sys.exit("fig07: %s lacks %s\n       has: %s"
                 % (SRC, ", ".join(miss), ", ".join(rdr.fieldnames or [])))
    for r in rdr:
        cpd = r["cpd"].strip()
        nm = (r["name"] or "").strip()
        nm = nm.rsplit("-c0", 1)[0].rsplit("-e0", 1)[0].strip() or cpd
        kls = RECLASS.get(cpd, (r["class"] or "other").strip())
        rows.append({"cpd": cpd, "name": nm, "cls": kls,
                     "grp": GROUP.get(kls, "other"),
                     "v": f(r["sum_smetana"]), "n": int(f(r.get("n_transfers")))})
if not rows:
    sys.exit("fig07: no rows in %s" % SRC)

for c, k in RECLASS.items():
    hit = next((r for r in rows if r["cpd"] == c), None)
    if hit:
        print("   reclassified %s (%s) -> %s" % (c, hit["name"], k))

total = sum(r["v"] for r in rows)
by = {}
for r in rows:
    by[r["cls"]] = by.get(r["cls"], 0.0) + r["v"]
print("   %d compound(s), total score %.3f" % (len(rows), total))
for k, v in sorted(by.items(), key=lambda kv: -kv[1]):
    print("     %-18s %8.3f   %5.1f%%" % (k, v, 100.0 * v / total))

idx = {r["cpd"]: r for r in rows}
aa = 100.0 * by.get("amino acid", 0.0) / total
arg4 = sum(idx[c]["v"] for c in ARGUED if c in idx)
h2 = idx.get("cpd11640", {}).get("v", 0.0)
print("   amino acids %.1f%%;  H2+formate+acetate+methanol %.3f (%.1f%%);"
      "  H2 alone %.3f (%.2f%%)"
      % (aa, arg4, 100.0 * arg4 / total, h2, 100.0 * h2 / total))

rows.sort(key=lambda r: -r["v"])
top = [r for r in rows if r["cpd"] not in ARGUED][:TOPN]
low = sorted((idx[c] for c in ARGUED if c in idx), key=lambda r: -r["v"])
missing = [c for c in ARGUED if c not in idx]
if missing:
    print("   !! not present in the table at all: %s" % ", ".join(missing))

untranslated = sorted({r["name"] for r in top + low if r["name"] not in GR})
if untranslated:
    print("   no Greek name for: %s" % ", ".join(untranslated))


def label(r):
    return GR.get(r["name"], r["name"])


# ---------------------------------------------------------------- layout
S.apply()
W = FIG_W_MM / 25.4
ROW = 0.185                       # inches per bar
GAP = 0.34                        # the rule between the two blocks
BARS = ROW * (len(top) + len(low)) + GAP
TOP_IN, BOT_IN = 0.62, 1.24
H = TOP_IN + BARS + BOT_IN
fig = plt.figure(figsize=(W, H))
LEFT = 1.62
ax = fig.add_axes([LEFT / W, BOT_IN / H, (W - LEFT - 0.62) / W, BARS / H])
S.recessive_grid(ax, axis="x")
S.despine(ax, keep=("bottom",))

VMAX = max(r["v"] for r in top + low)
order = top + [None] + low       # None marks where the rule goes
ypos, y = {}, len(top) + len(low)
for r in order:
    if r is None:
        y -= 0.0
        continue
    y -= 1
    ypos[r["cpd"]] = y + 0.5 - (0 if r in top else 0.5)

# draw with an explicit gap: top block sits above the rule, argued block below
ys = {}
for i, r in enumerate(top):
    ys[r["cpd"]] = len(top) + len(low) + 0.55 - i
for i, r in enumerate(low):
    ys[r["cpd"]] = len(low) - i - 0.45

for r in top + low:
    yy = ys[r["cpd"]]
    ax.barh(yy, r["v"], height=0.66, color=COL[r["grp"]], linewidth=0, zorder=3)
    ax.text(r["v"] + VMAX * 0.012, yy, S.num(r["v"], 2), va="center", ha="left",
            fontsize=S.PT_VALUE, color=S.INK_2, zorder=4)
    ax.text(-VMAX * 0.015, yy, label(r), va="center", ha="right",
            fontsize=S.PT_TICK, color=S.INK_2, zorder=4)

ax.set_ylim(-0.9, len(top) + len(low) + 1.1)
ax.set_xlim(0, VMAX * 1.14)
ax.set_yticks([])
ax.set_xlabel(EL["xlabel"])
_tk = [t for t in ax.get_xticks() if 0 <= t <= VMAX * 1.14]
ax.set_xticks(_tk)                       # fix the ticks before relabelling them,
ax.set_xticklabels([S.num(t, 0) for t in _tk])   # or matplotlib warns and drifts

# the rule, and what is below it
yr = len(low) + 0.30
ax.axhline(yr, color=S.RULE, lw=0.8, zorder=2)
ax.text(VMAX * 1.13, yr - 0.20, EL["argued"], ha="right", va="top",
        fontsize=S.PT_VALUE, color=S.INK_MUTED, zorder=4)

fig.text(0.012, 1 - 0.30 / H, EL["title"], ha="left", va="top",
         fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)

ax.legend(handles=[Patch(facecolor=COL["aa"], label=EL["aa"]),
                   Patch(facecolor=COL["c1"], label=EL["c1"]),
                   Patch(facecolor=COL["other"], label=EL["other"])],
          loc="lower right", bbox_to_anchor=(0.998, 0.012),
          handlelength=1.1, handleheight=1.1, borderpad=0.35,
          frameon=True, facecolor=S.SURFACE, edgecolor=S.RULE, framealpha=1.0)


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
caption_fit(fig, EL["caption"] % (S.num(aa, 1), S.num(100.0 * arg4 / total, 1),
                                  S.num(100.0 * h2 / total, 2)))

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))

S.save(fig, "fig07_what_is_exchanged", OUT)
