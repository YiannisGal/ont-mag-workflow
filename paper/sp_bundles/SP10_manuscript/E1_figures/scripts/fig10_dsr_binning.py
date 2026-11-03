#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig10_dsr_binning.py -- σε ποιο γονιδίωμα ανήκουν τα δύο contigs του dsrAB.

    python3 fig10_dsr_binning.py [57_dsr_bin dir] [outdir]

THE QUESTION. Steps 55-59 found the lagoon's sulfate reducers on two contigs that
no binner placed. The chapter's central premise -- sulfate reducers outcompeting
methanogens for hydrogen and acetate -- rests on them, so "unbinned" is not an
acceptable place to leave them. One of the two, LIQ_contig_19299, is most
parsimoniously part of ETO_LIQ_MAG_022, and the tetranucleotide test is one of
the five lines that say so. This panel is that test.

WHY THE YARDSTICK IS PER GENOME AND NOT A SINGLE THRESHOLD. A tetranucleotide
distance means nothing in the abstract: genomes differ in how internally
consistent they are. So each MAG's own contigs, every one of them at least 5 kb,
are measured against their own genome profile, and the 90th percentile of that
spread is what a contig BELONGING to that genome typically looks like. The bar
behind each row is that genome's own p90. A dot inside its own bar is compatible
with membership; a dot outside it is not. Twenty-nine genomes, twenty-nine
different bars, because that is what the data says the yardstick is.

IT RECOMPUTES, IT DOES NOT TRANSCRIBE. The 136 canonical 4-mers, the profile, the
Euclidean distance and the >= 5 kb spread rule are step 57's code, lifted
unchanged, reading step 57's own bin_paths.tsv and dsr_contigs.fasta. Figure and
step cannot drift apart because they are the same arithmetic.

SLOW BY HONEST NECESSITY: it reads all 29 genomes to rebuild the spreads, which
takes a minute or two. Caching the distances would be faster and would reintroduce
exactly the transcription problem this avoids.
"""
import sys, os, math, itertools
import etoliko_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

args = [a for a in sys.argv[1:] if not a.startswith("--")]
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 170.0
TOPN = 10                          # nearest genomes shown per contig
MINLEN_SPREAD = 5000               # step 57's rule
MINSIZE_BIN = 50000                # step 57's rule

EL = {
    "title": "",
    "xlabel": "tetranucleotide distance to genome",
    "band": "genome's own yardstick",
    "dot": "outside: incompatible",
    "in": "inside: compatible",
    "na": "no yardstick",
    "caption": (
        "Η τετρανουκλεοτιδική συχνότητα — η συχνότητα και των %d μη πλεοναζόντων "
        "τετραμερών σε μια αλληλουχία — είναι ιδιότητα του γονιδιώματος και όχι "
        "του μεμονωμένου γονιδίου: τμήματα DNA που προέρχονται από τον ίδιο "
        "οργανισμό έχουν παρόμοιο προφίλ, και σε αυτή την αρχή στηρίζεται η "
        "ομαδοποίηση των contigs σε γονιδιώματα. Εδώ τα δύο contigs που φέρουν "
        "τα γονίδια dsrA και dsrB, και που καμία μέθοδος ομαδοποίησης δεν "
        "τοποθέτησε, συγκρίνονται με καθένα από τα %d γονιδιώματα του καταλόγου· "
        "εικονίζονται τα %d πλησιέστερα και το εύρος των υπολοίπων σημειώνεται "
        "πάνω δεξιά. Η απόσταση είναι ευκλείδεια μεταξύ των δύο προφίλ. Επειδή "
        "μια τέτοια απόσταση δεν έχει απόλυτο νόημα — τα γονιδιώματα διαφέρουν "
        "ως προς την εσωτερική τους ομοιογένεια — το μέτρο ορίζεται χωριστά για "
        "καθένα: η γκρίζα μπάρα είναι το 90ό εκατοστημόριο των αποστάσεων που "
        "έχουν τα ΔΙΚΑ ΤΟΥ contigs, όσα ξεπερνούν τα %s kb, από το δικό του "
        "προφίλ — δηλαδή το πόσο απέχει ένα τμήμα που όντως του ανήκει. Σημείο "
        "μέσα στη μπάρα είναι συμβατό με ένταξη· έξω από αυτήν, ασύμβατο. "
        "Γονιδιώματα χωρίς αρκετά δικά τους contigs αυτού του μήκους δεν "
        "αποκτούν μέτρο και σημειώνονται με τετράγωνο: για αυτά το ερώτημα δεν "
        "κρίνεται.%s %s"),
}


def pick(*cands):
    here = os.path.dirname(os.path.abspath(__file__))
    for c in cands:
        for base in ("", here, os.path.join(here, "..")):
            p = os.path.join(base, c) if base else c
            if os.path.isdir(p):
                return os.path.normpath(p)
    return None


D57 = args[0] if args and os.path.isdir(args[0]) else pick("57_dsr_bin")
if not D57:
    sys.exit("fig10: cannot find 57_dsr_bin/. Run step 57, or pass the path.")
BINS_TSV = os.path.join(D57, "bin_paths.tsv")
QFA = os.path.join(D57, "dsr_contigs.fasta")
for p in (BINS_TSV, QFA):
    if not os.path.exists(p):
        sys.exit("fig10: %s missing -- step 57 writes it.\n       %s holds: %s"
                 % (p, D57, ", ".join(sorted(os.listdir(D57))[:15])))
print("   reading %s" % D57)

# ---------------------------------------------------------------- step 57's code
def rc(s):
    c = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(c.get(b, "N") for b in reversed(s))


KM, seen = [], set()
for t in itertools.product("ACGT", repeat=4):
    k = "".join(t)
    if k in seen:
        continue
    seen.add(k); seen.add(rc(k)); KM.append(k)
IDX = {}
for i, k in enumerate(KM):
    IDX[k] = i; IDX[rc(k)] = i
print("   %d canonical 4-mers" % len(KM))


def read_fasta(p):
    name, buf = None, []
    with open(p, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def prof(seq):
    v = [0.0] * len(KM); n = 0
    seq = seq.upper()
    for i in range(len(seq) - 3):
        j = IDX.get(seq[i:i + 4])
        if j is not None:
            v[j] += 1; n += 1
    return [x / n for x in v] if n else None


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def q(a, p):
    if not a:
        return float("nan")
    i = min(len(a) - 1, max(0, int(round(p * (len(a) - 1)))))
    return a[i]


def short(t):
    parts = [x for x in t.split(";") if len(x.split("__")) > 1 and x.split("__")[1]]
    return " / ".join(p.split("__")[1] for p in parts[-2:]) if parts else t


queries = [(nm, len(sq), prof(sq)) for nm, sq in read_fasta(QFA)]
print("   %d query contig(s): %s" % (len(queries), ", ".join(n for n, _, _ in queries)))

bins, skipped = [], []
for line in open(BINS_TSV, errors="replace"):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 3:
        continue
    b, path, tax = parts[0], parts[1], parts[2]
    if not os.path.exists(path):
        print("   !! missing bin fasta, skipped: %s" % path)
        continue
    seqs = list(read_fasta(path))
    whole = "".join(s for _, s in seqs)
    if len(whole) < MINSIZE_BIN:
        skipped.append((b, len(whole)))
        continue
    wp = prof(whole)
    spread = sorted(dist(prof(s), wp) for _, s in seqs if len(s) >= MINLEN_SPREAD)
    bins.append({"name": b, "tax": short(tax), "prof": wp, "spread": spread})
print("   %d genome(s) profiled" % len(bins))
if skipped:
    print("   %d excluded, under %d bp: %s"
          % (len(skipped), MINSIZE_BIN,
             ", ".join("%s (%s bp)" % (n, format(L, ",")) for n, L in skipped)))
if not bins:
    sys.exit("fig10: no bins profiled -- check the paths in %s" % BINS_TSV)

panels = []
for nm, L, p in queries:
    rows = sorted(((dist(p, b["prof"]), b) for b in bins), key=lambda r: r[0])
    top = rows[:TOPN]
    rest = rows[TOPN:]
    compat = [(d, b) for d, b in rows
              if b["spread"] and q(b["spread"], 0.90) > 0
              and d <= q(b["spread"], 0.90)]
    nat = [b["name"] for _, b in rows if not b["spread"] or q(b["spread"], 0.90) <= 0]
    print()
    print("   %s  (%s bp)" % (nm, format(L, ",")))
    for i, (d, b) in enumerate(top[:5], 1):
        p90 = q(b["spread"], 0.90) if b["spread"] else float("nan")
        print("     %d. %-26s %.4f   own p90 %.4f  %s"
              % (i, b["name"], d, p90, "COMPATIBLE" if d <= p90 else ""))
    if len(top) > 1 and top[0][0] > 0:
        print("     runner-up is %.2fx further than the nearest"
              % (top[1][0] / top[0][0]))
    print("     %d of %d genome(s) compatible with membership: %s"
          % (len(compat), len(rows), ", ".join(b["name"] for _, b in compat) or "none"))
    if nat:
        print("     no yardstick (too few own contigs >= %d bp): %s"
              % (MINLEN_SPREAD, ", ".join(nat)))
    panels.append({"q": nm, "L": L, "top": top, "compat": compat,
                   "rest": (rest[0][0], rest[-1][0]) if rest else None,
                   "ratio": (top[1][0] / top[0][0]) if len(top) > 1 and top[0][0] else None})

# ---------------------------------------------------------------- layout
S.apply()
W = FIG_W_MM / 25.4
ROW = 0.158
GAP = 0.46
TOP_IN, BOT_IN = 0.46, 2.36
H = TOP_IN + len(panels) * (TOPN * ROW + GAP) + BOT_IN
fig = plt.figure(figsize=(W, H))
LEFT = 2.62

XMAX = max(max(d for d, _ in p["top"]) for p in panels)
XMAX = max(XMAX, max((q(b["spread"], 0.90) for p in panels for _, b in p["top"]
                      if b["spread"]), default=0)) * 1.16

for k, pn in enumerate(panels):
    bot = BOT_IN + (len(panels) - 1 - k) * (TOPN * ROW + GAP)
    ax = fig.add_axes([LEFT / W, bot / H, (W - LEFT - 0.34) / W, TOPN * ROW / H])
    S.recessive_grid(ax, axis="x")
    S.despine(ax, keep=("bottom",))
    for i, (d, b) in enumerate(pn["top"]):
        y = TOPN - 1 - i
        p90 = q(b["spread"], 0.90) if b["spread"] else 0.0
        testable = p90 > 0
        inside = testable and d <= p90
        ax.barh(y, p90, height=0.60, color=S.GRID, linewidth=0, zorder=2)
        # FILLED OR OPEN, NOT A COLOUR. Blue and orange mean Bacteria and Archaea
        # in every other panel of this set; borrowing them here for "compatible"
        # would make one reader in ten read this figure as a domain figure.
        ax.scatter(d, y, s=30, zorder=5,
                   marker="o" if testable else "s",
                   facecolors=S.INK if inside else "none",
                   edgecolors=S.INK if inside else
                              (S.INK_MUTED if testable else S.RULE),
                   linewidths=0 if inside else 1.0)
        # id and taxon both in the label column, clear of the data
        ax.text(-(LEFT / (W - LEFT - 0.34)) + 0.012, y, b["name"].replace("ETO_", ""),
                ha="left", va="center", fontsize=S.PT_TICK,
                color=S.INK if inside else S.INK_2,
                transform=ax.get_yaxis_transform(), clip_on=False)
        ax.text(-0.016, y, b["tax"][:30], ha="right", va="center", fontsize=6.4,
                style="italic", color=S.INK_MUTED,
                transform=ax.get_yaxis_transform(), clip_on=False)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(-0.8, TOPN - 0.2)
    ax.set_yticks([])
    ax.set_title("%s   (%s bp)" % (pn["q"], format(pn["L"], ",").replace(",", ".")),
                 loc="left", fontsize=S.PT_LABEL, pad=5, color=S.INK)
    if k == len(panels) - 1:
        ax.set_xlabel(EL["xlabel"])
    else:
        ax.set_xticklabels([])
    if pn["rest"]:
        # on the title line, where nothing else is competing for the space
        ax.text(1.0, 1.02, "remaining %d: %s – %s"
                % (len(bins) - TOPN, S.num(pn["rest"][0], 4), S.num(pn["rest"][1], 4)),
                ha="right", va="bottom", fontsize=6.4, color=S.INK_MUTED,
                transform=ax.transAxes)

fig.text(0.012, 1 - 0.20 / H, EL["title"], ha="left", va="top",
         fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)

fig.legend(handles=[
    Patch(facecolor=S.GRID, label=EL["band"]),
    Line2D([], [], marker="o", linestyle="none", markersize=5.2,
           markerfacecolor="none", markeredgecolor=S.INK_MUTED,
           markeredgewidth=1.0, label=EL["dot"]),
    Line2D([], [], marker="o", linestyle="none", markersize=5.2,
           markerfacecolor=S.INK, markeredgewidth=0, label=EL["in"]),
    Line2D([], [], marker="s", linestyle="none", markersize=5.0,
           markerfacecolor="none", markeredgecolor=S.RULE,
           markeredgewidth=1.0, label=EL["na"])],
    loc="lower left", bbox_to_anchor=(0.012, (BOT_IN - 0.60) / H),
    ncols=4, handlelength=1.0, handleheight=1.0, borderpad=0.25,
    columnspacing=1.2, handletextpad=0.5, frameon=False)


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


verdict = []
for pn in panels:
    if len(pn["compat"]) == 1:
        nm_ = pn["compat"][0][1]["name"].replace("ETO_", "")
        r_ = (" και το αμέσως επόμενο γονιδίωμα είναι %s φορές πιο απομακρυσμένο"
              % S.num(pn["ratio"], 2)) if pn["ratio"] else ""
        verdict.append("Το %s είναι συμβατό μόνο με το %s%s."
                       % (pn["q"], nm_, r_))
    elif not pn["compat"]:
        verdict.append("Το %s δεν είναι συμβατό με κανένα γονιδίωμα του "
                       "καταλόγου." % pn["q"])
    else:
        verdict.append("Το %s είναι συμβατό με %d γονιδιώματα, οπότε η δοκιμή "
                       "αυτή δεν το τοποθετεί μόνη της."
                       % (pn["q"], len(pn["compat"])))

fig.canvas.draw()
SKIPTXT = ("" if not skipped else
           " Εξαιρέθηκαν %d γονιδιώματα μικρότερα από %s kb, για τα οποία το "
           "προφίλ δεν είναι αξιόπιστο." % (len(skipped),
                                            S.num(MINSIZE_BIN / 1000.0, 0)))
caption_fit(fig, EL["caption"] % (len(KM), len(bins), TOPN,
                                  S.num(MINLEN_SPREAD / 1000.0, 0),
                                  SKIPTXT, " ".join(verdict)))

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print()
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))

S.save(fig, "S_dsr_binning", OUT)
