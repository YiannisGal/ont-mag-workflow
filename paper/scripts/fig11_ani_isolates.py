#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig11_ani_isolates.py -- the catalogue's genomes are not the cultured ones.

    python3 fig11_ani_isolates.py [40_integration dir] [outdir]

Reads the run directly:
    skani_raw.tsv   the pair table (ANI + both aligned fractions)
    mag_list.txt    the genomes submitted as MAGs
    iso_list.txt    the cultured genomes submitted alongside them

THREE GROUPS, NOT TWO. iso_list.txt holds two different kinds of file: genomes
of strains isolated in this work, and PUBLIC REFERENCE GENOMES of named species
used to identify them. Counting them together turns two different claims --
"no MAG is one of our isolates" and "no MAG is one of these named species" --
into one meaningless number, so they are separated and every one of the six
resulting blocks is reported. The split is by filename (a Genus_species
binomial is a reference) and is PRINTED for checking; a file named refs.txt in
the same directory overrides it, one reference filename per line.

WHY THIS IS NOT A SCATTER OF ANI VALUES. One point per MAG at its best ANI
against any cultured genome cannot be drawn, because no such value exists for
any MAG. The method reports a pair only when it clears an internal screening
threshold, so those pairs were computed and returned nothing. A scatter would
have to invent a number -- zero, or the axis minimum -- and the invented number
would then be read as a measurement. The figure shows the absence itself.

IDENTITY IS NOT TRANSITIVE, AND THIS DATA PROVES IT. An earlier version merged
files into "one organism" by connected components at >=99 % ANI. That is wrong
whenever one file is a MIXTURE: a mixed assembly scores ~100 % against each of
its own fractions while the fractions score nothing against each other, so
A~B and A~C without B~C. Components are now tested for being cliques, and a
file adjacent to every other member of a non-clique component is reported as a
mixture and excluded from the count of independent organisms rather than being
allowed to glue them together. The aligned fractions are printed for each such
pair, because a subset aligns almost entirely to its parent while the parent
aligns only partly to the subset -- that asymmetry is what distinguishes
"contained in" from "identical to".

THE FLOOR COMES FROM THE DATA. The lowest ANI anywhere in the table is the
lowest value this run demonstrably reports; every pair with no value lies below
it. That bound needs no knowledge of the tool's defaults, so it cannot be wrong
about them.
"""
import sys, os, re, math, itertools
import etoliko_figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SRC = args[0] if args else "40_integration"
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 160.0

SPECIES = 95.0      # Jain et al. 2018 -- the species boundary
DUP     = 99.0      # above this, two cultured FILES are the same genome

EL = {
    "title": "Ταυτότητα γονιδιωμάτων: MAG, απομονωμένα στελέχη, γονιδιώματα αναφοράς",
    "blk": {"MM": "MAG × MAG",
            "IM": "MAG × στέλεχος",
            "MR": "MAG × αναφορά",
            "II": "στέλεχος × στέλεχος",
            "IR": "στέλεχος × αναφορά",
            "RR": "αναφορά × αναφορά"},
    "b_cmp": "συγκρίθηκαν",
    "b_val": "με τιμή",
    "b_max": "μέγ. ANI",
    "b_sp":  "≥ 95 %",
    "none":  "καμία",
    "xlabel": "ANI (%)",
    "sp":    "όριο είδους, 95 %",
    "ctrl":  "θετικός μάρτυρας",
    "k_mm":  "MAG – MAG",
    "k_ii":  "στέλεχος – στέλεχος",
    "k_ir":  "στέλεχος – αναφορά (ταυτοποίηση)",
    "k_rr":  "αναφορά – αναφορά",
    "k_mi":  "MAG – καλλιεργημένου: καμία τιμή",
}

# ---------------------------------------------------------------- inputs
if not os.path.isdir(SRC):
    here = os.path.dirname(os.path.abspath(__file__))
    for c in (os.path.join(here, SRC), os.path.join(here, "..", SRC),
              os.path.join(here, "..", "40_integration")):
        if os.path.isdir(c):
            SRC = os.path.normpath(c); print("   found %s" % SRC); break
    else:
        sys.exit("fig11: cannot find the 40_integration directory -- pass it "
                 "as the first argument")

P_TAB = os.path.join(SRC, "skani_raw.tsv")
P_MAG = os.path.join(SRC, "mag_list.txt")
P_ISO = os.path.join(SRC, "iso_list.txt")
P_REF = os.path.join(SRC, "refs.txt")
for p in (P_TAB, P_MAG, P_ISO):
    if not os.path.exists(p):
        sys.exit("fig11: missing %s -- this figure is drawn from the run, not "
                 "from transcribed numbers, so it refuses to guess." % p)


def names(path):
    return [os.path.basename(l.strip())
            for l in open(path, encoding="utf-8", errors="replace") if l.strip()]


MAGS = names(P_MAG)
CULT = names(P_ISO)

# --- reference or isolate?
BINOMIAL = re.compile(r"^[A-Z][a-z]+_[a-z]{3,}$")
if os.path.exists(P_REF):
    REFS = [g for g in CULT if g in set(names(P_REF))]
    print("   reference set taken from refs.txt")
else:
    REFS = [g for g in CULT if BINOMIAL.match(g.rsplit(".", 1)[0])]
ISOS = [g for g in CULT if g not in set(REFS)]
MSET, ISET, RSET = set(MAGS), set(ISOS), set(REFS)

print("   submitted: %d MAG files, %d cultured files" % (len(MAGS), len(CULT)))
print("   of the cultured files, %d read as isolates and %d as named references"
      % (len(ISOS), len(REFS)))
print("     isolates   : %s" % ", ".join(ISOS))
print("     references : %s" % (", ".join(REFS) or "(none)"))
if not REFS:
    print("   no reference genomes recognised; the figure falls back to three")
    print("   blocks. If that is wrong, list them in %s" % P_REF)

# ---------------------------------------------------------------- the table
# Both directions are emitted. ANI is symmetric to rounding; the ALIGNED
# FRACTIONS are not, and their asymmetry is what tells a subset from a twin, so
# both directions are kept for af and collapsed only for ANI.
ani, af = {}, {}
rows = 0
for ln in open(P_TAB, encoding="utf-8", errors="replace"):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 3 or f[0] == "Ref_file":
        continue
    try:
        v = float(f[2])
    except ValueError:
        continue
    rows += 1
    a, b = os.path.basename(f[0]), os.path.basename(f[1])
    if a == b:
        continue
    try:
        af[(a, b)] = (float(f[3]), float(f[4]))      # (fraction of a, of b)
    except (IndexError, ValueError):
        pass
    k = (a, b) if a <= b else (b, a)
    if k not in ani or v > ani[k]:
        ani[k] = v
print("   table: %d rows, %d distinct non-self pairs" % (rows, len(ani)))

unknown = sorted({g for k in ani for g in k} - MSET - ISET - RSET)
if unknown:
    sys.exit("fig11: %d genome(s) in the table are on no list (%s) -- refusing "
             "to draw, because a genome that belongs to no group cannot be "
             "counted into any block" % (len(unknown), ", ".join(unknown[:6])))


def grp(g):
    return "M" if g in MSET else ("I" if g in ISET else "R")


def block(k):
    return "".join(sorted(grp(k[0]) + grp(k[1])))


# ------------------------------------------------- independent cultured genomes
# Connected components at >= DUP, then each component TESTED FOR BEING A CLIQUE.
adj = {g: set() for g in ISOS}
for (a, b), v in ani.items():
    if a in ISET and b in ISET and v >= DUP:
        adj[a].add(b); adj[b].add(a)

seen, comps = set(), []
for g in ISOS:
    if g in seen:
        continue
    stack, comp = [g], []
    seen.add(g)
    while stack:
        x = stack.pop(); comp.append(x)
        for y in adj[x] - seen:
            seen.add(y); stack.append(y)
    comps.append(sorted(comp))

MIXED, ORGS = [], []
for comp in comps:
    if len(comp) == 1:
        ORGS.append(comp); continue
    clique = all(b in adj[a] for a, b in itertools.combinations(comp, 2))
    if clique:
        ORGS.append(comp); continue
    # not a clique: whatever is adjacent to every other member is a mixture
    # holding them, not a duplicate of them
    bridges = [x for x in comp if len(adj[x] & set(comp)) == len(comp) - 1]
    rest = [x for x in comp if x not in bridges]
    print("   ! %s is not one organism: %s"
          % (" + ".join(comp), "; ".join(
              "%s—%s %s" % (a, b, "yes" if b in adj[a] else "NO VALUE")
              for a, b in itertools.combinations(comp, 2))))
    for x in bridges:
        MIXED.append(x)
        for y in rest:
            k = (x, y) if x <= y else (y, x)
            pa = af.get((x, y)) or af.get((y, x))
            if pa:
                fx, fy = pa if (x, y) in af else (pa[1], pa[0])
                print("     aligned: %s %s%% of its length, %s %s%% of its "
                      "length -> %s is contained in %s"
                      % (x, S.num(fx, 1), y, S.num(fy, 1), y, x))
    # the bridge removed, each remaining file stands on its own
    for y in rest:
        ORGS.append([y])

N_ORG = len(ORGS)
print("   %d isolate file(s) -> %d independent organism(s)%s"
      % (len(ISOS), N_ORG,
         ("; %s excluded as mixed assembl%s"
          % (", ".join(MIXED), "y" if len(MIXED) == 1 else "ies")) if MIXED else ""))

# ------------------------------------------------- what each isolate matches
IDENT = []
for comp in ORGS:
    best = None
    for g in comp:
        for r in REFS:
            k = (g, r) if g <= r else (r, g)
            v = ani.get(k)
            if v is not None and v >= SPECIES and (best is None or v > best[1]):
                best = (r, v)
    IDENT.append((comp, best))
    print("     %-32s -> %s"
          % (" + ".join(comp),
             "%s  %s %%" % (best[0].rsplit(".", 1)[0], S.num(best[1], 2))
             if best else "no reference match"))
NAMED = [c for c, b in IDENT if b]
UNNAMED = [c for c, b in IDENT if not b]

# ---------------------------------------------------------------- blocks
n = {"M": len(MAGS), "I": len(ISOS), "R": len(REFS)}
CMP = {"MM": n["M"] * (n["M"] - 1) // 2, "II": n["I"] * (n["I"] - 1) // 2,
       "RR": n["R"] * (n["R"] - 1) // 2, "IM": n["M"] * n["I"],
       "MR": n["M"] * n["R"], "IR": n["I"] * n["R"]}
VAL = dict.fromkeys(CMP, 0)
MAXV = dict.fromkeys(CMP, None)
NSP = dict.fromkeys(CMP, 0)
for k, v in ani.items():
    bl = block(k)
    VAL[bl] += 1
    if MAXV[bl] is None or v > MAXV[bl]:
        MAXV[bl] = v
    if v >= SPECIES:
        NSP[bl] += 1

CLAIM = ("IM", "MR")            # the blocks the chapter's sentence rests on
ORDER = ["MM", "IM", "MR", "II", "IR", "RR"] if REFS else ["MM", "IM", "II"]
FLOOR = min(ani.values()) if ani else None
print("   blocks (compared / with a value / max / >=95%):")
for bl in ORDER:
    print("     %-22s %5d / %3d / %s / %d"
          % (EL["blk"][bl], CMP[bl], VAL[bl],
             "--" if MAXV[bl] is None else "%.2f" % MAXV[bl], NSP[bl]))
print("   lowest value anywhere in the table: %s"
      % ("--" if FLOOR is None else "%.2f %%" % FLOOR))
CLAIM_CMP = sum(CMP[b] for b in CLAIM if b in ORDER)
CLAIM_VAL = sum(VAL[b] for b in CLAIM if b in ORDER)

# ---------------------------------------------------------------- caption
cap = ["Κάθε γονιδίωμα του καταλόγου συγκρίθηκε με κάθε άλλο, με κάθε στέλεχος "
       "που απομονώθηκε σε αυτή την εργασία και με κάθε δημόσιο γονιδίωμα "
       "αναφοράς, σε έναν ενιαίο υπολογισμό όλων προς όλους. Το ANI (μέση "
       "νουκλεοτιδική ταυτότητα) μετρά πόσο όμοια είναι δύο γονιδιώματα· τιμή "
       "≥ 95 % σημαίνει το ίδιο είδος. Στελέχη και γονιδιώματα αναφοράς "
       "μετρώνται χωριστά, γιατί απαντούν σε δύο διαφορετικά ερωτήματα.",
       "Τα πλαίσια δίνουν χωριστά ΠΟΣΑ ΖΕΥΓΗ ΣΥΓΚΡΙΘΗΚΑΝ και ΠΟΣΑ ΕΠΕΣΤΡΕΨΑΝ "
       "ΤΙΜΗ: η μέθοδος αναφέρει ένα ζεύγος μόνο όταν ξεπερνά ένα εσωτερικό "
       "κατώφλι ομοιότητας, οπότε ένα ζεύγος χωρίς τιμή έχει υπολογιστεί και "
       "βρεθεί κάτω από αυτό — δεν έχει παραλειφθεί."]

if CLAIM_VAL == 0:
    cap.append("Και τα %s ζεύγη MAG – καλλιεργημένου γονιδιώματος υπολογίστηκαν, "
               "κανένα δεν επέστρεψε τιμή. Η χαμηλότερη τιμή που αποδίδει η "
               "μέθοδος οπουδήποτε εδώ είναι %s %%, άρα καθένα τους βρίσκεται "
               "κάτω από εκείνο το όριο, πολύ χαμηλότερα από το όριο του "
               "είδους: κανένα από τα %d γονιδιώματα του καταλόγου δεν "
               "αντιστοιχεί ούτε σε στέλεχος που απομονώθηκε ούτε σε είδος "
               "αναφοράς." % (S.num(CLAIM_CMP, 0), S.num(FLOOR, 2), len(MAGS)))
else:
    cap.append("Στα ζεύγη MAG – καλλιεργημένου γονιδιώματος επέστρεψαν τιμή %s "
               "από τα %s." % (S.num(CLAIM_VAL, 0), S.num(CLAIM_CMP, 0)))

cap.append("Το κάτω μέρος δείχνει κάθε ζεύγος με τιμή και λειτουργεί ως μάρτυρας "
           "του πάνω: η μέθοδος αποδίδει τιμές σε όλο το εύρος όταν δύο "
           "γονιδιώματα συγγενεύουν, οπότε τα κενά μπλοκ είναι ιδιότητα των "
           "γονιδιωμάτων και όχι αδυναμία της μέτρησης.")

if NAMED:
    cap.append("Ισχυρότερος μάρτυρας είναι οι ταυτοποιήσεις των ίδιων των "
               "στελεχών: %s. Η μέθοδος βρίσκει αντιστοιχίσεις σε επίπεδο "
               "είδους όταν υπάρχουν, και παρ' όλα αυτά δεν βρίσκει καμία για "
               "κανένα MAG."
               % "· ".join("%s = %s (%s %%)"
                           % (" + ".join(x.rsplit(".", 1)[0] for x in c),
                              b[0].rsplit(".", 1)[0].replace("_", " "),
                              S.num(b[1], 2))
                           for c, b in IDENT if b))
if UNNAMED:
    cap.append("Το %s δεν επέστρεψε τιμή έναντι κανενός γονιδιώματος, ούτε "
               "αναφοράς ούτε MAG, και παραμένει χωρίς αντιστοίχιση σε γνωστό "
               "είδος μέσα σε αυτόν τον υπολογισμό."
               % ", ".join(" + ".join(x.rsplit(".", 1)[0] for x in c)
                           for c in UNNAMED))
if MIXED:
    cap.append("ΤΑ ΑΡΧΕΙΑ ΔΕΝ ΕΙΝΑΙ ΟΡΓΑΝΙΣΜΟΙ: τα %d αρχεία στελεχών "
               "αντιστοιχούν σε %d ανεξάρτητους οργανισμούς. Το %s είναι μεικτή "
               "καλλιέργεια: εμφανίζεται με ~100 %% ταυτότητα προς καθένα από τα "
               "κλάσματα στα οποία διαχωρίστηκε, ενώ τα κλάσματα δεν "
               "επιστρέφουν τιμή μεταξύ τους — είναι διαφορετικοί οργανισμοί. Η "
               "ταυτότητα δεν είναι μεταβατική όταν το ένα αρχείο περιέχει το "
               "άλλο, γι' αυτό το μεικτό αρχείο δεν προσμετράται."
               % (len(ISOS), N_ORG,
                  ", ".join(x.rsplit(".", 1)[0] for x in MIXED)))

cap.append("ANI κάτω από το όριο του είδους σημαίνει ότι οι οργανισμοί δεν είναι "
           "οι ίδιοι· δεν σημαίνει ότι το στέλεχος απουσιάζει από τη "
           "λιμνοθάλασσα, παρά μόνο ότι δεν ανακτήθηκε ως γονιδίωμα από αυτά τα "
           "δεδομένα.")
CAPTION = " ".join(cap)

# ---------------------------------------------------------------- draw
S.apply()
W = FIG_W_MM / 25.4
ranked = sorted(ani.items(), key=lambda kv: kv[1], reverse=True)
PAIRLAB = ["%s  ×  %s" % (a.rsplit(".", 1)[0], b.rsplit(".", 1)[0])
           for (a, b), v in ranked]

# --- pass 1: measure. Nothing is placed until every band's height is known;
# sizing the figure first and adding the caption afterwards printed the caption
# straight through the panel.
probe = plt.figure(figsize=(W, 4.0))
probe.canvas.draw()


def text_w_in(f, s, pt):
    t = f.text(0.5, 0.5, s, fontsize=pt)
    w = t.get_window_extent(f.canvas.get_renderer()).width / f.dpi
    t.remove()
    return w


def caption_width(f, txt, lo_w=70, hi_w=150):
    lim = f.get_figwidth() * (1.0 - 0.012)
    best = lo_w
    for w in range(lo_w, hi_w + 1, 2):
        S.caption(f, txt, width=w)
        art = f.texts[-1]
        fits = art.get_window_extent(f.canvas.get_renderer()).x1 / f.dpi <= lim
        art.remove()
        if not fits:
            break
        best = w
    return best


LAB_W = max(text_w_in(probe, s, S.PT_VALUE) for s in PAIRLAB)
VAL_W = max(text_w_in(probe, S.num(v, 2), S.PT_VALUE) for v in ani.values())
CAP_W = caption_width(probe, CAPTION)
S.caption(probe, CAPTION, width=CAP_W)
CAP_H = probe.texts[-1].get_window_extent(probe.canvas.get_renderer()).height \
    / probe.dpi

# the key gets its own band: put inside the panel it lands on the lowest rows,
# which still carry a dot and its printed value
KEYS = [("MM", EL["k_mm"]), ("II", EL["k_ii"])]
if REFS:
    KEYS += [("IR", EL["k_ir"]), ("RR", EL["k_rr"])]
KEYS += [("CTRL", EL["ctrl"])]
KEY_W = sum(text_w_in(probe, lb, S.PT_TICK) + 0.34 for _, lb in KEYS)
LEG_NCOL = 3 if KEY_W > W * 0.98 else len(KEYS)
LEG_ROWS = int(math.ceil(len(KEYS) / float(LEG_NCOL)))
LEG_H = 0.20 * LEG_ROWS + 0.10
print("   key: %.0f mm of text -> %d column(s), %d row(s)"
      % (KEY_W * 25.4, LEG_NCOL, LEG_ROWS))

# the widest thing in a block box decides whether six boxes fit across
BOXW_IN = W / 3.0
need = max([text_w_in(probe, EL["b_cmp"], S.PT_TICK)
            + text_w_in(probe, S.num(max(CMP.values()), 0), S.PT_TICK)]
           + [text_w_in(probe, EL["blk"][b], S.PT_LABEL) for b in ORDER]) + 0.22
NCOL_A = 3 if need <= BOXW_IN else 2
NROW_A = int(math.ceil(len(ORDER) / float(NCOL_A)))
print("   blocks: widest content %.0f mm in a %.0f mm cell -> %d × %d grid"
      % (need * 25.4, BOXW_IN * 25.4, NROW_A, NCOL_A))
plt.close(probe)

GAP = 0.055
L_IN = 0.035 * W + LAB_W + GAP
R_IN = VAL_W + 0.10

ROW_IN = 0.118
M_TOP, M_BOT = 0.10, 0.10
TITLE_H, GAP1, GAP2, GAP_L, XLAB_H, GAP3 = 0.26, 0.14, 0.26, 0.12, 0.38, 0.22
A_ROW = 0.80
A_H = NROW_A * A_ROW
B_H = len(ranked) * ROW_IN + 0.34
H = (M_TOP + TITLE_H + GAP1 + A_H + GAP2 + LEG_H + GAP_L + B_H + XLAB_H
     + GAP3 + CAP_H + M_BOT)
print("   bands (in): title %.2f  blocks %.2f  key %.2f  pairs %.2f  caption %.2f"
      "  -> %.2f" % (TITLE_H, A_H, LEG_H, B_H, CAP_H, H))

fig = plt.figure(figsize=(W, H))
y_cap = M_BOT / H
y_B = (M_BOT + CAP_H + GAP3 + XLAB_H) / H
y_L = (M_BOT + CAP_H + GAP3 + XLAB_H + B_H + GAP_L) / H
y_A = (M_BOT + CAP_H + GAP3 + XLAB_H + B_H + GAP_L + LEG_H + GAP2) / H
y_T = (H - M_TOP - TITLE_H * 0.5) / H

fig.text(0.035, y_T, EL["title"], ha="left", va="center",
         fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)

# ---- (a) every block of the comparison, none aggregated away
axA = fig.add_axes([0.035, y_A, 0.945, A_H / H])
axA.set_xlim(0, NCOL_A); axA.set_ylim(0, NROW_A); axA.axis("off")
for i, bl in enumerate(ORDER):
    col, row = i % NCOL_A, i // NCOL_A
    y0 = (NROW_A - 1 - row)
    claim_empty = (bl in CLAIM and VAL[bl] == 0)
    axA.add_patch(Rectangle((col + 0.03, y0 + 0.05), 0.94, 0.90,
                            facecolor=S.SURFACE,
                            edgecolor=S.HILITE if claim_empty else S.RULE,
                            lw=1.5 if claim_empty else 0.8, zorder=2))
    axA.text(col + 0.5, y0 + 0.80, EL["blk"][bl], ha="center", va="center",
             fontsize=S.PT_LABEL, color=S.HILITE if claim_empty else S.INK,
             zorder=3)
    lines = [(EL["b_cmp"], S.num(CMP[bl], 0)),
             (EL["b_val"], EL["none"] if VAL[bl] == 0 else S.num(VAL[bl], 0)),
             (EL["b_max"], "—" if MAXV[bl] is None else "%s %%" % S.num(MAXV[bl], 2)),
             (EL["b_sp"], EL["none"] if NSP[bl] == 0 else S.num(NSP[bl], 0))]
    for j, (a, b) in enumerate(lines):
        yy = y0 + 0.60 - j * 0.145
        axA.text(col + 0.10, yy, a, ha="left", va="center", fontsize=S.PT_TICK,
                 color=S.INK_MUTED, zorder=3)
        axA.text(col + 0.90, yy, b, ha="right", va="center", fontsize=S.PT_TICK,
                 color=S.HILITE if (claim_empty and j in (1, 3)) else S.INK,
                 zorder=3)

# ---- (b) every pair that returned a value
axB = fig.add_axes([L_IN / W, y_B, (W - L_IN - R_IN) / W, B_H / H])
S.despine(axB, keep=("bottom",))
S.recessive_grid(axB, axis="x")
lo = math.floor((FLOOR - 0.6) * 2) / 2 if FLOOR is not None else 80.0
axB.set_xlim(lo, 100.0)
axB.set_ylim(len(ranked) - 0.5, -1.35)

# colour says WHICH SIDE (catalogue or cultured set), shape says WHAT KIND of
# comparison. Spending a new hue on the reference pairs would break the domain
# colours the rest of the chapter uses.
STYLE = {"MM": dict(marker="o", color=S.DOMAIN["Bacteria"], fill=True),
         "II": dict(marker="o", color=S.INK_2, fill=True),
         "IR": dict(marker="s", color=S.INK_2, fill=True),
         "RR": dict(marker="s", color=S.INK_2, fill=False)}
CTRL = None
for y, (((a, b), v), lab) in enumerate(zip(ranked, PAIRLAB)):
    st = STYLE.get(block((a, b)), STYLE["MM"])
    axB.plot([lo, v], [y, y], color=S.GRID, lw=0.7, zorder=1)
    if CTRL is None and block((a, b)) == "MM" and v >= SPECIES:
        CTRL = (y, v)
        axB.scatter([v], [y], s=116, facecolors="none", edgecolors=S.INK,
                    linewidths=1.0, zorder=3, clip_on=False)
    # clip_on=False: a pair at exactly 100 sits on the axis edge and would
    # otherwise be drawn as half a marker
    axB.scatter([v], [y], s=34, marker=st["marker"], zorder=4, clip_on=False,
                facecolors=st["color"] if st["fill"] else "none",
                edgecolors=S.SURFACE if st["fill"] else st["color"],
                linewidths=0.6 if st["fill"] else 1.0)
    axB.annotate(S.num(v, 2), (v, y), textcoords="offset points", xytext=(7, 0),
                 ha="left", va="center", fontsize=S.PT_VALUE, color=S.INK_2,
                 zorder=4, annotation_clip=False)
    axB.annotate(lab, (lo, y), textcoords="offset points", xytext=(-4, 0),
                 ha="right", va="center", fontsize=S.PT_VALUE, color=S.INK,
                 zorder=4, annotation_clip=False)

axB.axvline(SPECIES, color=S.INK, lw=0.9, ls=(0, (4, 2)), zorder=2)
axB.text(SPECIES, -1.28, EL["sp"], ha="center", va="bottom",
         fontsize=S.PT_VALUE, color=S.INK)
if FLOOR is not None:
    axB.axvline(FLOOR, color=S.RULE, lw=0.9, zorder=2)
axB.set_xlabel(EL["xlabel"])
axB.set_yticks([])

h = []
for tag, lb in KEYS:
    if tag == "CTRL":
        if CTRL is None:
            continue
        h.append(Line2D([], [], marker="o", linestyle="none", markersize=7.0,
                        markerfacecolor="none", markeredgecolor=S.INK,
                        markeredgewidth=1.0, label=lb))
    else:
        st = STYLE[tag]
        h.append(Line2D([], [], marker=st["marker"], linestyle="none",
                        markersize=5.2,
                        markerfacecolor=st["color"] if st["fill"] else "none",
                        markeredgecolor=S.SURFACE if st["fill"] else st["color"],
                        markeredgewidth=0.6 if st["fill"] else 1.0, label=lb))
axL = fig.add_axes([0.035, y_L, 0.945, LEG_H / H])
axL.axis("off")
axL.legend(handles=h, loc="center", ncol=LEG_NCOL, frameon=False,
           handlelength=1.0, borderpad=0.0, labelspacing=0.4,
           columnspacing=1.4, fontsize=S.PT_TICK)

S.caption(fig, CAPTION, y=y_cap, width=CAP_W)

fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered width: %.1f mm (target %.0f mm), height %.1f mm"
      % (bb.width * 25.4, FIG_W_MM, bb.height * 25.4))
if bb.width * 25.4 > FIG_W_MM + 0.5:
    print("   ! wider than the text block -- shorten the pair labels")

S.save(fig, "fig11_ani_isolates", OUT)
