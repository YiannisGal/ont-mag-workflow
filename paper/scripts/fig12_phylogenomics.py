#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig12_phylogenomics.py -- where the 29 genomes sit among GTDB references.

    python3 scripts/fig12_phylogenomics.py . figures

Writes TWO figures, archaea and bacteria, from one load of the data:
    fig12a_phylogenomics_archaea    ar53_ref   34 tips,  8 ours + 26 GTDB
    fig12b_phylogenomics_bacteria   bac120_ref 97 tips, 21 ours + 76 GTDB

IT USES ref_trees/, NOT trees/. 12_phylogeny/trees/ar53 contains our eight
archaea and NO reference genomes at all, so it shows relationships among them
and not their placement; the same tree with references is ref_trees/ar53_ref.
A MAG chapter's phylogenomic figure has to show placement.

Reference tips are named from *_ref_labels.tsv where that file provides a
label, and from the accession otherwise -- with the caption saying how many
fell back, rather than printing bare accessions silently.

Our genomes are identified by joining the tip names against MASTER.tsv, the
same authority fig14/fig15 use, because a tip label that names a genome has to
be one the curated table confirms.
"""
import sys, os, re, math, csv, collections
import etoliko_figstyle as S
import etoliko_tree as T
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS2 = args[0] if args else "."
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 160.0
SUPPORT_MIN_FRAC = 0.70
MIN_COLLAPSE = 3
ROW = 0.088

EL = {
    "ta": "Φυλογενωμική θέση των αρχαίων γονιδιωμάτων του καταλόγου",
    "tb": "Φυλογενωμική θέση των βακτηριακών γονιδιωμάτων του καταλόγου",
    "ours": "γονιδιώματα αυτής της εργασίας",
    "ref": "γονιδιώματα αναφοράς GTDB",
    "coll": "συμπτυγμένη ομάδα αναφοράς",
    "scale": "υποκ./θέση",
    "out": "εξωομάδα",
    "nogen": "(χωρίς γένος στη GTDB)",
}

P = os.path.join(WS2, "12_phylogeny", "ref_trees")
if not os.path.isdir(P):
    sys.exit("fig12: %s not found -- pass the workspace_2 root as the first "
             "argument" % P)

# ------------------------------------------------------------------- MASTER
MAG, GENUS, FAM, DOM, COMP, CONT, MIMAG = {}, {}, {}, {}, {}, {}, {}
mp = os.path.join(WS2, "47_master", "MASTER.tsv")
if not os.path.exists(mp):
    sys.exit("fig12: %s not found; it is what confirms a tip names a genome" % mp)
for r in csv.DictReader(open(mp, errors="replace"), delimiter="\t"):
    wn = (r.get("working_name") or "").strip()
    mg = (r.get("mag_id") or "").strip()
    if not mg:
        continue
    for k in (wn, mg, mg.replace("_", ""), wn.replace("_", "")):
        if k:
            MAG[k] = mg
    GENUS[mg] = (r.get("genus") or "").strip()
    FAM[mg] = (r.get("family") or "").strip()
    DOM[mg] = (r.get("domain") or "?").strip()
    for src, dst in (("completeness", COMP), ("contamination", CONT)):
        try:
            dst[mg] = float(r.get(src) or "nan")
        except ValueError:
            dst[mg] = float("nan")
    MIMAG[mg] = (r.get("mimag") or "").strip()
print("   MASTER.tsv: %d genome(s), %d archaea, %d bacteria"
      % (len(GENUS), sum(1 for v in DOM.values() if v == "Archaea"),
         sum(1 for v in DOM.values() if v == "Bacteria")))


def ours(tip):
    """the mag_id this tip names, or None. Tip names vary between runs, so
    every key MASTER knows is tried against the whole name and its fields."""
    cand = [tip] + re.split(r"[|~\s]+", tip)
    for c in cand:
        if c in MAG:
            return MAG[c]
    for k, v in MAG.items():
        if len(k) > 6 and k in tip:
            return v
    return None


def load(stem):
    d = {}
    tp = None
    for c in ("%s.contree" % stem, "%s.treefile" % stem):
        if os.path.exists(os.path.join(P, c)):
            tp = os.path.join(P, c); break
    if tp is None:
        return None
    root = T.parse(open(tp, errors="replace").read())
    d["tree"], d["path"] = root, tp
    # labels
    lab = {}
    lp = os.path.join(P, "%s_labels.tsv" % stem)
    if os.path.exists(lp):
        for ln in open(lp, errors="replace"):
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and f[0].lower() not in ("tip", "label", "id"):
                lab[f[0].strip()] = f[1].strip()
    d["labels"] = lab
    # outgroup, if the run recorded one
    og = []
    op = os.path.join(P, "%s_outgroup.txt" % stem.replace("_ref", ""))
    if not os.path.exists(op):
        op = os.path.join(P, "%s_outgroup.txt" % stem)
    if os.path.exists(op):
        og = [x.strip() for x in open(op, errors="replace").read().split() if x.strip()]
    d["outgroup"] = og
    # model / sites, for the caption
    info = {}
    ip = os.path.join(P, "%s.iqtree" % stem)
    if os.path.exists(ip):
        txt = open(ip, errors="replace").read()
        m = re.search(r"Input data:\s*(\d+)\s+sequences with\s+(\d+)\s+amino-acid sites", txt)
        if m:
            info["n"], info["sites"] = int(m.group(1)), int(m.group(2))
        m = re.search(r"Model of substitution:\s*(\S+)", txt)
        if m:
            info["model"] = m.group(1)
        m = re.search(r"ultrafast bootstrap \((\d+) replicates\)", txt)
        if m:
            info["boot"] = int(m.group(1))
    d["info"] = info
    return d


PANELS = []
for stem, key, title in (("ar53_ref", "archaea", EL["ta"]),
                         ("bac120_ref", "bacteria", EL["tb"])):
    d = load(stem)
    if d is None:
        print("   ! no tree for %s -- skipping that panel" % stem)
        continue
    tips = d["tree"].tips()
    mine = [t for t in tips if ours(t.name)]
    print("   %-12s %s: %d tips, %d of ours, %d reference(s)  %s"
          % (stem, os.path.basename(d["path"]), len(tips), len(mine),
             len(tips) - len(mine),
             "model %s, %s sites" % (d["info"].get("model", "?"),
                                     d["info"].get("sites", "?"))))
    if not mine:
        print("     ! no tip of this tree matches a genome in MASTER.tsv.")
        print("       Refusing to draw a placement figure in which our own")
        print("       genomes cannot be identified. First few tip names:")
        for t in tips[:6]:
            print("         %s" % t.name)
        continue
    nlab = sum(1 for t in tips if t.name in d["labels"])
    print("     %d of %d tip(s) have a label in %s_labels.tsv"
          % (nlab, len(tips), stem))
    d["stem"], d["key"], d["title"], d["mine"] = stem, key, title, mine
    PANELS.append(d)

if not PANELS:
    sys.exit("fig12: nothing drawable")

S.apply()
W = FIG_W_MM / 25.4
probe = plt.figure(figsize=(W, 4.0)); probe.canvas.draw()


def tw(s, pt):
    t = probe.text(0.5, 0.5, s, fontsize=pt)
    w = t.get_window_extent(probe.canvas.get_renderer()).width / probe.dpi
    t.remove()
    return w


def caption_width(f, txt, lo=70, hi=150):
    lim = f.get_figwidth() * (1.0 - 0.012)
    best = lo
    for w in range(lo, hi + 1, 2):
        S.caption(f, txt, width=w)
        a = f.texts[-1]
        ok = a.get_window_extent(f.canvas.get_renderer()).x1 / f.dpi <= lim
        a.remove()
        if not ok:
            break
        best = w
    return best


def cap_geom(txt):
    w = caption_width(probe, txt)
    S.caption(probe, txt, width=w)
    h = probe.texts[-1].get_window_extent(probe.canvas.get_renderer()).height
    probe.texts[-1].remove()
    return w, h / probe.dpi


def clean(s):
    """GTDB labels arrive as d__;p__;...;s__Name -- keep the informative tail"""
    s = s.strip().strip("'\"")
    if ";" in s:
        parts = [x.strip() for x in s.split(";") if x.strip()]
        for p in reversed(parts):
            v = re.sub(r"^[a-z]__", "", p)
            if v:
                s = v
                break
    return re.sub(r"\s+", " ", s.replace("_", " ")).strip()


M_TOP, M_BOT = 0.10, 0.10
TITLE_H, GAP1, GAP2, GAP3, LEG_H = 0.26, 0.16, 0.30, 0.20, 0.32

for d in PANELS:
    root, lab, key = d["tree"], d["labels"], d["key"]
    accent = S.DOMAIN["Archaea"] if key == "archaea" else S.DOMAIN["Bacteria"]

    # ---- root on the recorded outgroup where there is one. Rule 132: a tree
    # is not rooted until it is rooted, and the parser's arbitrary top is not
    # a root. If the outgroup is not a clade as parsed, say so and leave the
    # tree as IQ-TREE wrote it rather than inventing a rooting.
    # The OUTGROUP| marker is on the tip's LABEL in *_ref_labels.tsv, not on
    # the tip name in the Newick, which is an accession. Testing only the name
    # found nothing and the caption then said there was no outgroup, directly
    # contradicting the tips the figure was drawing.
    def _og(t):
        for v in (t.name, lab.get(t.name, "")):
            if v and v.upper().startswith("OUTGROUP"):
                return True
        return (t.name in d["outgroup"]
                or any(o and o in t.name for o in d["outgroup"]))

    og = [t for t in root.tips() if _og(t)]
    rooted, ogstate = "", ("none" if not og else "notclade")
    if og:
        node = T.mrca(root, [t.name for t in og])
        got = {t.name for t in node.tips()} if node else set()
        if node is not None and got == {t.name for t in og}:
            root = T.reroot_on(root, node)
            rooted = ", ".join(sorted(clean(lab.get(t.name, t.name)) for t in og))
            ogstate = "rooted"
            print("   %s: rooted on %d outgroup tip(s)" % (d["stem"], len(og)))
        else:
            print("   %s: the recorded outgroup is not a clade as parsed;"
                  " left as written" % d["stem"])
    OGN = {t.name for t in og}

    sup = [n.support for n in root.walk() if n.support is not None]
    SUP01 = (max(sup) <= 1.5) if sup else False
    SUPMIN = SUPPORT_MIN_FRAC if SUP01 else SUPPORT_MIN_FRAC * 100.0

    def sfmt(v):
        return S.num(v, 2) if SUP01 else "%.0f" % v

    def kind(t):
        if ours(t.name):
            return "ours"
        return "out" if t.name in OGN else "ref"

    # ---- collapse reference-only clades, never one holding our genomes
    COLL, DEPTH, cand = {}, {}, []

    def cdepth(n):
        best = [0.0]

        def rec(x, a):
            if not x.children:
                best[0] = max(best[0], a)
            for c in x.children:
                rec(c, a + (c.length or 0.0))
        rec(n, 0.0)
        return best[0]

    def mark(n):
        tps = n.tips()
        ks = {kind(t) for t in tps}
        if n is not root and len(tps) >= MIN_COLLAPSE and ks == {"ref"}:
            cand.append(n); return
        for c in n.children:
            mark(c)

    mark(root)
    for n in cand:
        tps = n.tips()
        COLL[id(n)] = tps
        DEPTH[id(n)] = cdepth(n)
        names = [clean(lab.get(t.name, t.name)) for t in tps]
        gen = collections.Counter(x.split()[0] for x in names if x)
        n.children = []
        n.name = "%s%s (n = %d)" % (gen.most_common(1)[0][0],
                                    " κ.ά." if len(gen) > 1 else "", len(tps))
    root = T.ladderize(root)
    tips = T.layout(root)
    print("   %s: %d row(s) after collapsing %d clade(s)"
          % (d["stem"], len(tips), len(cand)))

    def pretty(t):
        if id(t) in COLL:
            return t.name
        mg = ours(t.name)
        if mg:
            c = COMP.get(mg, float("nan"))
            # two genomes are novel below FAMILY and GTDB gives them no genus.
            # Printing "?" reads as a broken lookup; naming the family and
            # saying no genus was assigned is the same fact, stated.
            who = GENUS.get(mg) or ("%s %s" % (FAM.get(mg) or "—", EL["nogen"]))
            return "%s  %s%s" % (mg.replace("ETO_", ""), who,
                                 ("  %s%%" % S.num(c, 1)) if c == c else "")
        return clean(lab.get(t.name, t.name))

    LABELS = [pretty(t) for t in tips]
    maxlab = max(tw(s, S.PT_VALUE) for s in LABELS)

    i = d["info"]
    nfall = sum(1 for t in tips if id(t) not in COLL and not ours(t.name)
                and t.name not in lab)
    cap = [
        "Δέντρο μέγιστης πιθανοφάνειας από συνενωμένα γονιδιακά σύνολα "
        "δεικτών: %s αλληλουχίες σε %s θέσεις αμινοξέων, μοντέλο %s, %s "
        "επαναλήψεις ultrafast bootstrap. Με έντονα γράμματα τα %d "
        "γονιδιώματα αυτής της εργασίας, με το ποσοστό πληρότητάς τους· τα "
        "υπόλοιπα είναι γονιδιώματα αναφοράς της GTDB."
        % (i.get("n", "—"), i.get("sites", "—"), i.get("model", "—"),
           i.get("boot", "—"), len(d["mine"])),

        ("Το δέντρο ριζώνεται στην εξωομάδα (%s)." % rooted)
        if ogstate == "rooted" else
        (("Η στοίχιση περιλαμβάνει εξωομάδα (%s), αλλά στο δέντρο αυτό οι "
          "αλληλουχίες της ΔΕΝ σχηματίζουν ενιαία ομάδα, οπότε το δέντρο "
          "σχεδιάζεται όπως το έγραψε το IQ-TREE· η ρίζωση σε μια μη "
          "μονοφυλετική εξωομάδα θα ήταν επινοημένη και η θέση της ρίζας δεν "
          "πρέπει να ερμηνευτεί."
          % ", ".join(sorted({clean(lab.get(t.name, t.name)) for t in og})))
         if ogstate == "notclade" else
         "Δεν υπάρχει εξωομάδα σε αυτό το δέντρο, οπότε σχεδιάζεται όπως το "
         "έγραψε το IQ-TREE και η θέση της ρίζας δεν πρέπει να ερμηνευτεί."),

        "ΓΙΑΤΙ ΑΥΤΟ ΤΟ ΔΕΝΤΡΟ ΚΑΙ ΟΧΙ ΤΟ ΑΛΛΟ. Το αντίστοιχο δέντρο χωρίς "
        "γονιδιώματα αναφοράς δείχνει μόνο τις σχέσεις των δικών μας "
        "γονιδιωμάτων μεταξύ τους, όχι τη θέση τους· η τοποθέτηση απαιτεί "
        "αλληλουχίες αναφοράς στην ίδια στοίχιση, και αυτές είναι εδώ.",

        ("Η υποστήριξη δίνεται ως ultrafast bootstrap (%%) και εμφανίζεται "
         "όπου είναι ≥ %s." % ("%.0f" % SUPMIN)) if sup else
        "Το αρχείο δεν φέρει τιμές υποστήριξης και δεν σχεδιάζονται.",

        "Τα τρίγωνα είναι συμπτυγμένες ομάδες αναφοράς — κάθε μέγιστη ομάδα με "
        "≥ %d φύλλα που δεν περιέχει γονιδίωμα αυτής της εργασίας· καμία "
        "τοπολογία δεν άλλαξε και δίπλα σε κάθε τρίγωνο δίνεται το πλήθος των "
        "φύλλων." % MIN_COLLAPSE,

        (("Για μία αλληλουχία αναφοράς δεν υπάρχει όνομα στα αρχεία της "
          "ανάλυσης και εμφανίζεται με τον αριθμό καταχώρισής της.")
         if nfall == 1 else
         ("Για %d αλληλουχίες αναφοράς δεν υπάρχει όνομα στα αρχεία της "
          "ανάλυσης και εμφανίζονται με τον αριθμό καταχώρισής τους." % nfall))
        if nfall else "",

        "Η θέση στο δέντρο δείχνει συγγένεια. Η ταξινομική απόδοση των "
        "γονιδιωμάτων έγινε χωριστά με το GTDB-Tk και δεν διαβάζεται από αυτό "
        "το σχήμα.",
    ]
    CAP = " ".join(x for x in cap if x)
    CW, CH = cap_geom(CAP)

    # the panel height is rows x pitch; shrink the pitch only as far as
    # 7,5 pt labels tolerate, and say so when the floor is reached
    MAX_IN, FLOOR = 5.6, 0.064
    r_ = ROW
    if len(tips) * r_ > MAX_IN:
        r_ = max(FLOOR, MAX_IN / len(tips))
        print("   %s: %d rows will not fit at the default pitch; using %.3f in"
              % (d["stem"], len(tips), r_))
    B_H = len(tips) * r_ + 0.40
    H = M_TOP + TITLE_H + GAP1 + B_H + GAP2 + LEG_H + GAP3 + CH + M_BOT
    L_IN = 0.035 * W
    R_IN = maxlab + 0.16
    fig = plt.figure(figsize=(W, H))
    fig.text(0.035, (H - M_TOP - TITLE_H * 0.5) / H, d["title"], ha="left",
             va="center", fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)
    yB = (M_BOT + CH + GAP3 + LEG_H + GAP2) / H
    ax = fig.add_axes([L_IN / W, yB, (W - L_IN - R_IN) / W, B_H / H])
    ax.axis("off")
    ax.set_ylim(len(tips) - 0.4, -0.8)
    xmax = max(t.x + DEPTH.get(id(t), 0.0) for t in tips) or 1.0
    ax.set_xlim(0, xmax * 1.02)

    COL = {"ours": S.HILITE, "ref": S.INK_2, "out": S.INK_MUTED}

    def bc(n):
        ks = {kind(t) for t in n.tips()} if n.children else {kind(n)}
        if id(n) in COLL:
            ks = {"ref"}
        return COL[list(ks)[0]] if len(ks) == 1 else S.INK

    T.draw(ax, root, color=bc, lw=0.7, support_min=SUPMIN if sup else None,
           support_fmt=sfmt, support_size=5.4, support_color=S.INK_MUTED,
           support_halo=S.SURFACE)
    for t, l in zip(tips, LABELS):
        k = "ref" if id(t) in COLL else kind(t)
        xo = t.x
        if id(t) in COLL:
            dd = max(DEPTH[id(t)], xmax * 0.012)
            ax.add_patch(Polygon([[t.x, t.y], [t.x + dd, t.y - 0.34],
                                  [t.x + dd, t.y + 0.34]], closed=True,
                                 facecolor=S.GRID, edgecolor=S.RULE, lw=0.5,
                                 zorder=2))
            xo = t.x + dd
        ax.annotate(l, (xo, t.y), textcoords="offset points", xytext=(4, 0),
                    ha="left", va="center", fontsize=S.PT_VALUE,
                    color=S.HILITE if k == "ours" else S.INK,
                    fontweight="bold" if k == "ours" else "normal",
                    zorder=4, annotation_clip=False)
    raw = xmax / 4.0
    _m = 10 ** math.floor(math.log10(raw)) if raw > 0 else 0.1
    for _f in (1.0, 2.0, 5.0, 10.0):
        step = _f * _m
        if raw <= step:
            break
    T.scale_bar(ax, 0.0, len(tips) - 0.05, step,
                "%s %s" % (S.num(step, 2), EL["scale"]), color=S.INK_2,
                size=S.PT_VALUE)

    h = [Line2D([], [], color=COL["ours"], lw=2.0, label=EL["ours"]),
         Line2D([], [], color=COL["ref"], lw=2.0, label=EL["ref"]),
         Line2D([], [], marker=">", linestyle="none", markersize=6,
                markerfacecolor=S.GRID, markeredgecolor=S.RULE,
                label=EL["coll"])]
    if OGN:
        h.insert(2, Line2D([], [], color=COL["out"], lw=2.0, label=EL["out"]))
    axL = fig.add_axes([0.035, (M_BOT + CH + GAP3) / H, 0.945, LEG_H / H])
    axL.axis("off")
    axL.legend(handles=h, loc="center", ncol=2, frameon=False, handlelength=1.4,
               borderpad=0.0, labelspacing=0.35, columnspacing=1.6,
               fontsize=S.PT_TICK)
    S.caption(fig, CAP, y=M_BOT / H, width=CW)
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    stem = "fig12a_phylogenomics_archaea" if key == "archaea" \
        else "fig12b_phylogenomics_bacteria"
    print("   %s rendered %.1f × %.1f mm" % (stem, bb.width * 25.4,
                                             bb.height * 25.4))
    if bb.width * 25.4 > FIG_W_MM + 0.5:
        print("   ! wider than the text block -- shorten the tip labels")
    S.save(fig, stem, OUT)
    plt.close(fig)

plt.close(probe)
