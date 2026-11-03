#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig13_dsr_tree.py -- where the two unbinned dsrAB operons sit, and what they are.

    python3 fig13_dsr_tree.py [workspace_2 root] [outdir]

WHAT THE FIRST VERSION GOT WRONG, AND WHY THIS ONE HAS TWO PANELS.

(1) IT DREW TWO GENES AS ONE TREE. Consecutive reference accessions sat in
    opposite halves -- AFQ44482/AFQ44483, AAL57464/AAL57463, SDX68115/SDX68081 --
    which are adjacent genes of one operon, i.e. dsrA and dsrB. The deepest
    split in that alignment is the ancient A/B PARALOGY, not a relationship
    between organisms, which is also why the reverse-type references looked
    non-monophyletic: there are two rDsr clades, one per subunit. Rooting on
    them was meaningless. Each subunit is now its own panel, and the placement
    claim is therefore made TWICE, independently -- which is stronger evidence
    than the single tree ever was.

(2) IT HAD NO CONTEXT, SO IT ARGUED IN THE CAPTION INSTEAD. A reader met a tip
    called LIQ_contig_19299_14 with no way to know what it was, where it came
    from, or why it was being examined on its own. A caption that has to argue
    for its figure is a sign the figure is not carrying its own weight. Panel
    (a) now states the facts that make the question worth asking: no MAG in the
    catalogue carries dsrAB, these two short unbinned contigs hold the only
    assembled copy, and here is how long they are and how deeply covered.

SUBUNIT ASSIGNMENT IS READ FROM THE DATA, NOT FROM THE PICTURE. Reference
deflines say "alpha subunit" or "beta subunit"; step 34's marker_ref.faa labels
our own five sequences. The tree is what those labels are used to annotate, so
taking the labels from the tree would be circular.

GENE POSITIONS ARE NOT DRAWN, BECAUSE THEY ARE NOT ON DISK. Every protein
header is a bare id with no prodigal coordinates, so panel (a) shows gene ORDER
and says so, rather than implying positions it does not have.
"""
import sys, os, re, math, collections
import etoliko_figstyle as S
import etoliko_tree as T
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle, FancyBboxPatch

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS2 = args[0] if args else "."
OUT = args[1] if len(args) > 1 else "figures"
FIG_W_MM = 160.0
SUPPORT_MIN_FRAC = 0.70      # as a FRACTION; the absolute cutoff is derived
MIN_COLLAPSE = 3

RDSR = ("Allochromatium", "Chlorobium", "Chlorobaculum", "Thiobacillus",
        "Desulfurivibrio", "Halorhodospira", "Thermochromatium", "Thiocapsa",
        "Marichromatium", "Chromatium")

EL = {
    "title": "Τα δύο μη ταξινομημένα οπερόνια dsrAB: πού βρίσκονται και τι είναι",
    "pa": "α",
    "pb": "β",
    "ha": "Τα δύο contigs που φέρουν τα μόνα συναρμολογημένα γονίδια dsrAB",
    "nomag": "σε κανένα MAG",
    "reads": "αναγνώσεις",
    "order": "σειρά γονιδίων· οι θέσεις δεν είναι διαθέσιμες",
    "ours": "αλληλουχίες αυτής της εργασίας",
    "rdsr": "αντίστροφου τύπου (rDsr) — θειοξειδωτικά",
    "red": "αναγωγικού τύπου — θειοαναγωγικά",
    "coll": "συμπτυγμένη ομάδα αναφοράς",
    "scale": "υποκ./θέση",
}

P61 = os.path.join(WS2, "58_dsr_finish")
if not os.path.isdir(P61):
    sys.exit("fig13: %s not found -- pass the workspace_2 root as the first "
             "argument" % P61)

# ---------------------------------------------------------------- the tree
TP = None
for c in ("dsr.tree", "dsr.treefile", "dsr.contree"):
    if os.path.exists(os.path.join(P61, c)):
        TP = os.path.join(P61, c); break
if TP is None:
    sys.exit("fig13: no tree file in %s" % P61)
root = T.parse(open(TP, errors="replace").read())
all_tips = root.tips()
print("   tree: %s  (%d tips)" % (TP, len(all_tips)))

# ---------------------------------------------------------------- subunits
SUB = {}
for fn in ("refs_raw.faa", "all.faa", "refs.faa", "queries.faa"):
    p = os.path.join(P61, fn)
    if not os.path.exists(p):
        continue
    for ln in open(p, errors="replace"):
        if not ln.startswith(">"):
            continue
        head = ln[1:].rstrip("\n")
        m = re.search(r"\b(alpha|beta)\s+subunit", head, re.I)
        if not m:
            continue
        tag = "dsrA" if m.group(1).lower() == "alpha" else "dsrB"
        acc = head.split()[0].split("|")[-1]
        SUB[acc] = tag
        SUB[acc.split(".")[0]] = tag
qf = os.path.join(WS2, "34_unbinned", "marker_ref.faa")
nq = 0
if os.path.exists(qf):
    for ln in open(qf, errors="replace"):
        if ln.startswith(">"):
            f = ln[1:].strip().split("|")
            if len(f) == 2 and f[0] in ("dsrA", "dsrB"):
                SUB[f[1].replace("ASM_", "")] = f[0]
                nq += 1
print("   subunit labels: %d reference(s), %d of our sequences" % (len(SUB), nq))


def acc_of(name):
    return name.split("|")[-1]


def subunit(name):
    a = acc_of(name)
    return SUB.get(a) or SUB.get(a.split(".")[0]) or "?"


def genus(name):
    f = name.split("|")
    org = f[1] if len(f) > 2 else f[-1]
    return org.split("_")[0]


def kind(name):
    if name.startswith("QUERY") or "contig_" in name:
        return "ours"
    return "rdsr" if any(genus(name).startswith(r) for r in RDSR) else "red"


# Any sequence still without a subunit gets it from its OWN BEST HIT, whose
# subunit is on a defline. LIQ_contig_19299_15 is the real case: it is dsrA's
# second ORF fragment and was found in step 59, after step 34 wrote the labels,
# so it appears in no label file. Taking its subunit from where it sits in the
# tree would be circular -- the tree is what these labels annotate -- but
# taking it from its best reference is the same evidence used for every other
# sequence.
bp = os.path.join(P61, "blastp.tsv")
inferred = []
if os.path.exists(bp):
    best = {}
    for ln in open(bp, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        try:
            bits = float(f[-1])
        except ValueError:
            continue
        if f[0] not in best or bits > best[f[0]][1]:
            best[f[0]] = (f[1], bits)
    for t in all_tips:
        if subunit(t.name) != "?":
            continue
        hit = best.get(t.name) or best.get("QUERY|" + acc_of(t.name))
        if not hit:
            continue
        a = hit[0].split("|")[-1]
        tag = SUB.get(a) or SUB.get(a.split(".")[0])
        if tag:
            SUB[acc_of(t.name)] = tag
            inferred.append((acc_of(t.name), tag, a, hit[1]))
for nm, tag, a, bits in inferred:
    print("   %s -> %s, from its best hit %s (bits %.0f)" % (nm, tag, a, bits))

unlabelled = [t.name for t in all_tips if subunit(t.name) == "?"]
by_sub = collections.Counter(subunit(t.name) for t in all_tips)
print("   tips by subunit: %s" % dict(by_sub))
if unlabelled:
    print("   %d tip(s) have no subunit in any defline; first few: %s"
          % (len(unlabelled), ", ".join(unlabelled[:4])))

# ---- are the two subunits each a clade?
# AN UNROOTED TREE HAS NO PRIVILEGED ROOT, so this test is meaningless until the
# tree is rooted somewhere. Newick writes a trifurcation and the parser has to
# pick a top; if that arbitrary top happens to fall INSIDE one subunit, that
# subunit reads as paraphyletic while the other comes out clean. That is exactly
# the asymmetry the real tree produced -- dsrA monophyletic, dsrB not -- and it
# says nothing about the sequences.
#
# The two subunit sets are complementary, so rooting on whichever one IS a clade
# necessarily makes the other one a clade as well. If after that either still
# fails, the cause is a genuinely misplaced sequence and the step refuses, since
# drawing two panels would then hide a real conflict.
def subunit_clade(rt, tag):
    names = [t.name for t in rt.tips() if subunit(t.name) == tag]
    if not names:
        return names, None, False, set()
    node = T.mrca(rt, names)
    got = {t.name for t in node.tips()} if node else set()
    return names, node, (node is not None and got == set(names)), got


info = {tag: subunit_clade(root, tag) for tag in ("dsrA", "dsrB")}
for tag in ("dsrA", "dsrB"):
    names, node, ok, got = info[tag]
    print("   %s: %d tip(s), a clade as parsed: %s" % (tag, len(names), "yes" if ok else "no"))

clean = [t for t in ("dsrA", "dsrB") if info[t][2]]
if len(clean) == 1:
    tag = clean[0]
    other = "dsrB" if tag == "dsrA" else "dsrA"
    extra = sorted(info[other][3] - set(info[other][0]))
    print("   the parsed root sits inside %s (its smallest containing clade also"
          % other)
    print("   holds %d %s sequence(s)). Rerooting on %s, which is a clade."
          % (len(extra), tag, tag))
    root = T.reroot_on(root, info[tag][1])
    all_tips = root.tips()
    info = {t: subunit_clade(root, t) for t in ("dsrA", "dsrB")}
    for t in ("dsrA", "dsrB"):
        print("   %s after rerooting: %s" % (t, "a clade" if info[t][2] else "STILL NOT a clade"))
    clean = [t for t in ("dsrA", "dsrB") if info[t][2]]

PANELS = [(t, info[t][1]) for t in ("dsrA", "dsrB") if info[t][2]]
if len(PANELS) != 2:
    bad = [t for t in ("dsrA", "dsrB") if not info[t][2]]
    for t in bad:
        stray = sorted(info[t][3] - set(info[t][0]))
        print("   ! %s is not a clade even after rerooting. Sequences of the other"
              % t)
        print("     subunit inside its smallest containing clade: %s"
              % (", ".join(x.split("|")[-1] for x in stray[:8]) or "(none)"))
    print("   ! Refusing to draw. Two panels would present a paralogy split as")
    print("     an organismal one and hide whatever is actually misplaced.")
    sys.exit(1)

# ---- support values are on whichever scale the tree writer used
# IQ-TREE writes ultrafast bootstrap as 0-100; FastTree writes SH-like local
# support as 0-1. A fixed cutoff of 70 silently removes EVERY node of a
# FastTree tree, which is what the first render did -- it showed no support at
# all while the caption said support was shown. The scale is read off the data.
sup_vals = [n.support for n in root.walk() if n.support is not None]
SUP_MAX = max(sup_vals) if sup_vals else 0.0
SUP_01 = SUP_MAX <= 1.5
SUPPORT_MIN = SUPPORT_MIN_FRAC if SUP_01 else SUPPORT_MIN_FRAC * 100.0
print("   %d support value(s), max %.3f -> %s scale; showing >= %s"
      % (len(sup_vals), SUP_MAX, "0-1" if SUP_01 else "0-100",
         S.num(SUPPORT_MIN, 2) if SUP_01 else "%.0f" % SUPPORT_MIN))
if not sup_vals:
    print("   ! no support values in this tree; the caption will say so")


def sup_fmt(v):
    return S.num(v, 2) if SUP_01 else "%.0f" % v


# ---- what can actually be claimed about distance to rDsr
def patristic(a, b):
    """path length between two tips, via their last common ancestor"""
    pa, n = {}, a
    d = 0.0
    while n is not None:
        pa[id(n)] = d
        d += (n.length or 0.0)
        n = n.parent
    n, d = b, 0.0
    while n is not None:
        if id(n) in pa:
            return pa[id(n)] + d
        d += (n.length or 0.0)
        n = n.parent
    return float("inf")


NEAREST = {}
for tag, node in PANELS:
    tps = node.tips()
    refs = [t for t in tps if kind(t.name) != "ours"]
    for q in [t for t in tps if kind(t.name) == "ours"]:
        best = min(refs, key=lambda r: patristic(q, r)) if refs else None
        if best is None:
            continue
        rd = [r for r in refs if kind(r.name) == "rdsr"]
        nearest_rdsr = min(rd, key=lambda r: patristic(q, r)) if rd else None
        NEAREST[acc_of(q.name)] = (
            tag, best.name, kind(best.name), patristic(q, best),
            nearest_rdsr.name if nearest_rdsr else None,
            patristic(q, nearest_rdsr) if nearest_rdsr else float("inf"))
print("   nearest reference to each of our sequences:")
all_red = True
for q, (tag, nm, kd, d, rn, rdd) in sorted(NEAREST.items()):
    if kd != "red":
        all_red = False
    print("     %-22s %-6s %-42s %s  (nearest rDsr %s away)"
          % (q, tag, nm.split("|")[1] if nm.count("|") > 1 else nm,
             S.num(d, 3), S.num(rdd, 3) if rdd != float("inf") else "--"))
print("   every nearest neighbour is reductive-type: %s" % ("yes" if all_red else "NO"))

# the least well-supported placement is the shortest alignment, and the figure
# should say which one rather than let every tip look equally solid
SHORTEST = None
if os.path.exists(bp):
    qb = {}
    for ln in open(bp, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        try:
            qb[acc_of(f[0])] = max(qb.get(acc_of(f[0]), 0.0), float(f[-1]))
        except ValueError:
            pass
    ours_b = {k: v for k, v in qb.items() if k in NEAREST}
    if ours_b:
        SHORTEST = min(ours_b, key=lambda k: ours_b[k])
        print("   weakest placement (lowest best-hit bit score): %s (%.0f)"
              % (SHORTEST, ours_b[SHORTEST]))

# the margin, stated as a measurement rather than as "far"
MARGIN = ""
if NEAREST and all_red:
    dn = [v[3] for v in NEAREST.values()]
    dr = [v[5] for v in NEAREST.values() if v[5] != float("inf")]
    if dn and dr:
        MARGIN = ("Ποσοτικά: η απόσταση προς τον πλησιέστερο αναγωγικού τύπου "
                  "συγγενή είναι %s–%s υποκαταστάσεις ανά θέση, ενώ προς την "
                  "πλησιέστερη αλληλουχία rDsr %s–%s — δηλαδή %s έως %s φορές "
                  "μακρύτερα."
                  % (S.num(min(dn), 3), S.num(max(dn), 3),
                     S.num(min(dr), 3), S.num(max(dr), 3),
                     S.num(min(r / n for n, r in
                               zip([v[3] for v in NEAREST.values()],
                                   [v[5] for v in NEAREST.values()])), 1),
                     S.num(max(r / n for n, r in
                               zip([v[3] for v in NEAREST.values()],
                                   [v[5] for v in NEAREST.values()])), 1)))
        print("   margin: %s" % MARGIN)

# ---------------------------------------------------------------- panel (a)
CONTIGS = collections.OrderedDict()
for t in all_tips:
    if kind(t.name) != "ours":
        continue
    m = re.match(r"^((?:LIQ|SED)_contig_\d+)_(\d+)$", acc_of(t.name))
    if m:
        CONTIGS.setdefault(m.group(1), []).append((int(m.group(2)), t.name))
for c in CONTIGS:
    CONTIGS[c].sort()

LEN, DEPTH, READS, INBIN = {}, {}, {}, {}
fa = os.path.join(WS2, "55_dsr_id", "dsr_contigs.fasta")
if os.path.exists(fa):
    cur, n = None, 0
    for ln in open(fa, errors="replace"):
        if ln.startswith(">"):
            if cur in CONTIGS:
                LEN[cur] = n
            cur, n = ln[1:].split()[0], 0
        else:
            n += len(ln.strip())
    if cur in CONTIGS:
        LEN[cur] = n
cov = os.path.join(P61, "coverage.tsv")
if os.path.exists(cov):
    for ln in open(cov, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        for c in CONTIGS:
            if c in f[0]:
                try:
                    READS[c] = int(f[3]); DEPTH[c] = float(f[6])
                except ValueError:
                    pass
bins = os.path.join(WS2, "07_mags", "bins")
if os.path.isdir(bins):
    for c in CONTIGS:
        INBIN[c] = None
        for b in os.listdir(bins):
            if not b.endswith(".fa"):
                continue
            try:
                for ln in open(os.path.join(bins, b), errors="replace"):
                    if ln.startswith(">") and ln[1:].split()[0] == c:
                        INBIN[c] = b[:-3]; break
            except OSError:
                pass
            if INBIN[c]:
                break
for c in CONTIGS:
    print("   %s: %s bp, %s reads, %s×, %s"
          % (c, LEN.get(c, "?"), READS.get(c, "?"),
             ("%.2f" % DEPTH[c]) if c in DEPTH else "?",
             INBIN.get(c) or "in no MAG"))

# ---------------------------------------------------------------- draw
S.apply()
W = FIG_W_MM / 25.4
probe = plt.figure(figsize=(W, 4.0)); probe.canvas.draw()


def tw(s, pt):
    t = probe.text(0.5, 0.5, s, fontsize=pt)
    w = t.get_window_extent(probe.canvas.get_renderer()).width / probe.dpi
    t.remove()
    return w


def pretty(t, coll):
    if id(t) in coll:
        return t.name
    f = t.name.split("|")
    if kind(t.name) == "ours":
        return acc_of(t.name)
    org = f[1] if len(f) > 2 else f[-1]
    return "%s  %s" % (org.replace("_", " "), f[-1] if len(f) > 2 else "")


# collapse per panel
PDATA = []
for tag, node in PANELS:
    coll, depth = {}, {}

    def cdepth(n):
        best = [0.0]

        def rec(x, acc):
            if not x.children:
                best[0] = max(best[0], acc)
            for c in x.children:
                rec(c, acc + (c.length or 0.0))
        rec(n, 0.0)
        return best[0]

    cand = []

    def mark(n):
        ks = {kind(t.name) for t in n.tips()}
        if n is not node and ks == {"red"} and len(n.tips()) >= MIN_COLLAPSE:
            cand.append(n); return
        for c in n.children:
            mark(c)
    mark(node)
    for n in cand:
        coll[id(n)] = n.tips()
        depth[id(n)] = cdepth(n)
        lab = collections.Counter(genus(t.name) for t in n.tips())
        n.children = []
        n.name = "%s%s (n = %d)" % (lab.most_common(1)[0][0],
                                    " κ.ά." if len(lab) > 1 else "",
                                    len(coll[id(n)]))
    node = T.ladderize(node)
    tips = T.layout(node)
    PDATA.append((tag, node, tips, coll, depth))
    print("   %s panel: %d row(s), %d clade(s) collapsed" % (tag, len(tips), len(cand)))

LABELS = {}
maxlab = 0.0
for tag, node, tips, coll, depth in PDATA:
    LABELS[tag] = [pretty(t, coll) for t in tips]
    maxlab = max(maxlab, max(tw(s, S.PT_VALUE) for s in LABELS[tag]))


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


nours = sum(len(v) for v in CONTIGS.values())
rv = sorted({genus(t.name) for t in all_tips if kind(t.name) == "rdsr"})
cap = ["(α) Τα δύο contigs που φέρουν τα μόνα συναρμολογημένα γονίδια dsrAB "
       "αυτού του συνόλου δεδομένων. Κανένα από τα 29 γονιδιώματα του καταλόγου "
       "δεν φέρει dsrA ή dsrB, και κανένα από τα δύο contigs δεν ανήκει σε "
       "κάποιο MAG — γι' αυτό εξετάζονται χωριστά και όχι ως μέρος ενός "
       "γονιδιώματος. Δίνεται το μήκος, οι αναγνώσεις και το βάθος κάλυψης. "
       "Εμφανίζεται η ΣΕΙΡΑ των γονιδίων, όχι η θέση τους: οι συντεταγμένες δεν "
       "υπάρχουν στα αρχεία της ανάλυσης και δεν επινοούνται εδώ.",
       "(β) Μέγιστη πιθανοφάνεια σε στοιχισμένες πρωτεϊνικές αλληλουχίες. Η "
       "στοίχιση περιείχε ΚΑΙ ΤΙΣ ΔΥΟ υπομονάδες, DsrA και DsrB, οι οποίες "
       "είναι παράλογα με πολύ αρχαία διπλασιασμό· επομένως η βαθύτερη διχοτόμηση "
       "του ενιαίου δέντρου διαχωρίζει ΓΟΝΙΔΙΑ και όχι οργανισμούς. Σχεδιάζονται "
       "εδώ ως δύο χωριστά δέντρα, ένα ανά υπομονάδα. Η υπομονάδα κάθε "
       "αλληλουχίας διαβάστηκε από τις περιγραφές των αρχείων αναφοράς και όχι "
       "από τη θέση της στο δέντρο.",
       "ΓΙΑΤΙ ΧΡΕΙΑΖΕΤΑΙ ΤΟ ΔΕΝΤΡΟ. Η οικογένεια dsr περιέχει και το αντίστροφο "
       "ένζυμο (rDsr) των θειοξειδωτικών βακτηρίων, ομόλογο και με υψηλή "
       "βαθμολογία σε κάθε αναζήτηση ομοιότητας. Το ότι μια αλληλουχία μοιάζει "
       "με dsrAB δεν δείχνει αναγωγή θειικών· η θέση της στο δέντρο το δείχνει. "
       "Τα γένη αντίστροφου τύπου που υπάρχουν στη στοίχιση είναι: %s."
       % ", ".join(rv),
       ("ΤΟ ΑΠΟΤΕΛΕΣΜΑ ΔΙΝΕΤΑΙ ΔΥΟ ΦΟΡΕΣ: σε κάθε δέντρο χωριστά, καμία "
        "αλληλουχία αυτής της εργασίας δεν βρίσκεται μέσα στην ομάδα rDsr, και "
        "ο πλησιέστερος συγγενής αναφοράς καθεμιάς είναι %s αναγωγικού τύπου. "
        "Σημειώνεται ότι στο δέντρο DsrA η ομάδα rDsr ΔΕΝ είναι εξωομάδα αλλά "
        "εμφωλευμένη μέσα στην αναγωγικού τύπου ποικιλότητα — αναμενόμενο, "
        "καθώς το rDsr θεωρείται ότι προέρχεται από αυτήν — οπότε η διατύπωση "
        "εδώ είναι «δεν ανήκουν στην ομάδα rDsr». %s"
        % ("σε κάθε περίπτωση" if all_red else "στις περισσότερες περιπτώσεις",
           MARGIN)),
       "Τα τρίγωνα είναι συμπτυγμένες ομάδες αναφοράς — κάθε μέγιστη ομάδα με "
       "≥ %d φύλλα χωρίς αλληλουχία αυτής της εργασίας ή αντίστροφου τύπου. Η "
       "σύμπτυξη γίνεται με βάση το περιεχόμενο, όχι την ταξινόμηση· καμία "
       "τοπολογία δεν άλλαξε και δίπλα σε κάθε τρίγωνο δίνεται το πλήθος των "
       "φύλλων." % MIN_COLLAPSE,
       ("Η τοποθέτηση δεν είναι εξίσου ισχυρή για όλες τις αλληλουχίες: %s "
        "προέρχεται από θραύσμα γονιδίου και στοιχίζεται στο μικρότερο μήκος, "
        "οπότε η θέση του είναι η λιγότερο καλά υποστηριγμένη του σχήματος."
        % SHORTEST) if SHORTEST else "",
       "Η θέση στο δέντρο δείχνει συγγένεια, όχι δραστηριότητα: το ότι ένας "
       "οργανισμός φέρει αναγωγικού τύπου dsrAB σημαίνει ότι έχει τη γενετική "
       "ικανότητα αναγωγής θειικών, όχι ότι την εκφράζει εδώ."]
CAPTION = " ".join(x for x in cap if x)
CAP_W = caption_width(probe, CAPTION)
S.caption(probe, CAPTION, width=CAP_W)
CAP_H = probe.texts[-1].get_window_extent(probe.canvas.get_renderer()).height / probe.dpi
plt.close(probe)

ROW = 0.090
A_ROW = 0.46
A_H = len(CONTIGS) * A_ROW + 0.34
M_TOP, M_BOT = 0.10, 0.10
TITLE_H, GAP1, GAP2, GAP3, SCALE_H, LEG_H = 0.26, 0.16, 0.30, 0.20, 0.30, 0.30
PH = [len(p[2]) * ROW + 0.30 for p in PDATA]
H = (M_TOP + TITLE_H + GAP1 + A_H + GAP2 + sum(PH) + GAP2
     + SCALE_H + LEG_H + GAP3 + CAP_H + M_BOT)
L_IN = 0.035 * W
R_IN = maxlab + 0.16
print("   label column %.0f mm, tree %.0f mm" % (R_IN * 25.4, (W - L_IN - R_IN) * 25.4))

fig = plt.figure(figsize=(W, H))
fig.text(0.035, (H - M_TOP - TITLE_H * 0.5) / H, EL["title"], ha="left",
         va="center", fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)

COL = {"ours": S.HILITE, "rdsr": S.DOMAIN["Archaea"], "red": S.INK_2}

# ---- panel (a): the contigs
yA = (H - M_TOP - TITLE_H - GAP1 - A_H) / H
axA = fig.add_axes([0.035, yA, 0.945, A_H / H])
axA.set_xlim(0, 1); axA.set_ylim(len(CONTIGS), -0.55); axA.axis("off")
axA.text(0.0, -0.30, "(%s)  %s" % (EL["pa"], EL["ha"]), ha="left", va="center",
         fontsize=S.PT_LABEL, color=S.INK, fontweight="bold")
maxlen = max([LEN.get(c, 1) for c in CONTIGS] or [1])
for i, (c, genes) in enumerate(CONTIGS.items()):
    y = i + 0.18
    frac = LEN.get(c, 0) / float(maxlen) if maxlen else 0.5
    bw = 0.34 * max(frac, 0.12)
    axA.add_patch(FancyBboxPatch((0.0, y - 0.07), bw, 0.14,
                                 boxstyle="round,pad=0.002,rounding_size=0.01",
                                 facecolor=S.GRID, edgecolor=S.RULE, lw=0.6))
    axA.text(0.0, y - 0.20, c, ha="left", va="center", fontsize=S.PT_VALUE,
             color=S.INK, fontweight="bold")
    facts = "%s bp · %s %s · %s× · %s" % (
        S.num(LEN.get(c, 0), 0).replace(S.DECIMAL, ""), READS.get(c, "?"),
        EL["reads"], ("%.2f" % DEPTH[c]).replace(".", S.DECIMAL) if c in DEPTH else "?",
        INBIN.get(c) or EL["nomag"])
    axA.text(0.0, y + 0.17, facts, ha="left", va="center", fontsize=S.PT_TICK,
             color=S.INK_MUTED)
    x = 0.40
    for idx, nm in genes:
        tag = subunit(nm)
        lab = "%s  %s" % (acc_of(nm).rsplit("_", 1)[-1], tag if tag != "?" else "")
        wbox = 0.085
        axA.add_patch(Rectangle((x, y - 0.075), wbox, 0.15, facecolor=S.SURFACE,
                                edgecolor=S.HILITE, lw=1.0))
        axA.text(x + wbox / 2.0, y, lab, ha="center", va="center",
                 fontsize=S.PT_VALUE, color=S.HILITE)
        x += wbox + 0.012
    # the "order not position" caveat lives in the caption, not trailing the
    # chips -- inline it ran past the right edge of the text block

# ---- panel (b): one tree per subunit
ytop = yA - GAP2 / H
for k, (tag, node, tips, coll, depth) in enumerate(PDATA):
    hh = PH[k] / H
    ytop -= hh
    ax = fig.add_axes([L_IN / W, ytop, (W - L_IN - R_IN) / W, hh])
    ax.axis("off")
    ax.set_ylim(len(tips) - 0.4, -1.25)
    # a collapsed clade's triangle extends past its node by the clade's own
    # depth, and its label starts after the triangle. The axis has to span that
    # too, or the deepest collapsed clade's label runs off the text block --
    # which is exactly what it did.
    xmax = max(t.x + depth.get(id(t), 0.0) for t in tips) or 1.0
    ax.set_xlim(0, xmax * 1.02)

    def bc(n, coll=coll):
        ks = {kind(t.name) for t in n.tips()} if n.children else {kind(n.name)}
        if id(n) in coll:
            ks = {"red"}
        if ks == {"ours"}:
            return COL["ours"]
        if ks == {"rdsr"}:
            return COL["rdsr"]
        return S.INK if len(ks) > 1 else COL["red"]

    T.draw(ax, node, color=bc, lw=0.7, support_min=SUPPORT_MIN,
           support_fmt=sup_fmt, support_size=5.6, support_color=S.INK_MUTED,
           support_halo=S.SURFACE)
    ax.text(0.0, -1.05, "(%s%d)  %s" % (EL["pb"], k + 1, tag), ha="left",
            va="center", fontsize=S.PT_LABEL, fontweight="bold", color=S.INK)
    for t, lab in zip(tips, LABELS[tag]):
        kk = "red" if id(t) in coll else kind(t.name)
        xo = t.x
        if id(t) in coll:
            d = max(depth[id(t)], xmax * 0.012)
            ax.add_patch(Polygon([[t.x, t.y], [t.x + d, t.y - 0.34],
                                  [t.x + d, t.y + 0.34]], closed=True,
                                 facecolor=S.GRID, edgecolor=S.RULE, lw=0.5,
                                 zorder=2))
            xo = t.x + d
        ax.annotate(lab, (xo, t.y), textcoords="offset points", xytext=(4, 0),
                    ha="left", va="center", fontsize=S.PT_VALUE,
                    color=S.HILITE if kk == "ours" else S.INK,
                    fontweight="bold" if kk == "ours" else "normal",
                    zorder=4, annotation_clip=False)
    # one bar PER PANEL. The two subunit trees are scaled independently -- each
    # has its own xmax -- so a single bar under the lower one would be read as
    # applying to both and would be wrong for the upper.
    step = 10 ** round(math.log10(xmax / 4.0))
    if step > xmax / 2:
        step /= 2.0
    T.scale_bar(ax, 0.0, len(tips) - 0.05, step,
                "%s %s" % (S.num(step, 2), EL["scale"]), color=S.INK_2,
                size=S.PT_VALUE)
    ytop -= GAP2 / H if k == 0 else 0.0

h = [Line2D([], [], color=COL["ours"], lw=2.0, label=EL["ours"]),
     Line2D([], [], color=COL["rdsr"], lw=2.0, label=EL["rdsr"]),
     Line2D([], [], color=COL["red"], lw=2.0, label=EL["red"]),
     Line2D([], [], marker=">", linestyle="none", markersize=6,
            markerfacecolor=S.GRID, markeredgecolor=S.RULE, label=EL["coll"])]
axL = fig.add_axes([0.035, (M_BOT + CAP_H + GAP3) / H, 0.945, LEG_H / H])
axL.axis("off")
axL.legend(handles=h, loc="center", ncol=2, frameon=False, handlelength=1.4,
           borderpad=0.0, labelspacing=0.4, columnspacing=1.6,
           fontsize=S.PT_TICK)

S.caption(fig, CAPTION, y=M_BOT / H, width=CAP_W)
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   rendered %.1f × %.1f mm" % (bb.width * 25.4, bb.height * 25.4))
if bb.width * 25.4 > FIG_W_MM + 0.5:
    print("   ! wider than the text block -- shorten the tip labels")
S.save(fig, "fig13_dsr_tree", OUT)
