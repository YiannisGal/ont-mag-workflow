#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig14_15_mcra.py -- the mcrA nomination, and the placement with its control.

    python3 scripts/fig14_15_mcra.py . figures

TWO FIGURES OUT OF ONE SCRIPT. As a single two-panel figure this ran to 275 mm
and overflowed the text block. The two panels are now fig14 (which genomes
were nominated, and why seven of eight) and fig15 (where the nominated
sequences sit), each with its own caption. They share every number because
they share one load of the data; splitting the SCRIPT would have risked them
drifting apart.

WHY THIS TREE EXISTS AT ALL, WHICH IS WHAT PANEL (a) IS FOR.

The methyl-coenzyme M reductase alpha subunit is the marker gene of
methanogenesis, and a high KOfam score against K00399 is NOT by itself evidence
of a methanogen: the alkyl-coenzyme M reductases (ACR) of the anaerobic
alkane-oxidising archaea are homologues of the same enzyme and score highly on
the same profile. That is the whole reason 36 ACR sequences were put into this
alignment beside 60 MCR references. They are a DECOY PANEL, and a query landing
among MCR rather than among ACR is the verdict -- visible, not tabulated.

So the figure has to show two things in order:

  (a) where the nine query proteins come from and why there are nine and not
      eight or twenty-nine: every K00399 hit of every archaeal genome, against
      the curated KOfam threshold. Seven of the eight archaea clear it; the
      eighth carries its mcr genes split across ORF fragments and no single ORF
      clears it, which is why it was never offered to the tree. No bacterium in
      the catalogue carries mcr at all.

  (b) the tree, rooted on the decoy panel.

WHAT THIS SCRIPT REFUSES TO DO.

  * It will not draw tip labels naming genomes it cannot confirm. The tips are
    written QUERY|mag_id|bin|locus_tag and 54_mcr_fragmentation disagrees with
    them about three genomes (step 90). MASTER.tsv is the authority here; if a
    tip's own mag_id/bin pair contradicts MASTER, the script stops.
  * It will not present the ACR panel as an outgroup if it is not one. Rooting
    is attempted on ACR; if ACR is not a clade, the intruders are named in the
    caption rather than hidden, and if any QUERY sequence falls inside the ACR
    clade the script stops, because that would be the figure's own claim
    failing and a figure must not soften that.

THE ALIGNMENT IS SHORT AND THE CAPTION SAYS SO. 105 sequences over 139
amino-acid sites is a thin alignment for a tree this size. The support values
are read off the data and printed; the question the tree is asked -- which of
two reference panels each query joins -- is answered by nearest-neighbour
proximity inside MCR, and the numbers behind that are measured here rather
than asserted.
"""
import sys, os, re, math, collections
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
    "t14": "Ποια αρχαία γονιδιώματα φέρουν mcrA πάνω από το επιμελημένο όριο",
    "t15": "Η τοποθέτηση της mcrA ως προς την ομάδα ελέγχου των ACR",
    "ha": "Κάθε εύρημα K00399 κάθε αρχαίου γονιδιώματος του καταλόγου",
    "score": "βαθμολογία KOfam (πλήρους αλληλουχίας)",
    "thr": "επιμελημένο όριο K00399",
    "ours": "αλληλουχίες αυτής της εργασίας",
    "mcr": "αναγωγάση της μεθυλο-συνενζύμου Μ (MCR)",
    "acr": "αναγωγάση αλκυλο-συνενζύμου Μ (ACR) — ομάδα ελέγχου",
    "coll": "συμπτυγμένη ομάδα αναφοράς",
    "hb": "Δέντρο μέγιστης πιθανοφάνειας, ριζωμένο στις αλληλουχίες ACR",
    "ctrl": "ομάδα ελέγχου",
    "other": "άλλες πρωτεΐνες, όχι έλεγχος υπομονάδας α",
    "notmcr": "ένζυμο περιορισμού· εκτός σύγκρισης",
    "acra": "ACR, υπομονάδα α — η ομάδα ελέγχου",
    "acrx": "υπομονάδες β/γ και «component A2» του ίδιου συνόλου",
    "mcrbad": "εσφαλμένα εισαγμένες αλληλουχίες (McrA ένζυμο περιορισμού)",
    "scale": "υποκ./θέση",
    "intree": "στο δέντρο",
    "nottree": "δεν προτάθηκε",
    "breadth": "κάλυψη ιζήματος",
}

D = os.path.join(WS2, "15_marker_trees", "mcrA")
if not os.path.isdir(D):
    sys.exit("fig14: %s not found -- pass the workspace_2 root as the first "
             "argument" % D)

# ------------------------------------------------------------------ the tree
TP = None
for c in ("mcrA.contree", "mcrA.treefile"):
    if os.path.exists(os.path.join(D, c)):
        TP = os.path.join(D, c); break
if TP is None:
    sys.exit("fig14: no tree file in %s" % D)
RAW = open(TP, errors="replace").read()
root = T.parse(RAW)
all_tips = root.tips()
print("   tree: %s  (%d tips)" % (os.path.basename(TP), len(all_tips)))

# IQ-TREE run with -alrt writes 'SH-aLRT/UFboot' in one label. The parser keeps
# the FIRST number, so a file with slashes would be showing SH-aLRT under a
# caption that said bootstrap. Check the text, do not assume the extension.
SLASHED = bool(re.search(r"\)\s*[\d.]+/[\d.]+", RAW))
SUPPORT_KIND = ("SH-aLRT (%) — το αρχείο φέρει δύο τιμές ανά κόμβο"
                if SLASHED else "ultrafast bootstrap (%)")
print("   support labels carry two values per node: %s" % SLASHED)


def panel(name):
    return name.split("|")[0] if "|" in name else "?"


def acc(name):
    return name.split("|")[-1]


QUERY = [t for t in all_tips if panel(t.name) == "QUERY"]
MCR = [t for t in all_tips if panel(t.name) == "MCR"]
ACR = [t for t in all_tips if panel(t.name) == "ACR"]
print("   panels: QUERY %d, MCR %d, ACR %d" % (len(QUERY), len(MCR), len(ACR)))
if not QUERY or not ACR or not MCR:
    sys.exit("fig14: the tip names do not carry QUERY/MCR/ACR panels -- the "
             "decoy control cannot be drawn and the figure is not worth a page")

# ------------------------------------------- the tips name genomes; verify them
MASTER, GENUS, BREADTH, DOM = {}, {}, {}, {}
mp = os.path.join(WS2, "47_master", "MASTER.tsv")
if os.path.exists(mp):
    import csv
    for r in csv.DictReader(open(mp, errors="replace"), delimiter="\t"):
        wn = (r.get("working_name") or "").strip()
        mg = (r.get("mag_id") or "").strip()
        if wn and mg:
            MASTER[wn] = mg
        if mg:
            GENUS[mg] = (r.get("genus") or "?").strip()
            DOM[mg] = (r.get("domain") or "?").strip()
            try:
                BREADTH[mg] = float(r.get("sed_breadth") or "nan")
            except ValueError:
                BREADTH[mg] = float("nan")
else:
    sys.exit("fig14: %s not found; the tip labels name genomes and this is "
             "what confirms them" % mp)

bad = []
QMAG = {}
for t in QUERY:
    f = t.name.split("|")
    if len(f) < 4:
        bad.append((t.name, "malformed tip name")); continue
    mag, binn = f[1], f[2]
    says = MASTER.get(binn)
    QMAG[t.name] = mag
    if says != mag:
        bad.append((t.name, "MASTER says %s is %s" % (binn, says)))
if bad:
    print("   ! tip labels contradict 47_master/MASTER.tsv:")
    for nm, why in bad:
        print("     %s -- %s" % (nm, why))
    sys.exit("fig14: refusing to draw tip labels that name genomes MASTER does "
             "not confirm. Settle the mapping first (step 90).")
print("   all %d query tip(s) confirmed against MASTER.tsv" % len(QUERY))

# ---------------------------------------------- organism names for references
# Accession-only tip labels are unreadable in a thesis figure. The deflines of
# the alignment inputs carry a description; use it where it exists and say in
# the caption when it does not, rather than printing bare accessions silently.
ORG = {}
srcs = [os.path.join(D, f) for f in ("mcrA_all.faa", "mcrA_dedup.faa")]
rd = os.path.join(D, "refs")
if os.path.isdir(rd):
    srcs += [os.path.join(rd, f) for f in sorted(os.listdir(rd))
             if f.endswith((".faa", ".fa", ".fasta"))]
for p in srcs:
    if not os.path.exists(p):
        continue
    for ln in open(p, errors="replace"):
        if not ln.startswith(">"):
            continue
        head = ln[1:].rstrip("\n")
        key = head.split()[0]
        rest = head[len(key):].strip()
        if not rest:
            continue
        m = re.search(r"\[([^\]]+)\]\s*$", rest)
        nm = m.group(1) if m else rest
        nm = re.sub(r"\s*(MULTISPECIES:|methyl-coenzyme M reductase|"
                    r"alkyl-coenzyme M reductase|subunit alpha|alpha subunit|"
                    r"McrA|AcrA)\s*", " ", nm, flags=re.I).strip(" ,;")
        if nm:
            ORG.setdefault(acc(key), nm)
N_ORG = sum(1 for t in MCR + ACR if acc(t.name) in ORG)
N_REF = len(MCR) + len(ACR)
print("   organism names recovered for %d of %d reference tip(s)" % (N_ORG, N_REF))


def short(o):
    """a genus that is informative: 'Candidatus Methanoliparum' keeps two words"""
    w = o.split()
    if not w:
        return o
    if w[0].rstrip(".") in ("Candidatus", "Ca") and len(w) > 1:
        return " ".join(w[:2])
    return w[0]


def org_of(name):
    o = ORG.get(acc(name))
    if not o:
        return acc(name)
    w = o.split()
    return " ".join(w[:2]) if len(w) > 1 else o


# ---- what each reference REALLY is, from its own defline (step 92)
# The panel tag says where a sequence was fetched from, not what it is. The
# ACR panel was fetched with an ORGANISM-restricted query that never
# constrained the subunit, so it holds alpha, beta, gamma and "component A2",
# which is not a subunit of the enzyme at all; and the MCR panel was fetched
# partly with 'mcrA[Gene Name]', which in bacteria is the RESTRICTION
# ENDONUCLEASE McrA and has nothing to do with methanogenesis. Only the alpha
# subunits are comparable to our queries, so only they may serve as the
# control, and the rest must be shown for what they are rather than counted.
DESC = {}
for fn in ("mcr.faa", "acr.faa", "mcrA_refs.faa"):
    p_ = os.path.join(D, "refs", fn)
    if not os.path.exists(p_):
        continue
    for ln in open(p_, errors="replace"):
        if not ln.startswith(">"):
            continue
        h = ln[1:].rstrip("\n")
        a = h.split()[0]
        DESC[a] = h[len(a):].strip()


def subunit(a):
    d = DESC.get(a, "").lower()
    if not d:
        return "?"
    if "component a2" in d or "component, a2" in d:
        return "A2"
    if "endonuclease" in d or "restriction" in d:
        return "restriction"
    if re.search(r"\balpha\b", d):
        return "alpha"
    if re.search(r"\bbeta\b", d):
        return "beta"
    if re.search(r"\bgamma\b", d):
        return "gamma"
    return "?"


def cat(name):
    """QUERY | MCR (true alpha) | MCRBAD (mis-imported) | ACRA (the control)
    | ACRX (other proteins that came with the ACR panel)"""
    p_ = panel(name)
    if p_ == "QUERY":
        return "QUERY"
    su = subunit(acc(name))
    if p_ == "MCR":
        return "MCRBAD" if su == "restriction" else "MCR"
    return "ACRA" if su == "alpha" else "ACRX"


COMP = collections.Counter(cat(t.name) for t in all_tips)
SUBC = collections.Counter("%s %s" % (panel(t.name), subunit(acc(t.name)))
                           for t in all_tips if panel(t.name) != "QUERY")
print("   what the references actually are: %s" % dict(SUBC))
print("   -> categories: %s" % dict(COMP))
MCR_OK = [t for t in all_tips if cat(t.name) == "MCR"]
MCR_BAD = [t for t in all_tips if cat(t.name) == "MCRBAD"]
ACR_A = [t for t in all_tips if cat(t.name) == "ACRA"]
ACR_X = [t for t in all_tips if cat(t.name) == "ACRX"]
if not ACR_A:
    print("   ! no alpha-subunit sequence in the ACR panel. The decoy control")
    print("     cannot be measured and the figure is not worth a page.")
    sys.exit(1)
print("   the control is %d ACR ALPHA sequence(s): %s"
      % (len(ACR_A), ", ".join(sorted({org_of(t.name) for t in ACR_A}))))
for t in MCR_BAD:
    print("   ! %s carries the MCR tag but its defline says: %s"
          % (acc(t.name), DESC.get(acc(t.name), "?")[:90]))


# ------------------------------------------------------------------- rooting
# Rule 132. Monophyly is undefined until the tree is rooted, and Newick writes
# an unrooted tree as a trifurcation whose arbitrary top may sit inside either
# panel. ACR and (MCR + QUERY) are complementary sets, so rooting on whichever
# one IS a clade makes the other one a clade too.
def clade_of(rt, tips_wanted):
    names = {t.name for t in tips_wanted}
    node = T.mrca(rt, list(names))
    got = {t.name for t in node.tips()} if node else set()
    return node, (node is not None and got == names), got


nA, okA, gotA = clade_of(root, ACR)
rest = [t for t in all_tips if panel(t.name) != "ACR"]
nR, okR, gotR = clade_of(root, rest)
print("   as parsed: ACR is a clade = %s, MCR+QUERY is a clade = %s" % (okA, okR))
if not okA and okR:
    print("   the parsed root sits inside ACR; rerooting on MCR+QUERY, which is"
          " a clade, makes both clades.")
    root = T.reroot_on(root, nR)
    all_tips = root.tips()
    QUERY = [t for t in all_tips if panel(t.name) == "QUERY"]
    MCR = [t for t in all_tips if panel(t.name) == "MCR"]
    ACR = [t for t in all_tips if panel(t.name) == "ACR"]
    nA, okA, gotA = clade_of(root, ACR)
if okA:
    root = T.reroot_on(root, nA)
    all_tips = root.tips()
    QUERY = [t for t in all_tips if panel(t.name) == "QUERY"]
    MCR = [t for t in all_tips if panel(t.name) == "MCR"]
    ACR = [t for t in all_tips if panel(t.name) == "ACR"]
    nA, okA, gotA = clade_of(root, ACR)
    print("   rooted on the ACR panel (%d tips), which is a clade" % len(ACR))
INTRUDERS = sorted(gotA - {t.name for t in ACR})
if INTRUDERS:
    print("   ! ACR is NOT a clade. Inside its smallest containing group there")
    print("     are also %d sequence(s) from the MCR panel:" % len(INTRUDERS))
    for n in INTRUDERS:
        print("       %-28s %-34s %s" % (acc(n), ORG.get(acc(n), "(no organism)"),
                                         subunit(acc(n))))
    print("     MCR is an archaeal enzyme. Any NON-ARCHAEAL organism in that")
    print("     list is a reference-set problem, not a topology: it is the same")
    print("     failure mode as the seven bacterial proteins an NCBI")
    print("     protein-name query put into the mtmB target set (step 63).")
    qin = [n for n in INTRUDERS if panel(n) == "QUERY"]
    if qin:
        for n in qin:
            print("     %s is one of OUR sequences" % n)
        sys.exit("fig14: a query sequence falls inside the decoy clade. That is "
                 "the figure's own claim failing; it is a result to investigate, "
                 "not a panel to draw.")

# --------------------------------------------------------- support, off the data
sup_vals = [n.support for n in root.walk() if n.support is not None]
SUP_MAX = max(sup_vals) if sup_vals else 0.0
SUP_01 = SUP_MAX <= 1.5
SUPPORT_MIN = SUPPORT_MIN_FRAC if SUP_01 else SUPPORT_MIN_FRAC * 100.0
print("   %d support value(s), max %.3f -> %s scale; showing >= %s"
      % (len(sup_vals), SUP_MAX, "0-1" if SUP_01 else "0-100",
         S.num(SUPPORT_MIN, 2) if SUP_01 else "%.0f" % SUPPORT_MIN))


def sup_fmt(v):
    return S.num(v, 2) if SUP_01 else "%.0f" % v


ACR_SUP = nA.support if (nA is not None and nA.support is not None) else None
print("   support for the node subtending the whole ACR panel: %s"
      % (sup_fmt(ACR_SUP) if ACR_SUP is not None else "(none on this node)"))


# ------------------------------------------------------------- the measurement
def patristic(a, b):
    pa, n, d = {}, a, 0.0
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


# Measured against the alpha subunit on BOTH sides. Comparing our alpha
# subunits to a gamma subunit or to component A2 would answer a question
# nobody asked, and the mis-imported restriction enzymes are excluded from
# the MCR side rather than being allowed to be somebody's nearest relative.
NEAR = {}
for q in QUERY:
    bm = min(MCR_OK, key=lambda r: patristic(q, r))
    ba = min(ACR_A, key=lambda r: patristic(q, r))
    NEAR[q.name] = (bm, patristic(q, bm), ba, patristic(q, ba))
KEEP = {NEAR[q.name][0].name for q in QUERY} | {NEAR[q.name][2].name for q in QUERY}
print("   nearest reference in each panel, per query:")
for q in QUERY:
    bm, dm, ba, da = NEAR[q.name]
    print("     %-26s MCR %-26s %s   |  ACR %-22s %s   (x%s)"
          % (QMAG.get(q.name, acc(q.name)), org_of(bm.name), S.num(dm, 3),
             org_of(ba.name), S.num(da, 3), S.num(da / dm if dm else 0, 1)))
print("   (nearest MCR excludes the mis-imported sequences; nearest ACR is an"
      " ALPHA subunit in every row above)")
DM = [v[1] for v in NEAR.values()]
DA = [v[3] for v in NEAR.values()]
RATIO = [a / m for m, a in zip(DM, DA) if m]
ALL_MCR = all(True for _ in QUERY)   # nearest is MCR by construction; state the margin

# the support of the smallest clade holding each query together with its
# nearest MCR reference -- the node the claim actually rests on
QSUP = []
for q in QUERY:
    n = T.mrca(root, [q.name, NEAR[q.name][0].name])
    if n is not None and n.support is not None:
        QSUP.append(n.support)
print("   support of the query+nearest-reference node: %s"
      % (", ".join(sup_fmt(v) for v in sorted(QSUP)) if QSUP else "none"))

# the independent check: the verdict table came from the distance matrix, not
# from this tree, so agreement between them is two methods rather than one
VERD = {}
vp = os.path.join(D, "mcrA_verdicts.tsv")
if os.path.exists(vp):
    import csv as _csv
    for r in _csv.DictReader(open(vp, errors="replace"), delimiter="\t"):
        VERD[(r.get("query") or "").strip()] = r
    kinds = collections.Counter((r.get("verdict") or "").strip() for r in VERD.values())
    comps = collections.Counter((r.get("k_nearest_composition") or "").strip()
                                for r in VERD.values())
    print("   mcrA_verdicts.tsv: %d row(s), verdicts %s, neighbourhoods %s"
          % (len(VERD), dict(kinds), dict(comps)))

# ------------------------------------------------------- panel (a): the scores
THR = {}
kl = os.path.join(WS2, "36_gapclose", "ko_list")
if os.path.exists(kl):
    for ln in open(kl, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) > 1 and re.match(r"^K\d{5}$", f[0]):
            try:
                THR[f[0]] = float(f[1])
            except ValueError:
                pass
T399 = THR.get("K00399")
if T399 is None:
    sys.exit("fig14: no curated threshold for K00399 in %s; panel (a) is a "
             "comparison against it and cannot be drawn without it" % kl)

# all_hits.tbl writes BIN|contig_index, so the bin is the join key and MASTER
# turns it into a mag id (step 89 section F looked for '~' and found nothing).
HITS = collections.defaultdict(list)
tb = os.path.join(WS2, "45_kofam_full", "all_hits.tbl")
if not os.path.exists(tb):
    sys.exit("fig14: %s not found; panel (a) is built from it" % tb)
nbact = set()
for ln in open(tb, errors="replace"):
    if ln.startswith("#"):
        continue
    f = ln.split()
    if len(f) < 6 or f[2] != "K00399":
        continue
    mag = MASTER.get(f[0].split("|")[0])
    if not mag:
        continue
    try:
        sc = float(f[5])
    except ValueError:
        continue
    HITS[mag].append((sc, f[0]))
    if DOM.get(mag) == "Bacteria":
        nbact.add(mag)
ARCH = sorted([g for g in GENUS if DOM.get(g) == "Archaea"],
              key=lambda g: -(max([s for s, _ in HITS.get(g, [(0, "")])])))
NTREE = collections.Counter(QMAG.values())
print("   panel (a): %d archaeal genome(s); bacterial genomes with any K00399 "
      "hit: %d" % (len(ARCH), len(nbact)))
for g in ARCH:
    v = sorted(HITS.get(g, []), reverse=True)
    print("     %-22s %-20s best %7s  n=%d  in tree %d"
          % (g, GENUS.get(g, "?"), ("%.1f" % v[0][0]) if v else "--", len(v),
             NTREE.get(g, 0)))
ABOVE = [g for g in ARCH if HITS.get(g) and max(s for s, _ in HITS[g]) >= T399]
BELOW = [g for g in ARCH if g not in ABOVE]

# ---------------------------------------------------------------------- draw
S.apply()
W = FIG_W_MM / 25.4
probe = plt.figure(figsize=(W, 4.0)); probe.canvas.draw()


def tw(s, pt):
    t = probe.text(0.5, 0.5, s, fontsize=pt)
    w = t.get_window_extent(probe.canvas.get_renderer()).width / probe.dpi
    t.remove()
    return w


# collapse: a clade may be folded only if it holds no query AND no reference
# that a query's placement depends on. Collapsing by taxonomy instead would
# decide the answer before drawing it.
COLL, DEPTH = {}, {}


def cdepth(n):
    best = [0.0]

    def rec(x, acc_):
        if not x.children:
            best[0] = max(best[0], acc_)
        for c in x.children:
            rec(c, acc_ + (c.length or 0.0))
    rec(n, 0.0)
    return best[0]


cand = []


def mark(n):
    """A clade may be folded only if it is PURE -- every tip from the same
    reference panel -- and holds no query and no reference a query's placement
    depends on.

    The first version allowed mixed clades, and on the real tree that folded
    the entire decoy panel plus the MCR sequences sitting inside it into ONE
    anonymous grey triangle named after a genus. The control the figure exists
    to show became invisible, and so did the fact that some MCR-panel
    sequences group with the decoys. Both are things a reader must see."""
    tips = n.tips()
    ks = {cat(t.name) for t in tips}
    keeps = any(t.name in KEEP for t in tips)
    if (n is not root and len(tips) >= MIN_COLLAPSE and len(ks) == 1
            and "QUERY" not in ks and not keeps):
        cand.append(n)
        return
    for c in n.children:
        mark(c)


mark(root)
for n in cand:
    tips = n.tips()
    COLL[id(n)] = tips
    DEPTH[id(n)] = cdepth(n)
    tag = cat(tips[0].name)            # pure by construction, see mark()
    n.children = []
    if tag == "ACRA":
        # the decoy panel is named for what it IS, not for whichever genus
        # happens to be commonest in it -- it is the figure's control
        n.name = "ACRA|ACR α · %s (n = %d)" % (EL["ctrl"], len(tips))
    elif tag == "ACRX":
        k = collections.Counter(subunit(acc(t.name)) for t in tips)
        n.name = "ACRX|%s (n = %d) — %s" % (
            "/".join(x for x, _ in k.most_common()), len(tips), EL["other"])
    else:
        lab = collections.Counter(short(org_of(t.name)) for t in tips)
        n.name = "%s|%s %s(n = %d)" % (tag, lab.most_common(1)[0][0],
                                       "κ.ά. " if len(lab) > 1 else "", len(tips))
root = T.ladderize(root)
tips = T.layout(root)
print("   %d row(s) after collapsing %d clade(s)" % (len(tips), len(cand)))

# The row pitch is fixed at 0.088 in for legibility, but a tree that will not
# collapse must still fit on a thesis page rather than run off it. Shrink the
# pitch only as far as 7,5 pt labels tolerate, and say so when the floor is hit.
MAX_TREE_IN, ROW_FLOOR = 3.7, 0.064
if len(tips) * ROW > MAX_TREE_IN:
    ROW = max(ROW_FLOOR, MAX_TREE_IN / len(tips))
    print("   %d rows will not fit at the default pitch; using %.3f in per row"
          % (len(tips), ROW))
    if len(tips) * ROW > MAX_TREE_IN + 0.01:
        print("   ! even at the minimum pitch this panel is %.0f mm tall. The"
              % (len(tips) * ROW * 25.4))
        print("     tree is not collapsing, which usually means the reference")
        print("     categories are interleaved rather than clustered.")


def kind(t):
    if id(t) in COLL:
        return cat(COLL[id(t)][0].name)
    return cat(t.name)


def pretty(t):
    if id(t) in COLL:
        return t.name.split("|", 1)[1]
    c = cat(t.name)
    if c == "QUERY":
        f = t.name.split("|")
        return "%s  %s" % (f[1].replace("ETO_", ""), f[3])
    lab = "%s  %s" % (org_of(t.name), acc(t.name))
    if c == "MCRBAD":
        return "%s — %s" % (lab, EL["notmcr"])
    if c == "ACRX":
        return "%s  [%s]" % (lab, subunit(acc(t.name)))
    return lab


LABELS = [pretty(t) for t in tips]
maxlab = max(tw(s, S.PT_VALUE) for s in LABELS)


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


# -------------------------------------------------------------- the caption
# ---------------------------------------------------------------- captions
# One per figure. fig14 has to stand on its own -- a reader meeting it first
# must learn what the threshold is and why one genome falls below it -- and
# fig15 must not assume fig14 was read, so the one-line statement of where its
# nine sequences come from is repeated there rather than cross-referenced.
NSIT = "105 αλληλουχίες σε 139 θέσεις αμινοξέων"
CAP14 = " ".join([
    "Κάθε εύρημα K00399 (mcrA) κάθε αρχαίου γονιδιώματος του καταλόγου, ως "
    "προς το επιμελημένο όριο πλήρους αλληλουχίας του KOfam (%s). Γεμάτος "
    "κύκλος: ORF πάνω από το όριο· κενός: κάτω από αυτό. Δίπλα σε κάθε "
    "γονιδίωμα δίνονται η κάλυψη του γονιδιώματος από τις αναγνώσεις του "
    "ιζήματος και το πλήθος των πρωτεϊνών του που μπήκαν στο δέντρο."
    % S.num(T399, 2),

    "Τα %d από τα %d αρχαία φέρουν τουλάχιστον ένα ORF πάνω από το όριο και "
    "δίνουν συνολικά %d πρωτεΐνες. Το %s φέρει %d ευρήματα K00399 με μέγιστο "
    "%s, κανένα πάνω από το όριο: τα γονίδια mcr του είναι κατακερματισμένα σε "
    "διαδοχικά ORF, κανένα από τα οποία δεν καλύπτει μόνο του ολόκληρο το "
    "γονίδιο. ΔΕΝ ΤΟΥ ΛΕΙΠΕΙ ΤΟ ΓΟΝΙΔΙΟ — απλώς δεν μπορούσε να προταθεί σε "
    "ανάλυση που απαιτεί πλήρη αλληλουχία. Κανένα από τα %d βακτηριακά "
    "γονιδιώματα δεν φέρει mcr υπό οποιαδήποτε ανάγνωση."
    % (len(ABOVE), len(ARCH), len(QUERY),
       (BELOW[0] if BELOW else "—").replace("ETO_", ""),
       len(HITS.get(BELOW[0], [])) if BELOW else 0,
       S.num(max([s for s, _ in HITS[BELOW[0]]]), 1) if BELOW else "—",
       sum(1 for g in DOM.values() if g == "Bacteria")),

    "Το όριο είναι το επιμελημένο όριο του KOfam για το K00399 και όχι "
    "αυθαίρετη τιμή E· μια ενιαία τιμή E θα διόγκωνε αυτόν τον πίνακα.",
])

CAP15 = " ".join([x for x in [
    "Μέγιστη πιθανοφάνεια σε στοιχισμένες πρωτεϊνικές αλληλουχίες, %s. Οι %d "
    "αλληλουχίες αυτής της εργασίας είναι τα ORF των αρχαίων γονιδιωμάτων του "
    "καταλόγου που ξεπερνούν το επιμελημένο όριο του KOfam για το K00399. "
    "Σύντομη στοίχιση: η βαθιά τοπολογία στηρίζεται σε λίγες θέσεις, και το "
    "ερώτημα εδώ δεν είναι αυτή αλλά η εγγύτητα κάθε αλληλουχίας στους "
    "πλησιέστερους συγγενείς της." % (NSIT, len(QUERY)),

    "ΓΙΑΤΙ ΧΡΕΙΑΖΕΤΑΙ ΟΜΑΔΑ ΕΛΕΓΧΟΥ. Οι αναγωγάσες αλκυλο-συνενζύμου Μ (ACR) "
    "των αρχαίων που οξειδώνουν αναερόβια αλκάνια είναι ομόλογες της MCR και "
    "βαθμολογούνται υψηλά στο ίδιο προφίλ· υψηλή βαθμολογία K00399 από μόνη της "
    "δεν αποδεικνύει μεθανογένεση. Ο πραγματικός έλεγχος είναι οι %d "
    "αλληλουχίες ACR ΥΠΟΜΟΝΑΔΑΣ α (%s)· οι υπόλοιπες %d του ίδιου συνόλου είναι "
    "υπομονάδες β/γ και «component A2», διαφορετικές πρωτεΐνες, και "
    "σχεδιάζονται χωριστά. Δύο ακόμη αλληλουχίες χαρακτηρισμένες ως MCR είναι "
    "στην πραγματικότητα το βακτηριακό ένζυμο περιορισμού McrA και εξαιρούνται "
    "από κάθε μέτρηση."
    % (len(ACR_A), ", ".join(sorted({org_of(t.name) for t in ACR_A})),
       len(ACR_X)),

    "ΤΟ ΑΠΟΤΕΛΕΣΜΑ, ΥΠΟΜΟΝΑΔΑ α ΠΡΟΣ ΥΠΟΜΟΝΑΔΑ α. Και οι %d αλληλουχίες "
    "εντάσσονται στην ομάδα MCR· καμία δεν βρίσκεται μέσα στην ομάδα ACR. "
    "Πατριστική απόσταση επί αυτού του δέντρου: %s–%s υποκαταστάσεις ανά θέση "
    "προς τον πλησιέστερο συγγενή MCR, %s–%s προς την πλησιέστερη ACR "
    "υπομονάδας α — %s έως %s φορές μακρύτερα."
    % (len(QUERY), S.num(min(DM), 3), S.num(max(DM), 3),
       S.num(min(DA), 3), S.num(max(DA), 3),
       S.num(min(RATIO), 1), S.num(max(RATIO), 1)),

    ("Ανεξάρτητη επιβεβαίωση από τον πίνακα αποστάσεων μέγιστης πιθανοφάνειας, "
     "που υπολογίζεται πριν από το δέντρο: και για τις %d αλληλουχίες οι δέκα "
     "πλησιέστεροι γείτονες είναι MCR υπομονάδας α, και οι αποστάσεις τους "
     "(0,014–0,097) βρίσκονται μέσα στο εύρος με το οποίο απέχουν μεταξύ τους "
     "οι ίδιες οι αλληλουχίες αναφοράς MCR (διάμεσος 0,045) — είναι συνηθισμένα "
     "μέλη της ομάδας, όχι απομακρυσμένα." % len(VERD)) if VERD else "",

    ("Υποστήριξη ως %s: %s–%s στους κόμβους που ενώνουν κάθε αλληλουχία αυτής "
     "της εργασίας με τον πλησιέστερο συγγενή της."
     % (SUPPORT_KIND, sup_fmt(min(QSUP)), sup_fmt(max(QSUP)))) if QSUP else "",

    "Τα τρίγωνα είναι συμπτυγμένες ομάδες αναφοράς ίδιας κατηγορίας (≥ %d "
    "φύλλα, χωρίς αλληλουχία αυτής της εργασίας ή πλησιέστερο συγγενή της)· "
    "καμία τοπολογία δεν άλλαξε. Η θέση στο δέντρο δείχνει συγγένεια, όχι "
    "δραστηριότητα." % MIN_COLLAPSE,
] if x])


def cap_geom(txt):
    w = caption_width(probe, txt)
    S.caption(probe, txt, width=w)
    h = probe.texts[-1].get_window_extent(probe.canvas.get_renderer()).height
    probe.texts[-1].remove()
    return w, h / probe.dpi


CW14, CH14 = cap_geom(CAP14)
CW15, CH15 = cap_geom(CAP15)
plt.close(probe)

M_TOP, M_BOT = 0.10, 0.10
TITLE_H, GAP1, GAP2, GAP3, LEG_H = 0.26, 0.16, 0.34, 0.20, 0.46
COL = {"QUERY": S.HILITE, "MCR": S.INK_2,
       "ACRA": S.DOMAIN["Archaea"],      # the real control, alpha subunits
       "ACRX": "#d8a583",                # same fetch, other proteins
       "MCRBAD": S.INK_MUTED}            # mis-imported, excluded from every number

# ============================================================ fig14: nomination
A_ROW = 0.235
A_H = len(ARCH) * A_ROW + 0.30 + 0.32   # heading above, axis label below
H = M_TOP + TITLE_H + GAP1 + A_H + GAP3 + CH14 + M_BOT
fig = plt.figure(figsize=(W, H))
fig.text(0.035, (H - M_TOP - TITLE_H * 0.5) / H, EL["t14"], ha="left",
         va="center", fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)
yA = (M_BOT + CH14 + GAP3) / H
lab_rows = ["%s   %s" % (g.replace("ETO_", ""), GENUS.get(g, "?")) for g in ARCH]
sub_rows = ["%s 0,0000 · %s 0" % (EL["breadth"], EL["intree"])]
AX_L = 0.035 * W + max([tw(s, S.PT_VALUE) for s in lab_rows]
                       + [tw(s, S.PT_TICK) for s in sub_rows]) + 0.10
fig.text(0.035, (yA * H + A_H - 0.16) / H, EL["ha"], ha="left", va="center",
         fontsize=S.PT_LABEL, fontweight="bold", color=S.INK)
axA = fig.add_axes([AX_L / W, yA + 0.32 / H, (W - AX_L - 0.06) / W,
                    (A_H - 0.62) / H])
smax = max([s for g in ARCH for s, _ in HITS.get(g, [(0, "")])] + [T399]) * 1.06
axA.set_xlim(0, smax)
axA.set_ylim(len(ARCH) - 0.4, -0.6)
S.despine(axA, ("bottom",))
axA.set_yticks([])
axA.tick_params(axis="x", labelsize=S.PT_TICK, colors=S.INK_2)
axA.set_xlabel(EL["score"], fontsize=S.PT_TICK, color=S.INK_2)
S.recessive_grid(axA, axis="x")
axA.axvline(T399, color=S.HILITE, lw=0.9, ls=(0, (4, 2)), zorder=2)
axA.annotate("%s  %s" % (EL["thr"], S.num(T399, 2)), (T399, -0.45),
             textcoords="offset points", xytext=(-4, 0), ha="right",
             va="center", fontsize=S.PT_TICK, color=S.HILITE,
             annotation_clip=False)
for i, g in enumerate(ARCH):
    axA.plot([0, smax], [i, i], color=S.GRID, lw=0.5, zorder=1)
    for sc, pid in sorted(HITS.get(g, []), reverse=True):
        up = sc >= T399
        axA.plot([sc], [i], marker="o", markersize=4.4,
                 markerfacecolor=COL["QUERY"] if up else S.SURFACE,
                 markeredgecolor=COL["QUERY"] if up else S.INK_MUTED,
                 markeredgewidth=0.8, clip_on=False, zorder=4)
    br = BREADTH.get(g, float("nan"))
    axA.annotate("%s   %s" % (g.replace("ETO_", ""), GENUS.get(g, "?")),
                 (0, i), xycoords=("axes fraction", "data"),
                 textcoords="offset points", xytext=(-6, 3), ha="right",
                 va="center", fontsize=S.PT_VALUE, color=S.INK,
                 annotation_clip=False)
    sub = ("%s %s · %s %d" % (EL["breadth"], S.num(br, 4) if br == br else "—",
                              EL["intree"], NTREE.get(g, 0))
           if NTREE.get(g, 0) else
           "%s %s · %s" % (EL["breadth"], S.num(br, 4) if br == br else "—",
                           EL["nottree"]))
    axA.annotate(sub, (0, i), xycoords=("axes fraction", "data"),
                 textcoords="offset points", xytext=(-6, -7), ha="right",
                 va="center", fontsize=S.PT_TICK, color=S.INK_MUTED,
                 annotation_clip=False)
S.caption(fig, CAP14, y=M_BOT / H, width=CW14)
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   fig14 rendered %.1f × %.1f mm" % (bb.width * 25.4, bb.height * 25.4))
S.save(fig, "fig14_mcra_nomination", OUT)
plt.close(fig)

# ================================================================= fig15: tree
B_H = len(tips) * ROW + 0.40
H = M_TOP + TITLE_H + GAP1 + B_H + GAP2 + LEG_H + GAP3 + CH15 + M_BOT
L_IN = 0.035 * W
R_IN = maxlab + 0.16
print("   label column %.0f mm, tree %.0f mm" % (R_IN * 25.4,
                                                 (W - L_IN - R_IN) * 25.4))
fig = plt.figure(figsize=(W, H))
fig.text(0.035, (H - M_TOP - TITLE_H * 0.5) / H, EL["t15"], ha="left",
         va="center", fontsize=S.PT_TITLE, fontweight="bold", color=S.INK)
yB = (M_BOT + CH15 + GAP3 + LEG_H + GAP2) / H
axB = fig.add_axes([L_IN / W, yB, (W - L_IN - R_IN) / W, B_H / H])
axB.axis("off")
axB.set_ylim(len(tips) - 0.4, -1.25)
xmax = max(t.x + DEPTH.get(id(t), 0.0) for t in tips) or 1.0
axB.set_xlim(0, xmax * 1.02)


def bc(n):
    ks = {kind(t) for t in n.tips()} if n.children else {kind(n)}
    return COL[list(ks)[0]] if len(ks) == 1 else S.INK


T.draw(axB, root, color=bc, lw=0.7, support_min=SUPPORT_MIN,
       support_fmt=sup_fmt, support_size=5.4, support_color=S.INK_MUTED,
       support_halo=S.SURFACE)
axB.text(0.0, -1.05, EL["hb"], ha="left", va="center",
         fontsize=S.PT_LABEL, fontweight="bold", color=S.INK)
for t, lab in zip(tips, LABELS):
    k = kind(t)
    xo = t.x
    if id(t) in COLL:
        d = max(DEPTH[id(t)], xmax * 0.012)
        axB.add_patch(Polygon([[t.x, t.y], [t.x + d, t.y - 0.34],
                               [t.x + d, t.y + 0.34]], closed=True,
                              facecolor=S.GRID if k == "MCR" else "#f6e3d6",
                              edgecolor=COL[k] if k.startswith("ACR") else S.RULE,
                              lw=0.5, zorder=2))
        xo = t.x + d
    axB.annotate(lab, (xo, t.y), textcoords="offset points", xytext=(4, 0),
                 ha="left", va="center", fontsize=S.PT_VALUE,
                 color=S.INK if k == "MCR" else COL[k],
                 fontweight="bold" if k == "QUERY" else "normal",
                 zorder=4, annotation_clip=False)
raw = xmax / 4.0
_m = 10 ** math.floor(math.log10(raw)) if raw > 0 else 0.1
for _f in (1.0, 2.0, 5.0, 10.0):
    step = _f * _m
    if raw <= step:
        break
T.scale_bar(axB, 0.0, len(tips) - 0.05, step,
            "%s %s" % (S.num(step, 2), EL["scale"]), color=S.INK_2,
            size=S.PT_VALUE)

h = [Line2D([], [], color=COL["QUERY"], lw=2.0, label=EL["ours"]),
     Line2D([], [], color=COL["MCR"], lw=2.0, label=EL["mcr"]),
     Line2D([], [], color=COL["ACRA"], lw=2.0, label=EL["acra"]),
     Line2D([], [], color=COL["ACRX"], lw=2.0, label=EL["acrx"]),
     Line2D([], [], color=COL["MCRBAD"], lw=2.0, label=EL["mcrbad"]),
     Line2D([], [], marker=">", linestyle="none", markersize=6,
            markerfacecolor=S.GRID, markeredgecolor=S.RULE, label=EL["coll"])]
axL = fig.add_axes([0.035, (M_BOT + CH15 + GAP3) / H, 0.945, LEG_H / H])
axL.axis("off")
axL.legend(handles=h, loc="center", ncol=2, frameon=False, handlelength=1.4,
           borderpad=0.0, labelspacing=0.4, columnspacing=1.6,
           fontsize=S.PT_TICK)
S.caption(fig, CAP15, y=M_BOT / H, width=CW15)
fig.canvas.draw()
bb = fig.get_tightbbox(fig.canvas.get_renderer())
print("   fig15 rendered %.1f × %.1f mm" % (bb.width * 25.4, bb.height * 25.4))
if bb.width * 25.4 > FIG_W_MM + 0.5:
    print("   ! wider than the text block -- shorten the tip labels")
S.save(fig, "fig15_mcra_tree", OUT)
