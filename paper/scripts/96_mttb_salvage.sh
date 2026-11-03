#!/usr/bin/env bash
# 96_mttb_salvage.sh -- the mttB matrix is broken; is the QUESTION still
# answerable, and what survives in the meantime?
#
# WHAT STEP 95 ESTABLISHED, AND IT IS NOT A SUBTLE RESULT.
#
# Every one of the five mttB queries shares **ZERO** alignment columns with the
# reference the distance matrix calls its nearest, in BOTH panels. A maximum
# likelihood distance is estimated from the columns where both sequences carry
# a residue; with none, there is nothing to estimate, IQ-TREE writes its floor
# of 1e-6, and that floor then sorts FIRST when you take the minimum of a row.
#
#     "nearest neighbour" in that matrix means "most undefined".
#
# So the k-nearest-neighbour composition behind every CONFIRMED verdict in
# mttB_verdicts.tsv selects the LEAST comparable references, not the most.
# The verdicts are not merely weak, they are inverted.
#
# mcrA in the same run is the positive control and passes cleanly: 139 of 139
# shared columns for every query against its nearest reference, a minimum of 52
# against ANY reference, and 4 floor values in 5460 pairs (0,1%) against mttB's
# 96 in 3003 (3,2%).
#
# BUT ONE NUMBER IN STEP 95 SAYS THE QUESTION MAY STILL BE ANSWERABLE. Each
# query's shared-column count against ALL references runs
#
#     min 0,  median 128,  max 320
#
# and 320 is every residue the query has. So SOME references do overlap our
# sequences completely -- they are simply not the ones the broken minimum
# picks. Section A ranks every reference by shared columns and names the ones
# at the top. If those are the archaeal trimethylamine methyltransferases
# (Methanococcoides), the panel can be rebuilt and the tree redone. If even
# they share nothing, the problem is not the reference set and these five
# proteins need re-examining before any tree is drawn.
#
# Section B gives the fallback claim that survives either way: the curated
# KOfam score per genome for K14083, which needs no alignment at all.
#
# Section C closes a loose end. Step 95 section D found mcrA's verdict
# distances ARE patristic (9 of 9 agree) but mttB's off by a CONSTANT 0,0033 in
# all five rows. A constant offset is not a different metric, it is the same
# metric on a slightly different tree, and .treefile and .contree share a
# topology while differing in branch lengths. Tested here rather than assumed.
set -u
WS2="${WS2}"
OUT="$WS2/96_mttb_salvage"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re, csv, collections, statistics
WS2, OUT = sys.argv[1], sys.argv[2]
D = os.path.join(WS2, "15_marker_trees", "mttB")


def read_fa(p):
    seqs, name, buf = {}, None, []
    for ln in open(p, errors="replace"):
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].strip(), []
        else:
            buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


DESC = {}
rd = os.path.join(D, "refs")
if os.path.isdir(rd):
    for fn in sorted(os.listdir(rd)):
        if not fn.endswith((".faa", ".fa", ".fasta")):
            continue
        for ln in open(os.path.join(rd, fn), errors="replace"):
            if ln.startswith(">"):
                h = ln[1:].rstrip("\n")
                a = h.split()[0]
                DESC[a] = h[len(a):].strip()
                DESC[a.split("|")[-1]] = DESC[a]


def org(t):
    d = DESC.get(t.split("|")[-1], "") or DESC.get(t, "")
    m = re.search(r"\[([^\]]+)\]", d)
    if m:
        return m.group(1)
    m = re.search(r"_([A-Z]{4,6})\b", d)
    return d[:38] if d else t.split("|")[-1]


ARCH = ("methano", "methermicoccus", "archaeoglob", "halo", "thermoplasma")


def is_arch(o):
    o = (o or "").lower()
    return any(a in o for a in ARCH)


print("=" * 64)
print("A. which references actually overlap our sequences")
print("=" * 64)
aln = os.path.join(D, "mttB_trim.faa")
if not os.path.exists(aln):
    print("   MISSING %s" % aln); sys.exit(0)
A = read_fa(aln)
L = len(next(iter(A.values())))
GAPS = set("-.")
MASK = {k: [c not in GAPS for c in v] for k, v in A.items()}
qs = [k for k in A if k.startswith("QUERY")]
refs = [k for k in A if not k.startswith("QUERY")]
print("   %d columns; %d queries, %d references" % (L, len(qs), len(refs)))
print()


def span(k):
    m = MASK[k]
    first = next((i for i, x in enumerate(m) if x), -1)
    last = len(m) - 1 - next((i for i, x in enumerate(reversed(m)) if x), 0)
    return first, last, sum(m)


print("   the column block each sequence occupies (first..last, residues):")
for k in sorted(qs):
    f, l, n = span(k)
    print("      %-52s %4d..%-4d  %4d" % (k.split("|")[1] + " " + k.split("|")[3],
                                          f, l, n))
sp = [span(r) for r in refs]
print("      references: first column min %d median %d max %d; last column"
      " min %d median %d max %d"
      % (min(x[0] for x in sp), int(statistics.median([x[0] for x in sp])),
         max(x[0] for x in sp), min(x[1] for x in sp),
         int(statistics.median([x[1] for x in sp])), max(x[1] for x in sp)))
print()

BEST = {}
for q in sorted(qs):
    ov = sorted(((sum(1 for a, b in zip(MASK[q], MASK[r]) if a and b), r)
                 for r in refs), reverse=True)
    BEST[q] = ov
    f = q.split("|")
    print("   %s %s -- the 10 references sharing the most columns:"
          % (f[1].replace("ETO_", ""), f[3]))
    print("      %6s  %-12s %-34s %s" % ("shared", "panel", "organism", "archaeal"))
    for n, r in ov[:10]:
        o = org(r)
        print("      %6d  %-12s %-34s %s"
              % (n, r.split("|")[0], o[:34], "yes" if is_arch(o) else ""))
    narch = sum(1 for n, r in ov if n >= 100 and is_arch(org(r)))
    print("      references sharing >= 100 columns: %d, of which archaeal: %d"
          % (sum(1 for n, _ in ov if n >= 100), narch))
    print()

with open(os.path.join(OUT, "mttB_overlap_full.tsv"), "w") as fh:
    fh.write("query\treference\tpanel\torganism\tshared_columns\n")
    for q in sorted(qs):
        for n, r in BEST[q]:
            fh.write("%s\t%s\t%s\t%s\t%d\n" % (q, r, r.split("|")[0], org(r), n))
print("   -> 96_mttb_salvage/mttB_overlap_full.tsv")

# ---- a defensible nearest neighbour: only pairs with real overlap
md = os.path.join(D, "mttB.mldist")
if os.path.exists(md):
    toks = open(md, errors="replace").read().split()
    n = int(toks[0]); i = 1
    names, rows = [], []
    for _ in range(n):
        names.append(toks[i]); i += 1
        v = []
        while len(v) < n:
            v.append(float(toks[i])); i += 1
        rows.append(v)
    idx = {x: j for j, x in enumerate(names)}
    print()
    print("   nearest reference RESTRICTED to pairs sharing >= 100 columns:")
    print("   %-26s %-12s %9s  %-30s %s"
          % ("query", "panel", "distance", "organism", "archaeal"))
    for q in sorted(qs):
        if q not in idx:
            print("   %-26s not in the distance matrix (folded as a duplicate)"
                  % q.split("|")[1]); continue
        for pn in ("TMA_MttB", "BETAINE_SF"):
            cand = [(rows[idx[q]][idx[r]], r) for nsh, r in BEST[q]
                    if nsh >= 100 and r in idx and r.split("|")[0] == pn]
            if not cand:
                print("   %-26s %-12s  no reference in this panel shares >= 100"
                      " columns" % (q.split("|")[1].replace("ETO_", ""), pn))
                continue
            dist, who = min(cand)
            o = org(who)
            print("   %-26s %-12s %9.4f  %-30s %s"
                  % (q.split("|")[1].replace("ETO_", ""), pn, dist, o[:30],
                     "yes" if is_arch(o) else ""))
        print()

print("=" * 64)
print("B. what survives with no alignment at all: the KOfam scores")
print("=" * 64)
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
MASTER, GEN, DOM = {}, {}, {}
mp = os.path.join(WS2, "47_master", "MASTER.tsv")
for r in csv.DictReader(open(mp, errors="replace"), delimiter="\t"):
    wn = (r.get("working_name") or "").strip(); mg = (r.get("mag_id") or "").strip()
    if wn and mg:
        MASTER[wn] = mg
    if mg:
        GEN[mg] = (r.get("genus") or "?").strip()
        DOM[mg] = (r.get("domain") or "?").strip()
KOS = ("K14083", "K16176", "K16178", "K04480")
NAMES = {"K14083": "mttB trimethylamine MT", "K16176": "mtmB monomethylamine MT",
         "K16178": "mtbB dimethylamine MT", "K04480": "mtaB methanol MT"}
per = collections.defaultdict(lambda: collections.defaultdict(list))
tb = os.path.join(WS2, "45_kofam_full", "all_hits.tbl")
if os.path.exists(tb):
    for ln in open(tb, errors="replace"):
        if ln.startswith("#"):
            continue
        f = ln.split()
        if len(f) < 6 or f[2] not in KOS:
            continue
        mag = MASTER.get(f[0].split("|")[0])
        if not mag:
            continue
        try:
            per[mag][f[2]].append((float(f[5]), f[0]))
        except ValueError:
            pass
else:
    print("   MISSING %s" % tb)
for ko in KOS:
    t = THR.get(ko)
    hits = [(max(s for s, _ in per[g][ko]), g) for g in per if per[g].get(ko)]
    above = [(s, g) for s, g in hits if t and s >= t]
    print()
    print("   %s  %s   curated threshold %s" % (ko, NAMES[ko],
                                                ("%.2f" % t) if t else "?"))
    print("      %d genome(s) with any hit, %d above threshold"
          % (len(hits), len(above)))
    for s, g in sorted(hits, reverse=True):
        print("      %-22s %-20s %-10s %8.1f  %s"
              % (g, GEN.get(g, "?")[:20], DOM.get(g, "?"), s,
                 "ABOVE" if (t and s >= t) else "below"))

print()
print("=" * 64)
print("C. the constant 0,0033 offset in mttB's verdict distances")
print("=" * 64)


class N(object):
    __slots__ = ("name", "length", "children", "parent")

    def __init__(s, nm="", ln=0.0):
        s.name = nm; s.length = ln; s.children = []; s.parent = None

    def tips(s):
        if not s.children:
            return [s]
        o = []
        for c in s.children:
            o.extend(c.tips())
        return o


def parse(text):
    s = text.strip()
    if ";" in s:
        s = s[:s.rindex(";")]
    root = N(); cur = root; buf = ""; closed = None

    def lab(tok):
        nm, _, ln = tok.rpartition(":") if ":" in tok else (tok, "", "")
        try:
            x = float(ln)
        except ValueError:
            x = 0.0
        return nm.strip().strip("'\""), x

    def flush():
        if not buf.strip():
            return
        nm, ln = lab(buf)
        if closed is not None:
            closed.length = ln
            if nm and not re.match(r"^[\d.]+$", nm):
                closed.name = nm
        else:
            t = N(nm, ln); t.parent = cur; cur.children.append(t)
    for ch in s:
        if ch == "(":
            x = N(); x.parent = cur; cur.children.append(x); cur = x
            buf = ""; closed = None
        elif ch == ",":
            flush(); buf = ""; closed = None
        elif ch == ")":
            flush(); buf = ""; closed = cur; cur = cur.parent
        else:
            buf += ch
    flush()
    return root


def patr(a, b):
    pa, x, d = {}, a, 0.0
    while x is not None:
        pa[id(x)] = d; d += (x.length or 0.0); x = x.parent
    x, d = b, 0.0
    while x is not None:
        if id(x) in pa:
            return pa[id(x)] + d
        d += (x.length or 0.0); x = x.parent
    return float("nan")


vp = os.path.join(D, "mttB_verdicts.tsv")
if os.path.exists(vp):
    V = list(csv.DictReader(open(vp, errors="replace"), delimiter="\t"))
    for fn in ("mttB.contree", "mttB.treefile"):
        p = os.path.join(D, fn)
        if not os.path.exists(p):
            continue
        root = parse(open(p, errors="replace").read())
        byname = {t.name: t for t in root.tips()}
        diffs = []
        print("   --- %s ---" % fn)
        for r in V:
            qn = (r.get("query") or "").strip()
            ref = (r.get("closest_ref") or "").strip()
            try:
                want = float(r.get("dist_to_closest_ref") or "nan")
            except ValueError:
                continue
            qt = byname.get("QUERY|" + qn) or byname.get(qn)
            rt = byname.get(ref)
            if not (qt and rt):
                continue
            got = patr(qt, rt)
            diffs.append(want - got)
            print("      %-40s file %.4f  patristic %.4f  Δ %+.4f"
                  % (qn.split("|")[-1], want, got, want - got))
        if diffs:
            print("      Δ: min %+.4f  max %+.4f  spread %.5f"
                  % (min(diffs), max(diffs), max(diffs) - min(diffs)))
            if max(diffs) - min(diffs) < 1e-4:
                print("      -> a CONSTANT offset: the same metric on a tree")
                print("         whose branch lengths differ slightly.")
        print()
PYEOF

echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   A decides whether the mttB question is answerable at all. If the"
echo "   references that overlap our sequences include archaeal"
echo "   trimethylamine methyltransferases, the panel is rebuildable and the"
echo "   tree can be redone properly. If nothing overlaps, no reference set"
echo "   will help and the five proteins themselves need re-examining."
echo
echo "   B is the sentence that survives whatever A says: gene presence above"
echo "   a curated threshold, which needs no alignment."
echo
echo "   C closes the last open question about what the verdict tables report."
