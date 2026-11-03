#!/usr/bin/env bash
# 95_mttb_overlap.sh -- two corrections of mine, and the measurement that
# decides whether the mttB tree can carry a claim at all.
#
# CORRECTION 1. I BLAMED MY OWN READER AND THE READER WAS RIGHT.
#
# Step 94's validating reader reports mttB.mldist as 78 x 78, square, zero on
# the diagonal, symmetric, with no tokens left over -- and it returns the SAME
# near-zero distances step 93 did. The header says 78 and there are 78 rows;
# the count matched the tree's 85 tips minus the 7 identical sequences IQ-TREE
# folded together. Nothing was mis-parsed. The file genuinely contains
#
#     QUERY|...MAG018_18010   0.0000000 0.0603920 ... 0.0000010 1.4837095
#                                                     ^^^^^^^^^
# a distance of 1e-6 between one of our archaeal sequences and a bacterial
# reference. 1e-6 is IQ-TREE'S FLOOR: what it writes when the maximum
# likelihood distance cannot be estimated. It is not a measurement of
# similarity, it is the absence of one.
#
# CORRECTION 2. THE VERDICT TABLES ARE NOT READ FROM THE DISTANCE MATRIX.
#
# Step 94 flagged 7 of 9 mcrA rows as disagreeing with mcrA_verdicts.tsv --
# 0.0760 against 0.0660, 0.0972 against 0.0912, 0.1329 against 0.1308. Those
# are not rounding. But every one of them matches the PATRISTIC distance
# fig15 measures along the tree: 0,076, 0,097, 0,133. So
# dist_to_closest_ref is a patristic distance, and step 94 compared two
# different quantities and called the difference a disagreement. Section D
# proves that by recomputing patristic distances and expecting them to agree.
#
# WHAT ACTUALLY HAS TO BE MEASURED.
#
# A maximum likelihood distance between two aligned sequences is estimated from
# the columns where BOTH have a residue. Step 93 found that no column of the
# mttB alignment is occupied by 90% of the 40 TMA_MttB references, and none by
# 90% of all sequences. If our queries and their "nearest" references share
# almost no columns, the distance between them is undefined, IQ-TREE emits its
# floor, and every mttB number that rests on that matrix is meaningless --
# including the k-nearest-neighbour composition behind the CONFIRMED verdicts.
#
# Section A counts, for every query-reference pair, the number of alignment
# columns where both sequences have a residue. That is the quantity the
# distance is estimated from, and it is a fact about the alignment file.
# Section B does the same for mcrA, where it should come out large.
set -u
WS2="${WS2}"
OUT="$WS2/95_mttb_overlap"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re, collections, statistics
WS2, OUT = sys.argv[1], sys.argv[2]
M = os.path.join(WS2, "15_marker_trees")


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


def read_mldist(path):
    toks = open(path, errors="replace").read().split()
    n = int(toks[0]); i = 1
    names, rows = [], []
    for _ in range(n):
        names.append(toks[i]); i += 1
        vals = []
        while len(vals) < n:
            vals.append(float(toks[i])); i += 1
        rows.append(vals)
    return n, names, rows


def panel(t):
    return t.split("|")[0] if "|" in t else "?"


def accn(t):
    return t.split("|")[-1]


GAPS = set("-.")
FLOOR = 1e-6

for marker in ("mttB", "mcrA"):
    d = os.path.join(M, marker)
    print("=" * 64)
    print("%s" % marker)
    print("=" * 64)
    aln = None
    for c in ("%s_trim.faa" % marker, "%s_aln.faa" % marker):
        if os.path.exists(os.path.join(d, c)):
            aln = os.path.join(d, c); break
    if aln is None:
        print("   no alignment file"); continue
    A = read_fa(aln)
    L = len(next(iter(A.values())))
    print("   alignment %s: %d sequence(s) x %d columns"
          % (os.path.basename(aln), len(A), L))

    # residue masks, one bit per column
    MASK = {k: [c not in GAPS for c in v] for k, v in A.items()}
    NRES = {k: sum(m) for k, m in MASK.items()}
    qs = [k for k in A if panel(k) == "QUERY"]
    refs = [k for k in A if panel(k) != "QUERY"]
    print("   %d query, %d reference sequence(s)" % (len(qs), len(refs)))
    print("   residues per sequence: queries median %d, references median %d"
          % (statistics.median([NRES[k] for k in qs]) if qs else 0,
             statistics.median([NRES[k] for k in refs]) if refs else 0))

    md = os.path.join(d, "%s.mldist" % marker)
    DIST, NAMES = None, None
    if os.path.exists(md):
        n, NAMES, rows = read_mldist(md)
        idx = {x: i for i, x in enumerate(NAMES)}
        DIST = (idx, rows)
        flat = [rows[i][j] for i in range(n) for j in range(i + 1, n)]
        nfloor = sum(1 for x in flat if x <= FLOOR * 1.5)
        print("   matrix %d x %d: %d of %d pairwise value(s) are at IQ-TREE's"
              " floor (<= 1,5e-6) = %.1f%%"
              % (n, n, nfloor, len(flat), 100.0 * nfloor / len(flat)))
        big = sum(1 for x in flat if x > 2.0)
        print("   %d of %d (%.1f%%) exceed 2,0 substitutions per site"
              % (big, len(flat), 100.0 * big / len(flat)))

    print()
    print("   SHARED COLUMNS between each query and its nearest reference in")
    print("   each panel -- the positions the distance is estimated from:")
    print("   %-24s %-12s %8s %9s %9s" % ("query", "panel", "shared", "distance",
                                          "verdict"))
    bypanel = collections.defaultdict(list)
    for r in refs:
        bypanel[panel(r)].append(r)
    rows_out = []
    for q in sorted(qs):
        for pn in sorted(bypanel):
            cand = bypanel[pn]
            if DIST:
                idx, mat = DIST
                have = [r for r in cand if r in idx and q in idx]
                if not have:
                    continue
                best = min(have, key=lambda r: mat[idx[q]][idx[r]])
                dist = mat[idx[q]][idx[best]]
            else:
                best, dist = cand[0], float("nan")
            shared = sum(1 for a, b in zip(MASK[q], MASK[best]) if a and b)
            verdict = ("UNDEFINED" if dist <= FLOOR * 1.5 else
                       ("thin" if shared < 50 else "ok"))
            f = q.split("|")
            print("   %-24s %-12s %8d %9.6f %9s"
                  % ((f[1].replace("ETO_", "") if len(f) > 1 else q)[:24],
                     pn, shared, dist, verdict))
            rows_out.append((marker, q, pn, best, shared, dist))
        # and the overlap with EVERY reference, summarised
        ov = [sum(1 for a, b in zip(MASK[q], MASK[r]) if a and b) for r in refs]
        print("   %-24s %-12s shared columns with ALL references:"
              " min %d, median %d, max %d"
              % ("", "", min(ov), int(statistics.median(ov)), max(ov)))
        print()

    with open(os.path.join(OUT, "%s_overlap.tsv" % marker), "w") as fh:
        fh.write("marker\tquery\tpanel\tnearest\tshared_columns\tmldist\n")
        for r in rows_out:
            fh.write("%s\t%s\t%s\t%s\t%d\t%.7f\n" % r)
    print("   -> 95_mttb_overlap/%s_overlap.tsv" % marker)
    print()

# ---------------------------------------------------------------- section D
print("=" * 64)
print("D. is dist_to_closest_ref a PATRISTIC distance?")
print("=" * 64)
print("   Recomputed along the tree. If these agree, step 94's mcrA")
print("   'DISAGREE' rows were my comparison being wrong, not the files.")
print()


class Node(object):
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name="", length=0.0):
        self.name = name; self.length = length
        self.children = []; self.parent = None

    def tips(self):
        if not self.children:
            return [self]
        out = []
        for c in self.children:
            out.extend(c.tips())
        return out


def parse(text):
    s = text.strip()
    if ";" in s:
        s = s[:s.rindex(";")]
    root = Node(); cur = root; buf = ""; closed = None

    def lab(tok):
        if ":" in tok:
            nm, _, ln = tok.rpartition(":")
        else:
            nm, ln = tok, ""
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
            t = Node(nm, ln); t.parent = cur; cur.children.append(t)

    for ch in s:
        if ch == "(":
            n = Node(); n.parent = cur; cur.children.append(n); cur = n
            buf = ""; closed = None
        elif ch == ",":
            flush(); buf = ""; closed = None
        elif ch == ")":
            flush(); buf = ""; closed = cur; cur = cur.parent
        else:
            buf += ch
    flush()
    return root


def patristic(a, b):
    pa, n, dd = {}, a, 0.0
    while n is not None:
        pa[id(n)] = dd; dd += (n.length or 0.0); n = n.parent
    n, dd = b, 0.0
    while n is not None:
        if id(n) in pa:
            return pa[id(n)] + dd
        dd += (n.length or 0.0); n = n.parent
    return float("nan")


import csv
for marker in ("mcrA", "mttB"):
    d = os.path.join(M, marker)
    vp = os.path.join(d, "%s_verdicts.tsv" % marker)
    tp = None
    for c in ("%s.contree" % marker, "%s.treefile" % marker):
        if os.path.exists(os.path.join(d, c)):
            tp = os.path.join(d, c); break
    if not (os.path.exists(vp) and tp):
        print("   %s: missing verdicts or tree" % marker); continue
    root = parse(open(tp, errors="replace").read())
    byname = {t.name: t for t in root.tips()}
    print("   --- %s (%s) ---" % (marker, os.path.basename(tp)))
    agree = 0; total = 0
    for r in csv.DictReader(open(vp, errors="replace"), delimiter="\t"):
        qn = (r.get("query") or "").strip()
        ref = (r.get("closest_ref") or "").strip()
        try:
            want = float(r.get("dist_to_closest_ref") or "nan")
        except ValueError:
            continue
        qt = byname.get("QUERY|" + qn) or byname.get(qn)
        rt = byname.get(ref)
        if not (qt and rt):
            print("      %-42s tip not found in the tree" % accn(qn)); continue
        got = patristic(qt, rt)
        total += 1
        ok = abs(got - want) < 1e-3
        agree += 1 if ok else 0
        print("      %-42s file %.4f  patristic %.4f  %s"
              % (accn(qn), want, got, "AGREE" if ok else "DISAGREE"))
    print("      %d of %d agree" % (agree, total))
    print()
PYEOF

echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   If our mttB queries share only a handful of alignment columns with"
echo "   the references they are 'closest' to, then those distances are not"
echo "   small, they are UNDEFINED, and no statement about mttB placement can"
echo "   rest on that matrix or on the neighbourhood counts derived from it."
echo "   FINDINGS_REPORT section 5's MttB sentence then has to fall back to"
echo "   gene presence above the KOfam threshold, with no phylogeny."
echo
echo "   mcrA is the positive control in the same run: its queries should"
echo "   share most of the alignment with their nearest references."
echo
echo "   D settles what dist_to_closest_ref actually is, so that no future"
echo "   step compares it against a distance-matrix value again."
