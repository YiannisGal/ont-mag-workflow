#!/usr/bin/env bash
# 94_mldist_reader_fix.sh -- step 93's mttB numbers are wrong, and the file on
# disk says so.
#
# STEP 93 REPORTED, for every one of the five mttB queries against BOTH panels,
# a nearest reference at 0,0000. The project's own verdict table, written from
# the same distance matrix, says
#
#     ETO_LIQ_MAG_003|...|ETOLIQMAG003_01500   1.2509   TMA_MttB|CAQ5842095.1
#     ETO_LIQ_MAG_010|...|ETOLIQMAG010_13150   1.2830   TMA_MttB|CAQ5842095.1
#
# 0,0000 and 1,2509 cannot both be the distance between the same two
# sequences. One of them is a parsing error, and since the verdict table was
# written by the pipeline and my reader was written yesterday, mine is the one
# under suspicion. Two further symptoms point the same way:
#
#   * my reader found 78 rows in mttB.mldist for a tree with 85 tips;
#   * 15_marker_trees/mttB/ contains mttB.uniqueseq.phy, which IQ-TREE writes
#     only when it has removed IDENTICAL sequences -- and 85 - 78 = 7.
#
# So there are two candidate explanations and they are not the same thing.
# Either IQ-TREE genuinely dropped seven duplicate sequences and the matrix is
# an honest 78 x 78 -- in which case my distances should still have been right
# and something else is broken -- or the rows of the file WRAP across lines and
# my line-per-row reader silently mis-assembled the matrix.
#
# NOTHING FROM STEP 93 SECTIONS A-D FOR mttB MAY BE QUOTED UNTIL THIS RUNS.
# Its section E, which reads the alignments directly and never touches the
# distance matrix, is unaffected and stands.
#
# THE READER IS REPLACED WITH ONE THAT VALIDATES. It reads the count from the
# header, then consumes exactly that many numbers per row ACROSS LINE BREAKS,
# and then checks three things that any distance matrix must satisfy: it is
# square, its diagonal is zero, and it is symmetric. A matrix failing any of
# them is not used.
set -u
WS2="${WS2}"
OUT="$WS2/94_mldist_fix"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. what the files physically look like"
echo "=============================================================="
for m in mcrA mttB; do
    f="$WS2/15_marker_trees/$m/$m.mldist"
    [ -f "$f" ] || { echo "   MISSING $f"; continue; }
    echo "   --- $m.mldist: $(wc -l < "$f") line(s), $(wc -c < "$f") bytes ---"
    echo "   first line:"
    head -1 "$f" | sed 's/^/      /'
    echo "   second line, first 300 characters:"
    sed -n '2p' "$f" | cut -c1-300 | sed 's/^/      /'
    echo "   second line is $(sed -n '2p' "$f" | wc -w) whitespace-separated token(s)"
    echo "   third line, first 120 characters:"
    sed -n '3p' "$f" | cut -c1-120 | sed 's/^/      /'
    echo
done
echo "   and whether IQ-TREE removed identical sequences:"
for m in mcrA mttB; do
    u="$WS2/15_marker_trees/$m/$m.uniqueseq.phy"
    if [ -f "$u" ]; then
        echo "      $m: uniqueseq.phy EXISTS -- $(head -1 "$u")"
    else
        echo "      $m: no uniqueseq.phy (no identical sequences removed)"
    fi
    grep -iE "identical|duplicat|ignored|kept for" "$WS2/15_marker_trees/$m/$m.log" 2>/dev/null \
        | head -8 | sed 's/^/         /'
done

echo
echo "=============================================================="
echo "B. re-read with a validating reader, and redo the measurement"
echo "=============================================================="
"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re, collections, statistics
WS2, OUT = sys.argv[1], sys.argv[2]
M = os.path.join(WS2, "15_marker_trees")


def read_mldist(path):
    """PHYLIP-style lower/full matrix: a count, then per row a name followed by
    exactly `n` numbers, which MAY be wrapped over several lines. Reading one
    row per LINE is what a line-based reader gets wrong."""
    toks = open(path, errors="replace").read().split()
    if not toks:
        raise ValueError("empty file")
    n = int(toks[0])
    i = 1
    names, rows = [], []
    for _ in range(n):
        if i >= len(toks):
            raise ValueError("ran out of tokens after %d row(s)" % len(names))
        names.append(toks[i]); i += 1
        vals = []
        while len(vals) < n:
            if i >= len(toks):
                raise ValueError("row %s has only %d of %d value(s)"
                                 % (names[-1], len(vals), n))
            vals.append(float(toks[i])); i += 1
        rows.append(vals)
    leftover = len(toks) - i
    return n, names, rows, leftover


def validate(n, names, rows):
    bad = []
    if len(names) != n or len(rows) != n:
        bad.append("not %d rows" % n)
    for k, r in enumerate(rows):
        if len(r) != n:
            bad.append("row %d has %d values" % (k, len(r)))
    for k in range(min(n, len(rows))):
        if abs(rows[k][k]) > 1e-9:
            bad.append("diagonal %d is %.4f, not zero" % (k, rows[k][k]))
            break
    for a in range(min(n, len(rows))):
        for b in range(a + 1, min(n, len(rows))):
            if abs(rows[a][b] - rows[b][a]) > 1e-6:
                bad.append("not symmetric at (%d,%d): %.4f vs %.4f"
                           % (a, b, rows[a][b], rows[b][a]))
                break
        if bad and "symmetric" in bad[-1]:
            break
    return bad


ARCH = ("methano", "archaeoglob", "halo", "thermoplasma", "nitrososphaer",
        "sulfolobus", "thermococc", "pyrococc", "bathyarchae", "candidatus methano")


def is_arch(o):
    o = (o or "").lower()
    return any(a in o for a in ARCH)


def deflines(d):
    out = {}
    rd = os.path.join(d, "refs")
    if not os.path.isdir(rd):
        return out
    for fn in sorted(os.listdir(rd)):
        if not fn.endswith((".faa", ".fa", ".fasta")):
            continue
        for ln in open(os.path.join(rd, fn), errors="replace"):
            if ln.startswith(">"):
                h = ln[1:].rstrip("\n")
                a = h.split()[0]
                out[a] = h[len(a):].strip()
                out[a.split("|")[-1]] = out[a]
    return out


def org_of(d):
    m = re.search(r"\[([^\]]+)\]", d or "")
    return m.group(1) if m else ""


for marker in ("mcrA", "mttB"):
    d = os.path.join(M, marker)
    p = os.path.join(d, "%s.mldist" % marker)
    print("=" * 62)
    print("%s" % marker)
    print("=" * 62)
    if not os.path.exists(p):
        print("   MISSING %s" % p); continue
    try:
        n, names, rows, leftover = read_mldist(p)
    except ValueError as e:
        print("   !! cannot read: %s" % e); continue
    print("   header says %d sequence(s); read %d row(s); %d token(s) left over"
          % (n, len(names), leftover))
    bad = validate(n, names, rows)
    if bad:
        print("   !! the matrix does not validate: %s" % "; ".join(bad[:4]))
        print("      Not used. A distance matrix that is not square, zero on")
        print("      the diagonal and symmetric is not a distance matrix.")
        continue
    print("   validates: square, zero diagonal, symmetric")

    DESC = deflines(d)
    idx = {nm: i for i, nm in enumerate(names)}

    def panel(t):
        return t.split("|")[0] if "|" in t else "?"

    def accn(t):
        return t.split("|")[-1]

    pan = collections.defaultdict(list)
    for t in names:
        pan[panel(t)].append(t)
    for k in sorted(pan):
        print("   panel %-12s %3d" % (k, len(pan[k])))

    print()
    print("   reference-to-nearest-same-panel-reference, the scale of this"
          " alignment:")
    for k in sorted(pan):
        if k == "QUERY" or len(pan[k]) < 2:
            continue
        nn = [min(rows[idx[t]][idx[u]] for u in pan[k] if u != t) for t in pan[k]]
        z = sum(1 for x in nn if x == 0.0)
        print("      %-12s min %.4f  median %.4f  max %.4f   (%d at exactly 0)"
              % (k, min(nn), statistics.median(nn), max(nn), z))

    qs = [t for t in names if panel(t) == "QUERY"]
    if not qs:
        print("   no QUERY rows in this matrix"); print(); continue
    print()
    print("   %-26s %-12s %9s  %s" % ("query", "panel", "nearest", "which"))
    for q in qs:
        i = idx[q]
        f = q.split("|")
        lab = (f[1].replace("ETO_", "") + " " + f[3]) if len(f) > 3 else q
        for k in sorted(pan):
            if k == "QUERY":
                continue
            dist, who = min((rows[i][idx[u]], u) for u in pan[k])
            o = org_of(DESC.get(accn(who), "")) or accn(who)
            print("   %-26s %-12s %9.4f  %s" % (lab[:26], k, dist, o[:36]))
        arcs = [u for u in names if panel(u) != "QUERY"
                and is_arch(org_of(DESC.get(accn(u), "")))]
        if arcs:
            dist, who = min((rows[i][idx[u]], u) for u in arcs)
            print("   %-26s %-12s %9.4f  %s"
                  % ("", "ARCHAEAL", dist,
                     org_of(DESC.get(accn(who), ""))[:36]))
        print()

    # the comparison the verdict table makes, reproduced
    vp = os.path.join(d, "%s_verdicts.tsv" % marker)
    if os.path.exists(vp):
        print("   against %s_verdicts.tsv, which the pipeline wrote from this"
              " same matrix:" % marker)
        import csv
        for r in csv.DictReader(open(vp, errors="replace"), delimiter="\t"):
            qn = (r.get("query") or "").strip()
            cand = [t for t in names if t.endswith(qn) or qn.endswith(accn(t))
                    or accn(t) == accn(qn)]
            ref = (r.get("closest_ref") or "").strip()
            try:
                want = float(r.get("dist_to_closest_ref") or "nan")
            except ValueError:
                want = float("nan")
            if cand and ref in idx:
                got = rows[idx[cand[0]]][idx[ref]]
                ok = abs(got - want) < 5e-4
                print("      %-46s file %.4f  reread %.4f  %s"
                      % (accn(qn), want, got, "AGREE" if ok else "DISAGREE"))
            else:
                print("      %-46s cannot locate %s in the matrix"
                      % (accn(qn), ref if ref not in idx else qn))
    print()
PYEOF

echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   If the re-read now agrees with the verdict tables, step 93's mttB"
echo "   numbers were a reader bug and the real question -- whether our"
echo "   archaeal MttB sequences have any close relative in either panel --"
echo "   is answered by section B's table."
echo
echo "   If the matrix does not validate, no number is taken from it at all"
echo "   and the mttB placement has to be redone from the alignment."
echo
echo "   Step 93 section E is untouched either way: it reads the alignment"
echo "   files directly and never opens a distance matrix."
