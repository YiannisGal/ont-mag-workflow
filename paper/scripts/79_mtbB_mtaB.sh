#!/usr/bin/env bash
# 79_mtbB_mtaB.sh -- the last two methylamine/methanol markers, measured the
# way mtmB was.
#
# NO NEW SEARCH IS NEEDED, AND THAT IS THE POINT. Step 61 searched with the
# whole MttB superfamily in one panel -- mtmB as target, mtbB, mtaB, mttB and
# mtmC as decoys able to WIN a read. A read was assigned to its best hit across
# all of them. So the competitive assignment for mtbB and mtaB has already been
# computed; it was simply never read off, because step 61 counted one label.
# Re-running the search would only reproduce the same table more slowly and
# would risk a different reference set, which is the exact failure that made
# the old counts unusable.
#
# WHAT MADE THE OLD NUMBERS UNREPORTABLE. mtbB at 4.9x and mtaB at 2.2x the
# sediment mcrA count violate a biological ceiling: these are single-copy-ish
# genes in organisms that all carry mcrA, so their reads cannot outnumber mcrA
# reads several-fold. Those counts came from panels that did not contain their
# own siblings, so sibling sequence was being scored as them.
#
# THE STEP CARRIES ITS OWN POSITIVE CONTROL. Before quoting anything new it
# recomputes mtmB from the same table and checks it reproduces the published
# 18 sediment reads and a ratio of 0.75 against a prediction of 0.78. If mtmB
# does not come back, the recount is wrong and no other number from it is
# trustworthy -- the step says so and stops.
set -u
WS2="${WS2}"
P61="$WS2/61_mtmB"
P63="$WS2/63_mtmB_relabel"
OUT="$WS2/79_mtbB_mtaB"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }
MINLEN=30          # step 63's own alignment-length cutoff, reused verbatim

echo "=============================================================="
echo "A. what step 61 left behind, and which column is which"
echo "=============================================================="
for f in "$P61"/*.tsv "$P61"/*.faa "$P61"/*.txt "$P63"/*.txt "$P63"/*.tsv; do
  [ -e "$f" ] || continue
  printf '   %9s  %s\n' "$(du -h "$f" | cut -f1)" "${f#$WS2/}"
done
echo
for tag in SED LIQ; do
  H="$P61/${tag}_hits.tsv"
  [ -s "$H" ] || { echo "   !! $H missing -- cannot proceed"; exit 1; }
  echo "   ${tag}_hits.tsv, first two rows:"
  head -2 "$H" | cat -A 2>/dev/null | sed 's/\$$//' | cut -c1-160 | sed 's/^/     /'
done
echo
echo "   column identification (rule 114: identify before comparing):"
"$PY" - "$P61/SED_hits.tsv" <<'PYEOF'
import sys
rows = []
for i, ln in enumerate(open(sys.argv[1], errors="replace")):
    f = ln.rstrip("\n").split("\t")
    rows.append(f)
    if i > 4000:
        break
n = min(len(r) for r in rows)
print("     %d columns" % n)
for c in range(n):
    vals = [r[c] for r in rows]
    ints = sum(1 for v in vals if v.lstrip("-").isdigit())
    flo = 0
    for v in vals:
        try:
            float(v); flo += 1
        except ValueError:
            pass
    uniq = len(set(vals))
    kind = "int" if ints == len(vals) else ("num" if flo == len(vals) else "text")
    print("     col %d  %-5s  distinct %6d  e.g. %s"
          % (c, kind, uniq, vals[0][:46]))
print()
print("     expected: col 0 read id, col 1 reference, col 3 alignment length,")
print("     col 5 bit score. If the shapes above disagree, STOP and fix the")
print("     indices below before believing any count.")
PYEOF
echo

echo "=============================================================="
echo "B. recount, one best hit per read, every label"
echo "=============================================================="
"$PY" - "$P61" "$P63" "$MINLEN" "$OUT" <<'PYEOF'
import sys, os, collections
P61, P63, MINLEN, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]

bact = set()
bp = os.path.join(P63, "bacterial.txt")
if os.path.exists(bp):
    bact = {l.strip() for l in open(bp) if l.strip()}
    print("   %d bacterial reference id(s) will be relabelled DECOYbact" % len(bact))
else:
    print("   !! %s missing." % bp)
    print("      Step 63 removed seven BACTERIAL proteins that an NCBI protein-name")
    print("      query had imported into the target set; they held 67 sediment reads.")
    print("      Without that list the counts below would repeat the error this")
    print("      whole thread exists to fix. Refusing to continue.")
    sys.exit(1)


def label(ref):
    """step 63's own rule: the panel tag is the first |-field, unless the
    reference is one of the bacterial imports."""
    p = ref.split("|")
    return "DECOYbact" if p[-1] in bact else p[0]


def best(tag):
    """one read, one hit -- the highest-scoring reference across the WHOLE
    panel. That competition is what makes the number mean anything."""
    path = os.path.join(P61, "%s_hits.tsv" % tag)
    top = {}
    kept = short = 0
    for ln in open(path, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        try:
            ln_aln, sc = int(f[3]), float(f[5])
        except ValueError:
            continue
        if ln_aln < MINLEN:
            short += 1
            continue
        kept += 1
        cur = top.get(f[0])
        if cur is None or sc > cur[1]:
            top[f[0]] = (f[1], sc)
    return top, kept, short


res = {}
for tag in ("SED", "LIQ"):
    top, kept, short = best(tag)
    c = collections.Counter(label(r) for r, _ in top.values())
    res[tag] = c
    print("   %s: %d alignments kept, %d below %d aa, %d reads assigned"
          % (tag, kept, short, MINLEN, len(top)))

labels = sorted(set(res["SED"]) | set(res["LIQ"]))
print()
print("   %-18s %8s %8s" % ("label", "SED", "LIQ"))
for l in labels:
    print("   %-18s %8d %8d" % (l, res["SED"].get(l, 0), res["LIQ"].get(l, 0)))

with open(os.path.join(OUT, "counts.tsv"), "w") as fh:
    fh.write("label\tSED\tLIQ\n")
    for l in labels:
        fh.write("%s\t%d\t%d\n" % (l, res["SED"].get(l, 0), res["LIQ"].get(l, 0)))
print()
print("   wrote %s" % os.path.join(OUT, "counts.tsv"))

# ---- the positive control, before anything new is quoted
mt = res["SED"].get("mtmB", 0)
mc = res["SED"].get("mcrA", 0)
print()
print("   CONTROL -- does this recount reproduce the published mtmB result?")
print("     mtmB sediment reads : %d   (published 18)" % mt)
print("     mcrA sediment reads : %d   (published 24)" % mc)
if mc:
    print("     mtmB / mcrA         : %.2f (published 0.75)" % (mt / float(mc)))
ok = (mt == 18 and mc == 24)
print("     %s" % ("[OK] reproduced -- the other labels can be read off the same table"
                   if ok else
                   "[--] NOT reproduced. The recount differs from the published run;\n"
                   "          every other number below is suspect. Do not quote them."))
# a failed control must be visible in the exit status, not only in the text
sys.exit(0 if ok else 2)
PYEOF
RC=$?
echo

echo "=============================================================="
echo "C. is each label's reference set able to carry a count?"
echo "=============================================================="
# rule 108: a fixed share threshold is a different test at every set size, so
# promiscuity is tested at 4x the EXPECTED share 1/n. And a label with very few
# references cannot be compared to one with many -- that has to be visible.
REFS=$(ls "$P61"/*.faa 2>/dev/null | head -1)
if [ -n "${REFS:-}" ]; then
  echo "   panel: ${REFS#$WS2/}"
  "$PY" - "$REFS" "$P61" "$MINLEN" "$P63" <<'PYEOF'
import sys, os, collections
refs, P61, MINLEN, P63 = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
bact = {l.strip() for l in open(os.path.join(P63, "bacterial.txt"))} \
    if os.path.exists(os.path.join(P63, "bacterial.txt")) else set()

per = collections.Counter()
alen = collections.defaultdict(list)
cur = None
seq = 0
for ln in open(refs, errors="replace"):
    if ln.startswith(">"):
        if cur:
            alen[cur].append(seq)
        i = ln[1:].split()[0]
        cur = "DECOYbact" if i.split("|")[-1] in bact else i.split("|")[0]
        per[cur] += 1
        seq = 0
    else:
        seq += len(ln.strip())
if cur:
    alen[cur].append(seq)

print("   %-18s %6s %10s" % ("label", "refs", "median aa"))
for l in sorted(per):
    m = sorted(alen[l])[len(alen[l]) // 2] if alen[l] else 0
    print("   %-18s %6d %10d" % (l, per[l], m))
print()
print("   A label with a handful of references cannot absorb reads the way one")
print("   with dozens can. Any ratio between two labels of very different set")
print("   size carries that bias and must be reported with both n values.")

# per-reference share within each label, flagged at 4x expected
for tag in ("SED", "LIQ"):
    top = {}
    for ln in open(os.path.join(P61, "%s_hits.tsv" % tag), errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        try:
            L, s = int(f[3]), float(f[5])
        except ValueError:
            continue
        if L < MINLEN:
            continue
        c = top.get(f[0])
        if c is None or s > c[1]:
            top[f[0]] = (f[1], s)
    byref = collections.Counter(r for r, _ in top.values())
    bylab = collections.defaultdict(collections.Counter)
    for r, n in byref.items():
        lab = "DECOYbact" if r.split("|")[-1] in bact else r.split("|")[0]
        bylab[lab][r] += n
    print()
    print("   %s -- references taking more than 4x their expected share:" % tag)
    hit = False
    for lab, cnt in sorted(bylab.items()):
        tot = sum(cnt.values())
        n = per.get(lab, 0)
        if not tot or not n:
            continue
        exp = 1.0 / n
        for r, c in cnt.most_common(3):
            share = c / float(tot)
            if share > 4 * exp and share > 0.20:
                hit = True
                print("     %-12s %-34s %4d/%4d = %.0f%%  (expected %.0f%%)"
                      % (lab, r.split("|")[-1][:34], c, tot, 100 * share, 100 * exp))
    if not hit:
        print("     none")
PYEOF
else
  echo "   !! no reference panel FASTA found in 61_mtmB/ -- section C skipped."
  echo "      Without it the per-label reference counts cannot be stated, and a"
  echo "      ratio between two labels of unknown set size is not interpretable."
fi
echo

echo "=============================================================="
echo "D. the ceiling, calibrated from this catalogue (rule 110)"
echo "=============================================================="
# The expectation is NOT 1:1. Methylotrophic methanogens carry paralogs of the
# methylamine methyltransferases while mcrA is single copy, so the ceiling has
# to come from our own gene counts. The same arithmetic reproduced mtmB at
# 0.78 against a measured 0.75; if it does that again here it can be trusted
# for the other two.
MASTER=$(ls "$WS2"/47_master/MASTER.tsv 2>/dev/null | head -1)
KO="$WS2/45_kofam_full/all_hits.tbl"
echo "   catalogue gene counts above curated thresholds are needed for:"
echo "     K00399 mcrA   K16176 mtmB   K16178 mtbB   K04480 mtaB"
if [ -s "$KO" ]; then
  echo "   from ${KO#$WS2/} (this is the FULL search; 11_kofam/hits holds only a"
  echo "   targeted profile set and absence there is not evidence of absence)"
  for k in K00399 K16176 K16178 K04480; do
    n=$(awk -v k="$k" '$1!~/^#/ && $3==k {print $1}' "$KO" 2>/dev/null | sort -u | wc -l)
    printf '     %-8s %4d protein(s) with any hit\n' "$k" "$n"
  done
  echo "   NOTE: those are ANY-hit counts. The reportable number is hits above"
  echo "   the KOfam curated threshold for each KO; step 45's own threshold"
  echo "   table is the place to take it from, and the mtmB value it must"
  echo "   reproduce is 11 proteins against 9 mcrA."
else
  echo "   !! $KO not found -- cannot calibrate. Do not quote a ratio."
fi
echo

echo "=============================================================="
echo "WHAT TO DO WITH THIS"
echo "=============================================================="
echo "If B's control reproduced mtmB = 18 and mcrA = 24:"
echo "   the mtbB and mtaB rows in B are measured the same way as the published"
echo "   mtmB number and become quotable, with their reference-set sizes from C"
echo "   stated alongside. Section 15 item 5 of the findings report can then be"
echo "   narrowed to whatever still fails D's ceiling."
echo "If it did not:"
echo "   nothing here is quotable. The difference between A and B has to be"
echo "   found before anything goes in the chapter."
exit $RC
