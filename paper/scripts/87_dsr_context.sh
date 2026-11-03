#!/usr/bin/env bash
# 87_dsr_context.sh -- the facts panel (a) of fig13 needs, and the A/B split.
#
# TWO THINGS CAME OUT OF THE FIRST RENDER.
#
# 1. THE dsr TREE CONTAINS BOTH SUBUNITS. Consecutive accessions sit in
#    opposite halves -- AFQ44482/AFQ44483, AAL57464/AAL57463, BAD05998/BAD05997,
#    AAK83214/AAK83213, SDX68115/SDX68081 -- which are adjacent genes of one
#    operon, i.e. dsrA and dsrB. The deepest split in that tree is the ancient
#    A/B PARALOGY, not a relationship between organisms, and that is why the
#    reverse-type references are not monophyletic: there are two rDsr clades,
#    one per subunit. Rooting on them was the wrong operation. The figure must
#    show two labelled halves, and this step decides WHICH half is which from
#    the blastp table rather than from my reading of a picture.
#
# 2. THE FIGURE HAS NO CONTEXT. A reader meets a tip called
#    LIQ_contig_19299_14 with no way to know it is an unbinned contig from the
#    culture assembly carrying the only assembled dsrAB in the dataset. "Why
#    are these sequences examined on their own" is answered by the fact that
#    NO MAG CONTAINS THEM, and the figure never says so.
#
# Read-only. Writes two tables with a fixed schema that fig13 then consumes, so
# the figure never has to guess at a column.
set -u
WS2="${WS2}"
OUT="$WS2/87_dsr_context"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }
CONTIGS="LIQ_contig_12880 LIQ_contig_19299"

echo "=============================================================="
echo "A. which query is dsrA and which is dsrB"
echo "=============================================================="
# The subunit is decided by what each query's BEST reference is, and the
# reference sets are named in the blastp table. Nothing here is inferred from
# the tree -- the tree is what this will be used to label.
B="$WS2/58_dsr_finish/blastp.tsv"
if [ -s "$B" ]; then
  echo "   from ${B#$WS2/}, best hit per query:"
  "$PY" - "$B" "$OUT" <<'PYEOF'
import sys, os, collections
src, OUT = sys.argv[1], sys.argv[2]
best = {}
for ln in open(src, errors="replace"):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 8:
        continue
    q, r = f[0], f[1]
    try:
        bits = float(f[-1])
    except ValueError:
        continue
    if q not in best or bits > best[q][1]:
        best[q] = (r, bits, f[2])
rows = []
for q in sorted(best):
    r, bits, pid = best[q]
    acc = r.split("|")[-1]
    org = r.split("|")[1] if r.count("|") >= 2 else r
    print("     %-34s %-46s %6s%%  bits %s"
          % (q.split("|")[-1], org.replace("_", " ")[:46], pid, bits))
    rows.append((q.split("|")[-1], org, acc, pid, bits))
with open(os.path.join(OUT, "query_best.tsv"), "w") as fh:
    fh.write("query\torganism\taccession\tpident\tbits\n")
    for r in rows:
        fh.write("%s\t%s\t%s\t%s\t%s\n" % r)
print()
print("   ACCESSION PAIRING -- the evidence that this is a two-subunit tree.")
print("   Adjacent accessions in one genome are adjacent genes in one operon;")
print("   for dsr that is dsrA and dsrB. Pairs found among the references:")
PYEOF
else
  echo "   !! $B missing"
fi
echo

echo "=============================================================="
echo "B. the two contigs: size, depth, and what else is on them"
echo "=============================================================="
for c in $CONTIGS; do
  echo "--- $c"
  # length, from whichever assembly fasta holds it
  for fa in "$WS2"/03_assembly/*/assembly.fasta "$WS2"/03_assembly/*.fasta \
            "$WS2"/49b_contigs/*.fasta "$WS2"/55_dsr_id/*.fasta \
            "$WS2"/57_dsr_bin/dsr_contigs.fasta; do
    [ -s "$fa" ] || continue
    L=$("$PY" - "$fa" "$c" <<'PYEOF'
import sys
fa, want = sys.argv[1], sys.argv[2]
n, on = 0, False
for ln in open(fa, errors="replace"):
    if ln.startswith(">"):
        if on:
            break
        on = ln[1:].split()[0] == want
    elif on:
        n += len(ln.strip())
print(n if n else "")
PYEOF
)
    [ -n "$L" ] && { echo "     length $L bp   (from ${fa#$WS2/})"; break; }
  done
  # coverage / reads
  grep -h "$c" "$WS2"/58_dsr_finish/coverage.tsv 2>/dev/null | head -3 \
    | awk -F'\t' '{printf "     coverage row: reads %s, covbases %s, meandepth %s\n",$4,$5,$7}'
  # is it in any bin?  the answer the figure exists to state
  hit=$(grep -l "^>$c\b" "$WS2"/07_mags/bins/*.fa 2>/dev/null | head -3)
  if [ -n "$hit" ]; then
    echo "     IN A BIN: $(echo "$hit" | sed 's#.*/##' | tr '\n' ' ')"
  else
    echo "     in no MAG  <-- this is why it is examined on its own"
  fi
done
echo

echo "=============================================================="
echo "C. gene coordinates on those contigs"
echo "=============================================================="
# panel (a) draws the genes to scale, so it needs start/end/strand. Prodigal
# writes them into the protein fasta header; several steps may hold one.
for faa in "$WS2"/55_dsr_id/*.faa "$WS2"/58_dsr_finish/*.faa \
           "$WS2"/59_dsr_close/*.faa "$WS2"/49b_contigs/*.faa; do
  [ -s "$faa" ] || continue
  n=$(grep -c "^>\(LIQ_contig_12880\|LIQ_contig_19299\)" "$faa" 2>/dev/null || echo 0)
  [ "$n" -eq 0 ] && continue
  echo "--- ${faa#$WS2/}  ($n matching header(s))"
  grep "^>\(LIQ_contig_12880\|LIQ_contig_19299\)" "$faa" | cut -c1-150 | sed 's/^/     /'
done
echo
echo "   If no header above carries '# start # end # strand', prodigal"
echo "   coordinates are not on disk for these contigs and panel (a) must"
echo "   draw the genes as ordered blocks WITHOUT a bp scale, saying so,"
echo "   rather than inventing positions."
echo

echo "=============================================================="
echo "D. the catalogue-level fact panel (a) has to state"
echo "=============================================================="
KO="$WS2/45_kofam_full/all_hits.tbl"
THR="$WS2/36_gapclose/ko_list"
if [ -s "$KO" ]; then
  for k in K11180 K11181; do
    n=$(awk -v k="$k" '$1!~/^#/ && $3==k {print $1}' "$KO" 2>/dev/null | sort -u | wc -l)
    printf '   %-8s proteins with any hit across all genomes: %s\n' "$k" "$n"
  done
fi
echo "   (the report's claim is 0 of 29 MAGs carry dsrAB; any non-zero count"
echo "    here is hits in the UNBINNED set and must not be read as MAG genes)"
echo
echo "   Desulfobacterota share of the sediment community, for the same panel:"
grep -m2 -i "Desulfobacterota" "$WS2"/22_singlem/SED_profile.tsv 2>/dev/null \
  | cut -c1-140 | sed 's/^/     /'
