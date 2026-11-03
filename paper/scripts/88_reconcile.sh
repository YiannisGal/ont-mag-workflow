#!/usr/bin/env bash
# 88_reconcile.sh -- a number in the report does not match the file, and two
# greps of mine were wrong.
#
# 1. THE DISCREPANCY. FINDINGS_REPORT sections 7 and 12 state
#      "SingleM: Desulfobacterota 7.01% of the sediment community, third phylum
#       behind Chloroflexota (21.61%) and Bacillota (13.65%)"
#    Step 87 read the same file and got
#      sediment_reads  23.41  Root; d__Bacteria; p__Desulfobacterota
#    23.41% would make it the MOST abundant phylum in the sediment, not the
#    third. One of those is wrong and the difference is not cosmetic: section 7
#    uses it to argue that sulfate reducers are abundant in situ, which is the
#    competition premise of the chapter.
#
#    Possible causes, all testable here: two different files; two different
#    columns in one file (SingleM writes coverage AND relative abundance);
#    rows at different taxonomic ranks being compared; or a transcription error
#    into the report. This step prints the whole phylum block of every profile
#    on disk, with the column headers, and lets the numbers speak.
#
# 2. MY TWO GREPS IN STEP 87 SECTION C.
#      n=$(grep -c ... || echo 0)   prints 0 TWICE when grep matches nothing
#                                   and exits 1 -- hence "[: 0\n0: integer
#                                   expression expected"
#      grep "^>LIQ_contig_..."      anchors at the start of the header, but
#                                   these headers are ">QUERY|LIQ_contig_12880_2"
#                                   so the anchor guaranteed zero matches
#    The coordinates may well be on disk. This looks properly.
set -u
WS2="${WS2}"
OUT="$WS2/88_reconcile"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

echo "=============================================================="
echo "A. every SingleM / taxonomy profile of the SEDIMENT, in full"
echo "=============================================================="
for f in "$WS2"/22_singlem/SED_profile.tsv "$WS2"/22_singlem/*SED*.tsv \
         "$WS2"/24_deep_taxonomy/singlem_ranks.tsv \
         "$WS2"/24_deep_taxonomy/*SED*.tsv; do
  [ -s "$f" ] || continue
  echo "--- ${f#$WS2/}   ($(wc -l < "$f") lines)"
  echo "    header:"
  head -1 "$f" | sed 's/^/      /'
  echo "    PHYLUM-level rows (exactly two taxonomy fields: domain; phylum):"
  "$PY" - "$f" <<'PYEOF'
import sys
p = sys.argv[1]
rows = []
for i, ln in enumerate(open(p, errors="replace")):
    f = ln.rstrip("\n").split("\t")
    tax = f[-1]
    if "p__" not in tax or "c__" in tax:
        continue
    # the numeric columns, whatever they are called
    nums = []
    for x in f[:-1]:
        try:
            nums.append(float(x))
        except ValueError:
            pass
    rows.append((nums, tax.split(";")[-1].strip()))
rows.sort(key=lambda r: -(r[0][0] if r[0] else 0))
tot = sum(r[0][0] for r in rows if r[0])
for nums, name in rows[:14]:
    print("      %-34s %s" % (name, "  ".join("%9.3f" % n for n in nums)))
print("      %-34s %9.3f  <- sum of the first numeric column" % ("TOTAL", tot))
PYEOF
  echo
done

echo "=============================================================="
echo "B. what does the report's 7.01 correspond to, if anything?"
echo "=============================================================="
# search every file on disk for the literal values the report quotes, so the
# source of each one is identified rather than guessed
for v in 7.01 23.41 21.61 13.65; do
  echo "--- value $v"
  grep -rl --include="*.tsv" --include="*.txt" --include="*.json" \
       -F "$v" "$WS2"/22_singlem "$WS2"/24_deep_taxonomy "$WS2"/34_unbinned \
       "$WS2"/50b_harvest* "$WS2"/REPORT_INPUTS*.txt 2>/dev/null | head -4 \
    | while read -r h; do
        echo "     ${h#$WS2/}"
        grep -m2 -F "$v" "$h" | cut -c1-130 | sed 's/^/        /'
      done
done
echo

echo "=============================================================="
echo "C. gene coordinates -- searched properly this time"
echo "=============================================================="
# unanchored, and the count does not use "|| echo 0"
for faa in "$WS2"/55_dsr_id/*.faa "$WS2"/57_dsr_bin/*.faa \
           "$WS2"/58_dsr_finish/*.faa "$WS2"/59_dsr_close/*.faa \
           "$WS2"/49b_contigs/*.faa "$WS2"/34_unbinned/*.faa; do
  [ -s "$faa" ] || continue
  n=$(grep -c "LIQ_contig_12880\|LIQ_contig_19299" "$faa" 2>/dev/null) || n=0
  [ "${n:-0}" -eq 0 ] && continue
  echo "--- ${faa#$WS2/}  ($n header(s))"
  grep "LIQ_contig_12880\|LIQ_contig_19299" "$faa" | grep "^>" \
    | cut -c1-165 | sed 's/^/     /'
done
echo
echo "   prodigal writes '# start # end # strand # ID=...' into the header."
echo "   If that is present above, panel (a) draws the genes to scale."
echo "   If not, it draws them as ordered blocks and says the positions are"
echo "   not available -- it does not invent them."
echo

echo "=============================================================="
echo "D. the dsrA / dsrB assignment, stated"
echo "=============================================================="
echo "   from step 87's best-hit table:"
echo "     LIQ_contig_12880_2   D. meridiei  89.5%  bits 676   full-length"
echo "     LIQ_contig_12880_3   D. meridiei  88.5%  bits 659   full-length"
echo "     LIQ_contig_19299_13  D. putei     78.0%  bits 605   full-length"
echo "     LIQ_contig_19299_14  D. alkaliphila 82.9%  bits 373  HALF a gene"
echo "     LIQ_contig_19299_15  D. aeronauticus 75.5%  bits 234  HALF a gene"
echo
echo "   Which of each pair is A and which is B is decided by the SUBUNIT of"
echo "   the best reference, not by position in a tree. Print it:"
for acc in AFQ44482 AFQ44483 SHE97896 AAL57464 AAL57463 SHK50728; do
  hit=$(grep -m1 -h "$acc" "$WS2"/58_dsr_finish/*.faa 2>/dev/null | cut -c1-150)
  [ -n "$hit" ] && echo "     $acc  $hit" || echo "     $acc  (no header on disk)"
done
