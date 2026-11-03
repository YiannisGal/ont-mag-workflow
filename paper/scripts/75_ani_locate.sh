#!/usr/bin/env bash
# 75_ani_locate.sh -- find the ANI-against-cultured-isolates table.
#
# READ-ONLY except for its own output directory. Nothing is recomputed here;
# this only locates the file that already holds the result the chapter claims
# (0 of 29 MAGs above 95 % ANI to any cultured isolate, control 99,75 %) and
# prints enough of its shape that the figure script can be written against the
# real columns instead of a guess.
set -u
WS2="${WS2}"
OUT="$WS2/75_ani_locate"
mkdir -p "$OUT" || exit 1

echo "=============================================================="
echo "A. files whose NAME looks like an ANI result"
echo "=============================================================="
find "$WS2" -maxdepth 4 -type f \
     \( -iname "*ani*" -o -iname "*fastani*" -o -iname "*skani*" \
        -o -iname "*Ndb*" -o -iname "*Cdb*" -o -iname "*mash*" \
        -o -iname "*dist*" \) \
     ! -path "*/tmp/*" ! -name "*.log" ! -name "*.sh" 2>/dev/null \
  | while read -r f; do
      printf '%10s  %s\n' "$(du -h "$f" | cut -f1)" "${f#$WS2/}"
    done | sort -k2
echo

echo "=============================================================="
echo "B. tabular files whose HEADER carries an ANI-like column"
echo "=============================================================="
# name alone is unreliable -- a table can hold ANI under any filename, and a
# file called *ani* can hold something else. The header is the evidence.
find "$WS2" -maxdepth 4 -type f \( -name "*.tsv" -o -name "*.csv" -o -name "*.txt" \) \
     -size -200M ! -path "*/tmp/*" 2>/dev/null \
  | while read -r f; do
      h=$(head -c 4000 "$f" 2>/dev/null | head -1)
      case "$h" in
        *ANI*|*ani*|*identity*|*Identity*|*ident*)
          printf '%s\n    header: %s\n' "${f#$WS2/}" "$(printf '%s' "$h" | cut -c1-200)" ;;
      esac
    done
echo

echo "=============================================================="
echo "C. headerless ANI outputs (fastANI writes none)"
echo "=============================================================="
# fastANI emits: query <TAB> reference <TAB> ANI <TAB> fragments_mapped <TAB>
# fragments_total  -- no header line at all, so section B cannot see it. The
# signature is five tab-separated fields with a 70-100 number in field 3.
find "$WS2" -maxdepth 4 -type f \( -name "*.tsv" -o -name "*.txt" -o -name "*.out" \) \
     -size -200M ! -path "*/tmp/*" 2>/dev/null \
  | while read -r f; do
      awk -F'\t' 'NR==1 {
             if (NF==5 && $3+0 >= 70 && $3+0 <= 100 && $3 ~ /[0-9]\./) print "HIT";
             exit }' "$f" 2>/dev/null | grep -q HIT && {
        printf '%s\n    first line: %s\n' "${f#$WS2/}" "$(head -1 "$f" | cut -c1-200)"
        printf '    lines: %s\n' "$(wc -l < "$f")"
      }
    done
echo

echo "=============================================================="
echo "D. steps that MENTION an ANI tool, and where they wrote"
echo "=============================================================="
# If the table has been deleted or moved, the script that made it still says
# which tool ran and which directory it wrote into.
grep -l -i -E "fastANI|skani|dRep|ANIclustermap|average nucleotide|--ani" \
     "$WS2"/*.sh "$WS2"/scripts/*.sh 2>/dev/null \
  | while read -r s; do
      echo "--- ${s#$WS2/}"
      grep -n -i -E "fastANI|skani|dRep|ANI|OUT=|OUTDIR=" "$s" | head -14
    done
echo

echo "=============================================================="
echo "E. GTDB-Tk's own ANI screen, if it ran"
echo "=============================================================="
# gtdbtk classify writes ani_closest.tsv / *.ani_summary.tsv against its
# reference set. Those references are GENOMES IN GTDB, which is not the same
# question as "a cultured isolate" -- GTDB is full of MAGs. Worth finding, but
# it answers a different question and must not be used for this panel without
# filtering to isolate-derived references.
find "$WS2" -maxdepth 5 -type f \( -name "*ani_closest*" -o -name "*ani_summary*" \
     -o -name "*.classify.tree" -o -name "gtdbtk.*.summary.tsv" \) 2>/dev/null \
  | sed "s|$WS2/||"
echo

echo "=============================================================="
echo "F. the control the chapter quotes (99,75 %)"
echo "=============================================================="
# Whatever file holds the result must contain a value very close to 99.75 --
# that is the positive control, a genome matched against itself or against its
# own known relative. Finding 99.7x is how we confirm we found the right table.
find "$WS2" -maxdepth 4 -type f \( -name "*.tsv" -o -name "*.csv" -o -name "*.txt" \) \
     -size -200M ! -path "*/tmp/*" 2>/dev/null \
  | while read -r f; do
      if grep -q -E "99\.7[0-9]|99,7[0-9]" "$f" 2>/dev/null; then
        echo "${f#$WS2/}"
        grep -n -E "99\.7[0-9]|99,7[0-9]" "$f" | head -3 | sed 's/^/    /'
      fi
    done
echo

echo "=============================================================="
echo "SUMMARY"
echo "=============================================================="
echo "If sections A-C found a table: note its path, its header line, and"
echo "   head -5  of it. The figure reads that file directly."
echo "If they found nothing: the ANI comparison has not been run in workspace_2,"
echo "   and the chapter's number came from somewhere else; the step that computes"
echo "   it has to be written, rather than drawing a figure from a"
echo "   number with no output behind it."
