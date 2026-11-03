#!/usr/bin/env bash
# 60_readlevel_recon.sh -- how does the read-level marker pipeline actually work,
#                          and hand back the figure house style.
#
#   bash scripts/run.sh scripts/60_readlevel_recon_XXXXXXX.sh
#
# WHY DISCOVERY FIRST. The mtmB (K16176) read search has to produce a number
# that can sit in the SAME TABLE as the mcrA, mtbB, mtaB and dsrAB read counts.
# That means the identical reference-construction rule, the identical diamond
# or hmmsearch cutoffs, the identical gene-length normalisation and the
# identical per-million denominator. Guessing any one of those makes a number
# that looks comparable and is not -- which is exactly the failure that put
# mtbB at 4.9x and mtaB at 2.2x the sediment mcrA count and made both
# unreportable.
#
# There is also a live discrepancy to settle here. The sediment library is
# 7,528,854 reads as used, but step 34's per-million denominators were
# 6,827,767. Both numbers are in the record. This prints what the pipeline
# actually counted so the mtmB row uses the right one.
#
# AND IT DUMPS THE FIGURE HOUSE STYLE, so the remaining five figures can be
# written against the real API rather than an invented one.
#
# Runs in seconds. Reads nothing but scripts and existing outputs. Writes only
# 60_recon/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/60_recon"; mkdir -p "$OUT"
exec > >(tee "$OUT/60.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }

# dump a file with a cap, marking truncation honestly
dump(){   # $1 path  $2 head lines  $3 tail lines
  local f="$1" h="${2:-200}" t="${3:-0}" n
  if [ ! -e "$f" ]; then echo "    [MISSING] $f"; return; fi
  if [ ! -s "$f" ]; then echo "    [EMPTY]   $f"; return; fi
  n=$(wc -l < "$f")
  printf '    <<< %s  (%s lines, %s bytes, md5 %s)\n' \
    "${f#$WS2/}" "$n" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -c1-8)"
  if [ "$n" -le $((h+t)) ]; then sed 's/^/    /' "$f"
  else
    head -"$h" "$f" | sed 's/^/    /'
    printf '    ... [%d lines omitted] ...\n' $((n-h-t))
    [ "$t" -gt 0 ] && tail -"$t" "$f" | sed 's/^/    /'
  fi
  printf '    >>> end\n\n'
}

echo "60 -- read-level pipeline recon + figure style"
echo "date : $(date -Is)"
echo "out  : $OUT"

# ================================================================ A
say "A. THE READ-LEVEL SCRIPTS THEMSELVES"
for pat in '31*read*' '34*unbinned*' '*readlevel*' '*read_level*'; do
  find "$WS2/scripts" "$WS2" -maxdepth 2 -type f -name "$pat.sh" 2>/dev/null
done | sort -u | while read -r s; do
  [ -s "$s" ] && dump "$s" 400 0
done

# ================================================================ B
say "B. HOW THE MARKER REFERENCE SETS WERE BUILT"
for D in 31_readlevel 34_unbinned; do
  sub "$D/"
  [ -d "$WS2/$D" ] || { no "no such directory"; continue; }
  ls -la "$WS2/$D" 2>/dev/null | sed 's/^/    /'
  for f in "$WS2/$D"/*.faa "$WS2/$D"/*.fasta "$WS2/$D"/*.txt "$WS2/$D"/*.tsv; do
    [ -s "$f" ] || continue
    case "$f" in
      *.faa|*.fasta) printf '\n    %s : %s sequences\n' "${f#$WS2/}" "$(grep -c '^>' "$f")"
                     grep '^>' "$f" | head -12 | sed 's/^/      /'
                     printf '      (headers above are the naming convention to copy)\n' ;;
      *) dump "$f" 40 8 ;;
    esac
  done
done

# ================================================================ C
say "C. THE CUTOFFS AND DENOMINATORS, AS ACTUALLY USED"
sub "every diamond / hmmsearch invocation in the read-level scripts"
grep -rn -E 'diamond|hmmsearch|--id |--query-cover|--subject-cover|-e |--evalue|--max-target' \
  "$WS2/scripts"/3*.sh 2>/dev/null | head -40 | sed 's/^/    /'

sub "every per-million / normalisation expression"
grep -rn -E '1e6|1000000|per_million|perM|/ *READS|TOTAL_READS|DENOM|normal' \
  "$WS2/scripts"/3*.sh 2>/dev/null | head -30 | sed 's/^/    /'

sub "the read files the pipeline points at, and their actual read counts"
for f in "$WS2"/00_input/*.fastq.gz; do
  [ -s "$f" ] || continue
  printf '    %-44s %s bytes\n' "${f#$WS2/}" "$(stat -c %s "$f")"
done
echo "    (not counting reads here -- that is minutes of gunzip; the numbers"
echo "     already recorded are 7,528,854 sediment as used and 6,827,767 as"
echo "     step 34 denominated. Section A should say which file each used.)"

# ================================================================ D
say "D. THE EXISTING MARKER TABLE -- what mtmB has to sit beside"
for f in "$WS2"/34_unbinned/marker_levels.tsv "$WS2"/34_unbinned/VERDICT.txt \
         "$WS2"/34_unbinned/assembly_hits.tsv "$WS2"/31_readlevel/*.tsv; do
  [ -s "$f" ] && dump "$f" 80 10
done

sub "is K16176 anywhere in this workspace already?"
grep -rl 'K16176' "$WS2" --include='*.tsv' --include='*.txt' --include='*.sh' 2>/dev/null \
  | head -20 | sed 's/^/    /'
echo "    --- and the methylamine KOs for reference:"
echo "        K16176 mtmB  monomethylamine methyltransferase   <- the target"
echo "        K16177 mtmC  monomethylamine corrinoid protein   <- DECOY, not target"
echo "        K16178 mtbB  dimethylamine methyltransferase"
echo "        K14083 mttB  trimethylamine methyltransferase"
echo "        K04480 mtaB  methanol methyltransferase"
echo "        K00399 mcrA  the biological ceiling"

# ================================================================ E
say "E. FIGURE HOUSE STYLE -- hand it back in full"
for f in "$WS2"/scripts/etoliko_figstyle*.py "$WS2"/scripts/fig0*.py "$WS2"/figs/*.py; do
  [ -s "$f" ] && dump "$f" 500 0
done
sub "anything already rendered"
find "$WS2" -maxdepth 3 -type f \( -name 'fig*.pdf' -o -name 'fig*.png' -o -name 'fig*.svg' \) 2>/dev/null \
  | head -20 | while read -r f; do printf '    %-52s %s bytes\n' "${f#$WS2/}" "$(stat -c %s "$f")"; done

say "F. WHAT COMES NEXT"
cat <<'EOT'
  Section A tells me the exact search and normalisation rule, section B the
  reference-naming convention, section C the denominator. Script 61 will then
  run K16176 against the sediment and culture reads with:

    - MtmC (K16177) in the reference set as a DECOY, so reads that are really
      the corrinoid protein are pulled away from mtmB instead of counted as it
    - the same positive controls the other read searches used (rpoB, glnA, tuf,
      EF-1alpha), because a zero is only interpretable when a control passes
    - a promiscuity check on every reference, the dsrB lesson: two references
      were taking 55% and 40% of the sediment dsrB reads and inflated it 18-fold
    - the mcrA ceiling stated with the result. Every methanogen carries mcrA;
      only some carry mtmB. An mtmB count above the mcrA count is a failed
      search, not a finding.

  Section E gives me the figure house style so the five remaining panels match
  fig02 exactly rather than approximately.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
