#!/usr/bin/env bash
# 52_mcra_recheck.sh -- why is ETO_LIQ_MAG_017 missing from the mcrA tree, and add it.
#
#   bash scripts/52_mcra_recheck.sh            # diagnose only
#   REBUILD=1 bash scripts/52_mcra_recheck.sh  # diagnose, then rebuild the tree
#
# THE PROBLEM. Step 15 placed mcrA proteins against a curated reference set and
# returned CONFIRMED for 9 proteins across 7 genomes. ETO_LIQ_MAG_017
# (LIQ_das_vamb_004, Methanogenium) is not among them -- yet the genome-wide
# KOfam run calls FOUR K00399 copies in it, and MASTER.tsv carries mcr = 3.
# Either step 15 never offered those proteins to the tree, or it offered them
# and they failed. Those are different facts and the chapter needs the right one.
#
# WHY DIAGNOSE BEFORE FIXING. Adding a sequence to a tree and reading the answer
# you wanted is how a false positive gets into a thesis. So section A asks what
# step 15 actually did, section B asks whether the proteins are real by their
# own KO's curated threshold, and the rebuild in section C only runs if you ask
# for it -- and it writes to a NEW directory. 15_marker_trees is never touched.
#
# SAFETY. Read-only except 52_mcra_recheck/. No conda env is modified; tools are
# invoked through the project's own run_tool resolver.

set -uo pipefail

WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=/path/to/workspace_2"; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/52_mcra_recheck"
MT="$WS2/15_marker_trees/mcrA"
TARGET_BIN="LIQ_das_vamb_004"
TARGET_MAG="ETO_LIQ_MAG_017"
mkdir -p "$OUT"
LOG="$OUT/52.log"
exec > >(tee "$LOG") 2>&1

hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){  printf '  [OK]   %s\n' "$*"; }
no(){  printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }

echo "52 -- mcrA re-check for $TARGET_MAG ($TARGET_BIN)"
echo "date : $(date -Is)"
echo "ws2  : $WS2"
echo "out  : $OUT   (nothing outside this directory is written)"

# ---------------------------------------------------------------- config
for f in 00_config.sh 00_config_patch.sh env_resolved.sh; do
  [ -f "$WS2/$f" ] && . "$WS2/$f" >/dev/null 2>&1
done
have(){ command -v "$1" >/dev/null 2>&1; }
tool(){  # tool <name> <args...>   -- project resolver first, PATH second
  local t="$1"; shift
  if declare -F run_tool >/dev/null 2>&1; then run_tool "$t" "$@" && return 0; fi
  if have "$t"; then "$t" "$@"; return $?; fi
  return 127
}

# ================================================================ A
say "A. WHAT DID STEP 15 ACTUALLY PUT IN THE mcrA TREE?"

if [ ! -s "$MT/mcrA_all.faa" ]; then
  bad "no $MT/mcrA_all.faa -- cannot tell what step 15 used."
  echo "STEP_STATUS: FAILED"; exit 1
fi

NQ=$(grep -c '^>' "$MT/mcrA_all.faa")
ok "mcrA_all.faa: $NQ sequences (queries + references)"

sub "A1. every NON-reference header (these are the catalogue's own proteins)"
grep '^>' "$MT/mcrA_all.faa" | grep -v '^>MCR|' | sed 's/^>//' | sort
NOWN=$(grep '^>' "$MT/mcrA_all.faa" | grep -vc '^>MCR|')
echo
echo "  catalogue proteins offered to the tree: $NOWN"
echo "  proteins that came back CONFIRMED     : $( [ -s "$MT/mcrA_verdicts.tsv" ] && echo $(( $(wc -l < "$MT/mcrA_verdicts.tsv") - 1 )) || echo '?' )"

sub "A2. is $TARGET_MAG / $TARGET_BIN in that file at all?"
if grep -q -e "$TARGET_MAG" -e "$TARGET_BIN" "$MT/mcrA_all.faa"; then
  ok "PRESENT -- step 15 DID offer it to the tree."
  echo "  its headers:"
  grep '^>' "$MT/mcrA_all.faa" | grep -e "$TARGET_MAG" -e "$TARGET_BIN" | sed 's/^/    /'
  echo
  echo "  >> Then its absence from mcrA_verdicts.tsv is a VERDICT, not an omission:"
  echo "     the placement ran and the k-nearest rule did not return CONFIRMED."
  echo "     Section C will show where it actually sits in the tree."
  OFFERED=yes
else
  no "ABSENT -- step 15 never offered it to the tree."
  echo
  echo "  >> So this is a NOMINATION failure, not a phylogenetic result. Nothing"
  echo "     about this genome's mcrA has been tested by placement, and the"
  echo "     chapter must not say it was."
  OFFERED=no
fi

sub "A3. which genomes ARE represented, counted"
grep '^>' "$MT/mcrA_all.faa" | grep -v '^>MCR|' | sed 's/^>//' |
  awk -F'|' '{print $1}' | sort | uniq -c | sort -rn

sub "A4. how step 15 chose its queries"
for s in "$WS2/scripts/15_marker_trees.sh"; do
  [ -s "$s" ] || { no "$s not found"; continue; }
  echo "  nomination lines from $(basename "$s"):"
  grep -n -i -E 'mcrA|K00399|hmmsearch|product|grep .*methyl|query' "$s" | head -25 | sed 's/^/    /'
done

# ================================================================ B
say "B. ARE THE FOUR K00399 COPIES REAL, BY THE KO's OWN THRESHOLD?"

FAA=""
for c in "$WS2/06_binqc/checkm2/protein_files/$TARGET_BIN.faa" \
         "$WS2/09_annotation/$TARGET_BIN/$TARGET_BIN.faa"; do
  [ -s "$c" ] && { FAA="$c"; break; }
done
if [ -z "$FAA" ]; then
  bad "no proteome found for $TARGET_BIN -- looked in 06_binqc/checkm2/protein_files and 09_annotation"
  echo "STEP_STATUS: FAILED"; exit 1
fi
ok "proteome: $FAA  ($(grep -c '^>' "$FAA") proteins)"

grep '^>' "$FAA" | sed 's/^>//' | awk '{print $1}' | sort > "$OUT/ids_$TARGET_BIN.txt"

sub "B1. K00399 / K00401 / K00402 hits belonging to this genome"
printf '%-28s %-8s %10s %10s  %s\n' TARGET KO EVALUE SCORE SOURCE > "$OUT/mcr_hits.tsv"
FOUND=0
while read -r TBL; do
  [ -s "$TBL" ] || continue
  awk -v ids="$OUT/ids_$TARGET_BIN.txt" -v src="$TBL" '
    BEGIN{ while((getline l < ids)>0) keep[l]=1 }
    /^#/ {next}
    { ko=$3; if(ko!="K00399" && ko!="K00401" && ko!="K00402") next
      if(!( $1 in keep)) next
      printf "%-28s %-8s %10s %10s  %s\n", $1, ko, $5, $6, src }
  ' "$TBL" >> "$OUT/mcr_hits.tsv"
done < <(find "$WS2/11_kofam" "$WS2/45_kofam_full" -type f -name '*.tbl' -size -2G 2>/dev/null | sort)

cat "$OUT/mcr_hits.tsv"
FOUND=$(( $(wc -l < "$OUT/mcr_hits.tsv") - 1 ))
echo
if [ "$FOUND" -gt 0 ]; then ok "$FOUND mcr hit row(s) for this genome"
else bad "no mcr hits found for this genome in any KOfam .tbl -- then MASTER's mcr=3 needs explaining"; fi

sub "B2. the curated thresholds those scores must clear"
KOLIST=$(find "$WS2" -maxdepth 4 -type f -name 'ko_list*' 2>/dev/null | head -1)
if [ -n "$KOLIST" ]; then
  ok "ko_list: $KOLIST"
  awk -F'\t' 'NR==1 || $1=="K00399" || $1=="K00401" || $1=="K00402"' "$KOLIST" | cut -f1-5 | sed 's/^/    /'
else
  no "no ko_list found -- compare the scores above against the other genomes' instead:"
fi

sub "B3. the same three KOs in the seven genomes that DID make the tree"
printf '%-26s %-8s %10s\n' GENOME KO BEST_SCORE
for g in LIQ_das_vamb_003 LIQ_das_vamb_011 LIQ_das_semibin_028 LIQ_das_vamb_005 \
         LIQ_das_vamb_017 LIQ_das_metabat_047 LIQ_das_metabat_030 "$TARGET_BIN"; do
  gf="$WS2/06_binqc/checkm2/protein_files/$g.faa"
  [ -s "$gf" ] || continue
  grep '^>' "$gf" | sed 's/^>//' | awk '{print $1}' | sort > "$OUT/.ids.tmp"
  find "$WS2/11_kofam" "$WS2/45_kofam_full" -type f -name '*.tbl' -size -2G 2>/dev/null |
  while read -r TBL; do
    awk -v ids="$OUT/.ids.tmp" '
      BEGIN{ while((getline l < ids)>0) keep[l]=1 }
      /^#/{next}
      { if($3!="K00399") next; if(!($1 in keep)) next; if($6+0>m) m=$6+0 }
      END{ if(m>0) print m }' "$TBL"
  done | sort -rn | head -1 | awk -v g="$g" '{printf "%-26s %-8s %10s\n", g, "K00399", $1}'
done
rm -f "$OUT/.ids.tmp"

# ================================================================ C
say "C. REBUILD -- only if you asked for it"

if [ "${REBUILD:-0}" != "1" ]; then
  echo "  Not rebuilding. Read sections A and B first, then:"
  echo "     REBUILD=1 bash scripts/52_mcra_recheck.sh"
  echo
  echo "  If A2 said PRESENT, a rebuild adds nothing -- the placement already ran"
  echo "  and the answer is in the tree. If A2 said ABSENT and B shows the four"
  echo "  proteins clear their KO's curated threshold, the rebuild is the right fix."
  hr; echo "STEP_STATUS: DIAGNOSED"; exit 0
fi

if [ "$OFFERED" = "yes" ]; then
  bad "refusing to rebuild: $TARGET_MAG was ALREADY in the tree (section A2)."
  echo "  Rebuilding would not change its placement, only re-run it. Read the"
  echo "  existing tree instead: $MT/mcrA.treefile"
  hr; echo "STEP_STATUS: NO-ACTION"; exit 0
fi
if [ "$FOUND" -eq 0 ]; then
  bad "refusing to rebuild: no mcr proteins found for this genome (section B1)."
  hr; echo "STEP_STATUS: NO-ACTION"; exit 0
fi

sub "C1. pull the K00399 proteins and label them step-15 style"
awk 'NR>1 && $2=="K00399"{print $1}' "$OUT/mcr_hits.tsv" | sort -u > "$OUT/add_ids.txt"
NADD=$(wc -l < "$OUT/add_ids.txt")
ok "$NADD protein(s) to add"; sed 's/^/    /' "$OUT/add_ids.txt"

awk -v ids="$OUT/add_ids.txt" -v mag="$TARGET_MAG" -v bin="$TARGET_BIN" '
  BEGIN{ while((getline l < ids)>0) keep[l]=1 }
  /^>/{ id=substr($1,2); p=(id in keep)
        if(p) print ">" mag "|" bin "|" id; next }
  p{print}
' "$FAA" > "$OUT/add.faa"
ok "wrote $OUT/add.faa ($(grep -c '^>' "$OUT/add.faa") seqs)"

cat "$MT/mcrA_all.faa" "$OUT/add.faa" > "$OUT/mcrA_all_plus.faa"
ok "combined set: $(grep -c '^>' "$OUT/mcrA_all_plus.faa") sequences"

sub "C2. align, trim, place -- same tools step 15 used"
tool mafft --auto --thread 8 "$OUT/mcrA_all_plus.faa" > "$OUT/mcrA_aln.faa" 2> "$OUT/mafft.err"
[ -s "$OUT/mcrA_aln.faa" ] && ok "mafft: $(grep -c '^>' "$OUT/mcrA_aln.faa") seqs aligned" || { bad "mafft failed -- see $OUT/mafft.err"; tail -5 "$OUT/mafft.err"; echo "STEP_STATUS: FAILED"; exit 1; }

tool trimal -in "$OUT/mcrA_aln.faa" -out "$OUT/mcrA_trim.faa" -gt 0.5 2> "$OUT/trimal.err"
[ -s "$OUT/mcrA_trim.faa" ] && ok "trimal done" || { no "trimal failed, using the untrimmed alignment"; cp "$OUT/mcrA_aln.faa" "$OUT/mcrA_trim.faa"; }

tool iqtree3 -s "$OUT/mcrA_trim.faa" -m LG+G4 -B 1000 -T 8 --prefix "$OUT/mcrA" -redo > "$OUT/iqtree.out" 2>&1
if [ -s "$OUT/mcrA.treefile" ]; then ok "tree: $OUT/mcrA.treefile"
else bad "iqtree3 failed"; tail -15 "$OUT/iqtree.out"; echo "STEP_STATUS: FAILED"; exit 1; fi

sub "C3. re-run step 15's OWN verdict rule -- not a new one"
VS="$WS2/scripts/15b_tree_verdicts.py"
if [ ! -s "$VS" ]; then
  bad "$VS not found. The tree is built but NO verdict is issued:"
  echo "     inventing a different k-nearest rule here would not be comparable"
  echo "     to the other seven genomes. Locate 15b_tree_verdicts.py and re-run."
  hr; echo "STEP_STATUS: TREE-ONLY"; exit 0
fi
echo "  usage line(s) from $(basename "$VS"):"
grep -n -E 'argv|argparse|add_argument|usage' "$VS" | head -12 | sed 's/^/    /'
PY=$(command -v python3 || command -v python)
echo
echo "  attempting: $PY $VS $OUT/mcrA.treefile"
"$PY" "$VS" "$OUT/mcrA.treefile" > "$OUT/verdicts_new.tsv" 2> "$OUT/verdicts.err"
if [ -s "$OUT/verdicts_new.tsv" ]; then
  ok "new verdicts written"
else
  no "that call signature did not work -- see $OUT/verdicts.err"
  head -5 "$OUT/verdicts.err" | sed 's/^/    /'
  echo "     The tree is built and correct; only the verdict call needs the right"
  echo "     arguments (see the usage lines above)."
  hr; echo "STEP_STATUS: TREE-ONLY"; exit 0
fi

sub "C4. OLD vs NEW -- nothing may change silently"
echo "  --- step 15 verdicts (unchanged file) ---"
cat "$MT/mcrA_verdicts.tsv"
echo
echo "  --- rebuilt verdicts ---"
cat "$OUT/verdicts_new.tsv"
echo
echo "  Any genome whose verdict MOVED between those two tables is a change this"
echo "  rebuild caused, and has to be explained before either table is quoted."

hr; echo "STEP_STATUS: COMPLETE  $(date -Is)"
