#!/usr/bin/env bash
# 53_mcra_nomination.sh -- find MAG_017's mcr genes properly, and find out why
#                          step 15 never nominated them.
#
#   bash scripts/53_mcra_nomination.sh
#
# WHY THIS EXISTS. Script 52 established the important thing: ETO_LIQ_MAG_017 was
# NEVER OFFERED to the mcrA tree, so nothing about its mcrA has been tested by
# placement. But 52's section B then reported "no mcr hits found for this genome"
# -- and section B3 reported the same for LIQ_das_vamb_003, a genome we know
# carries mcrA and that IS in the tree. A search that returns zero for a positive
# control is a broken search, not a result. 52 section B is WITHDRAWN.
#
# THE BUG. 52 built one protein-ID list from the CheckM2 proteome
# (06_binqc/checkm2/protein_files/*.faa, ids like LIQ_contig_20103_12) and
# matched it against every .tbl it could find. But the tree's own queries are
# named ETOLIQMAG003_04260 -- Bakta locus tags from 09_annotation. The two KOfam
# runs in this workspace used DIFFERENT proteomes with DIFFERENT id spaces, and
# one id list cannot match both. This script discovers each file's id space
# before matching anything, and carries three id sets per genome.
#
# WHAT IT DECIDES. Whether MAG_017's mcr proteins clear K00399's curated
# threshold of 775.53. If they do, the nomination step missed them and the tree
# should be rebuilt. If they do not, their absence from the tree is correct --
# and the "4 mcrA copies" figure in 41_env_checks is a sub-threshold count that
# must not be quoted as genes.
#
# Read-only except 53_mcra_nomination/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/53_mcra_nomination"; mkdir -p "$OUT"
exec > >(tee "$OUT/53.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }

echo "53 -- mcrA nomination audit"
echo "date : $(date -Is)"
echo "out  : $OUT"

TARGET_MAG=ETO_LIQ_MAG_017
TARGET_BIN=LIQ_das_vamb_004
LOCUS=ETOLIQMAG017_          # Bakta locus prefix for this MAG
# the seven genomes that DID reach the tree, as calibration
CAL_BINS="LIQ_das_vamb_003 LIQ_das_vamb_011 LIQ_das_semibin_028 LIQ_das_vamb_005 LIQ_das_vamb_017 LIQ_das_metabat_047 LIQ_das_metabat_030"
CAL_LOCI="ETOLIQMAG003_ ETOLIQMAG005_ ETOLIQMAG007_ ETOLIQMAG010_ ETOLIQMAG011_ ETOLIQMAG016_ ETOLIQMAG018_"

# ================================================================ A
say "A. WHAT ID SPACE DOES EACH KOfam OUTPUT ACTUALLY USE?"
echo "Before matching anything, look. This is the step 52 skipped."

find "$WS2/11_kofam" "$WS2/45_kofam_full" -type f -name '*.tbl' 2>/dev/null |
  sort > "$OUT/tbl_files.txt"
ok "$(wc -l < "$OUT/tbl_files.txt") .tbl file(s) found"

sub "A1. first data row of each .tbl, with its column count"
while read -r T; do
  [ -s "$T" ] || { printf '  %-58s EMPTY\n' "${T#$WS2/}"; continue; }
  R=$(grep -v '^#' "$T" 2>/dev/null | head -1)
  [ -z "$R" ] && { printf '  %-58s no data rows\n' "${T#$WS2/}"; continue; }
  printf '\n  %s\n' "${T#$WS2/}"
  printf '    cols=%s\n' "$(echo "$R" | awk '{print NF}')"
  printf '    $1=%s\n    $3=%s\n    $5=%s  $6=%s\n' \
    "$(echo "$R" | awk '{print $1}')" "$(echo "$R" | awk '{print $3}')" \
    "$(echo "$R" | awk '{print $5}')" "$(echo "$R" | awk '{print $6}')"
done < "$OUT/tbl_files.txt"

sub "A2. header style of the two proteome sets"
for f in "$WS2/09_annotation/$TARGET_BIN/$TARGET_BIN.faa" \
         "$WS2/06_binqc/checkm2/protein_files/$TARGET_BIN.faa" \
         "$WS2/45_kofam_full/all_proteins.faa"; do
  if [ -s "$f" ]; then
    printf '\n  %s\n' "${f#$WS2/}"
    grep -m3 '^>' "$f" | sed 's/^/    /'
  else printf '\n  %s  [MISSING]\n' "${f#$WS2/}"; fi
done

# ================================================================ B
say "B. WHERE DID STEP 15's QUERIES COME FROM?"
CAND=$(find "$WS2" -maxdepth 4 -type f -name 'mcrA_K00399_candidates.faa' 2>/dev/null | head -1)
if [ -n "$CAND" ]; then
  ok "candidate file: ${CAND#$WS2/}  ($(grep -c '^>' "$CAND") seqs)"
  sub "B1. which genomes are in it"
  grep '^>' "$CAND" | sed 's/^>//' | awk -F'|' '{print $1"|"$2}' | sort | uniq -c | sort -rn
  sub "B2. is $TARGET_MAG in it?"
  if grep -q -e "$TARGET_MAG" -e "$TARGET_BIN" -e "$LOCUS" "$CAND"; then
    ok "PRESENT in the candidate file -- so the loss happened AFTER nomination"
    grep '^>' "$CAND" | grep -e "$TARGET_MAG" -e "$LOCUS" | sed 's/^/    /'
  else
    no "ABSENT from the candidate file -- the loss happened AT nomination"
  fi
else
  no "mcrA_K00399_candidates.faa not found; searching for how it was made"
fi

sub "B3. the lines in 11_metabolism.sh / 11b_kofam_archaea.sh that build it"
for s in "$WS2/scripts/11b_kofam_archaea.sh" "$WS2/scripts/11_metabolism.sh"; do
  [ -s "$s" ] || continue
  echo "  -- $(basename "$s")"
  grep -n -E 'candidates|K00399|threshold|thresh|awk .*\$6|--cut|score' "$s" | head -18 | sed 's/^/    /'
done

# ================================================================ C
say "C. THE MCR HITS FOR $TARGET_MAG, MATCHED IN EVERY ID SPACE"

# build three id sets
: > "$OUT/ids_all.txt"
BAK="$WS2/09_annotation/$TARGET_BIN/$TARGET_BIN.faa"
CK="$WS2/06_binqc/checkm2/protein_files/$TARGET_BIN.faa"
for f in "$BAK" "$CK"; do
  [ -s "$f" ] && grep '^>' "$f" | sed 's/^>//' | awk '{print $1}' >> "$OUT/ids_all.txt"
done
# contig names, so a prefixed id still matches
if [ -s "$CK" ]; then
  grep '^>' "$CK" | sed 's/^>//' | awk '{print $1}' | sed 's/_[0-9]*$//' | sort -u > "$OUT/contigs.txt"
fi
sort -u "$OUT/ids_all.txt" -o "$OUT/ids_all.txt"
ok "id set: $(wc -l < "$OUT/ids_all.txt") protein ids, $( [ -s "$OUT/contigs.txt" ] && wc -l < "$OUT/contigs.txt" || echo 0) contigs, locus prefix $LOCUS"

sub "C1. every K00399 / K00401 / K00402 row that belongs to this genome"
printf '%-34s %-8s %12s %10s  %s\n' TARGET KO EVALUE SCORE FILE
while read -r T; do
  [ -s "$T" ] || continue
  awk -v ids="$OUT/ids_all.txt" -v ctg="$OUT/contigs.txt" -v loc="$LOCUS" -v src="${T#$WS2/}" '
    BEGIN{ while((getline l < ids)>0) keep[l]=1
           while((getline c < ctg)>0) cg[c]=1 }
    /^#/{next}
    {
      ko=$3; if(ko!="K00399" && ko!="K00401" && ko!="K00402") next
      t=$1; hit=0
      if(t in keep) hit=1
      else if(index(t,loc)==1) hit=1
      else { for(c in cg) if(index(t,c)>0) { hit=1; break } }
      if(!hit) next
      printf "%-34s %-8s %12s %10s  %s\n", t, ko, $5, $6, src
    }' "$T"
done < "$OUT/tbl_files.txt" | sort -k2,2 -k4,4gr | tee "$OUT/mcr_hits_017.tsv"
N017=$(wc -l < "$OUT/mcr_hits_017.tsv")
echo; ok "$N017 row(s) for $TARGET_MAG"

sub "C2. the same three KOs in the seven genomes that DID reach the tree (calibration)"
printf '%-24s %-10s %10s %10s %s\n' GENOME KO BEST_SCORE THRESHOLD VERDICT
set -- $CAL_LOCI
for b in $CAL_BINS; do
  L="$1"; shift
  bf="$WS2/09_annotation/$b/$b.faa"; cf="$WS2/06_binqc/checkm2/protein_files/$b.faa"
  : > "$OUT/.ids"; for f in "$bf" "$cf"; do [ -s "$f" ] && grep '^>' "$f" | sed 's/^>//' | awk '{print $1}' >> "$OUT/.ids"; done
  [ -s "$cf" ] && grep '^>' "$cf" | sed 's/^>//' | awk '{print $1}' | sed 's/_[0-9]*$//' | sort -u > "$OUT/.ctg" || : > "$OUT/.ctg"
  BEST=$(while read -r T; do
      [ -s "$T" ] || continue
      awk -v ids="$OUT/.ids" -v ctg="$OUT/.ctg" -v loc="$L" '
        BEGIN{ while((getline l < ids)>0) keep[l]=1; while((getline c < ctg)>0) cg[c]=1 }
        /^#/{next}
        { if($3!="K00399") next
          t=$1; hit=0
          if(t in keep) hit=1; else if(index(t,loc)==1) hit=1
          else { for(c in cg) if(index(t,c)>0){hit=1;break} }
          if(hit && $6+0>m) m=$6+0 }
        END{ if(m>0) printf "%.1f", m }' "$T"
    done < "$OUT/tbl_files.txt" | sort -rn | head -1)
  [ -z "$BEST" ] && BEST=0
  V=$(awk -v s="$BEST" 'BEGIN{print (s>=775.53)?"above threshold":"BELOW 775.53"}')
  printf '%-24s %-10s %10s %10s %s\n' "$b" K00399 "$BEST" 775.53 "$V"
done
rm -f "$OUT/.ids" "$OUT/.ctg"

BEST017=$(awk '$2=="K00399"{if($4+0>m)m=$4+0}END{printf "%.1f", m+0}' "$OUT/mcr_hits_017.tsv")
printf '%-24s %-10s %10s %10s %s\n' "$TARGET_BIN" K00399 "$BEST017" 775.53 \
  "$(awk -v s="$BEST017" 'BEGIN{print (s>=775.53)?"above threshold":"BELOW 775.53"}')"

# ================================================================ D
say "D. WHERE 41_env_checks GOT ITS 'K00399 = 4'"
AK="$WS2/41_env_checks/archaea_kofam_markers.tsv"
if [ -s "$AK" ]; then
  echo "  the row as written:"
  head -1 "$AK" | sed 's/^/    /'; grep "$TARGET_BIN" "$AK" | sed 's/^/    /'
  echo
  echo "  how 41 counted (from its script):"
  for s in "$WS2/scripts/41_env_and_checks.sh" "$WS2/scripts/41b_fixes.sh" "$WS2/scripts/41f_smoke_and_kofam.sh"; do
    [ -s "$s" ] || continue
    grep -n -E 'archaea_kofam_markers|K00399|threshold|\$6' "$s" | head -8 | sed "s|^|    $(basename "$s"):|"
  done
else
  no "$AK missing"
fi

# ================================================================ E
say "E. VERDICT"
if [ "$N017" -eq 0 ]; then
  bad "still zero rows for $TARGET_MAG after matching all three id spaces."
  echo "   If section C2 DID find scores for the seven tree genomes, then the"
  echo "   search works and this genome genuinely has no K00399 hit in these"
  echo "   tables -- in which case 41_env_checks' 'K00399 = 4' came from"
  echo "   somewhere else and section D says where."
  echo "   If C2 is also empty, the search is still wrong and nothing is decided."
elif awk -v s="$BEST017" 'BEGIN{exit !(s>=775.53)}'; then
  ok "$TARGET_MAG has K00399 hit(s) ABOVE the curated threshold (best $BEST017 >= 775.53)."
  echo "   >> Step 15 should have nominated it and did not. The tree is"
  echo "      incomplete and REBUILD is justified:"
  echo "         REBUILD=1 bash scripts/52_mcra_recheck_03732e7.sh"
  echo "      (52 will now find the proteins, because this script wrote"
  echo "       $OUT/mcr_hits_017.tsv with the correct ids.)"
else
  ok "$TARGET_MAG's best K00399 score is $BEST017, BELOW the curated 775.53."
  echo "   >> Its absence from the tree is CORRECT, not an omission. And the"
  echo "      '4 mcrA copies' in 41_env_checks is a count of sub-threshold hits."
  echo "      Two consequences for the chapter:"
  echo "        1. mcr is confirmed in 7 archaeal genomes, not 8. MASTER's mcr"
  echo "           column for this genome counts hits, not genes above threshold."
  echo "        2. Alkyl-CoM reductase scores against K00399 too. A sub-threshold"
  echo "           K00399 hit in a methanogen genome is exactly the ACR/MCR"
  echo "           ambiguity step 15 was built to avoid."
fi
hr; echo "STEP_STATUS: COMPLETE  $(date -Is)"
