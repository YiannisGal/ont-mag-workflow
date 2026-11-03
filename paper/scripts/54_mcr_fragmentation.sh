#!/usr/bin/env bash
# 54_mcr_fragmentation.sh -- is ETO_LIQ_MAG_017's mcr split by frameshifts?
#
#   bash scripts/run.sh scripts/54_mcr_fragmentation_XXXXXXX.sh
#
# WHY THIS EXISTS. Script 53 settled the nomination question: MAG_017 has no
# single ORF above K00399's curated threshold of 775.53, so step 15 was right
# not to offer it to the tree. But 53 also showed the hits arriving in pairs
# with CONSECUTIVE gene numbers on ONE contig -- 08195 + 08190, 11865 + 11870 --
# which looks like one gene broken into two ORFs by an internal frameshift, a
# known ONT artefact, rather than a genome that lacks mcr. That reading is
# currently an INFERENCE from gene numbering. This script tests it against
# coordinates.
#
# WHAT IT DOES. For all eight archaeal genomes, it takes every K00399/K00401/
# K00402 hit, looks the protein up in the CheckM2 prodigal proteome to get
# contig, start, end and strand, and asks per KO:
#   - is there a single ORF above threshold?          (intact)
#   - are there >=2 hits, same contig, same strand, adjacent gene numbers,
#     separated by a small intergenic gap?            (split candidate)
# It calibrates protein length against the above-threshold ORFs of the seven
# genomes that DID reach the tree, so "too short by half" is measured, not
# asserted.
#
# WHAT IT CANNOT DO. It cannot prove a frameshift. Coordinates + strand +
# adjacency + summed length near the reference length is strong circumstantial
# evidence; the proof would be aligning the two fragments to one reference MCR
# alpha and seeing them tile it end to end. If this script says SPLIT, that
# alignment is the follow-up and it is cheap.
#
# Read-only except 54_mcr_fragmentation/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/54_mcr_fragmentation"; mkdir -p "$OUT"
exec > >(tee "$OUT/54.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }

echo "54 -- mcr fragmentation test"
echo "date : $(date -Is)"
echo "out  : $OUT"

# bin -> MAG, the eight archaea. MAG_017 last so it reads as the test case.
BINS="LIQ_das_vamb_003:ETO_LIQ_MAG_003 \
LIQ_das_vamb_005:ETO_LIQ_MAG_005 \
LIQ_das_vamb_011:ETO_LIQ_MAG_007 \
LIQ_das_semibin_028:ETO_LIQ_MAG_010 \
LIQ_das_vamb_017:ETO_LIQ_MAG_011 \
LIQ_das_metabat_047:ETO_LIQ_MAG_016 \
LIQ_das_metabat_030:ETO_LIQ_MAG_018 \
LIQ_das_vamb_004:ETO_LIQ_MAG_017"

# ================================================================ A
say "A. INPUTS -- find them, do not assume them"

KOLIST="$(find "$WS2" -type f -name ko_list 2>/dev/null | head -1)"
if [ -n "$KOLIST" ]; then
  ok "ko_list: ${KOLIST#$WS2/}"
  for k in K00399 K00401 K00402; do
    line="$(awk -v k="$k" '$1==k{print; exit}' "$KOLIST")"
    if [ -n "$line" ]; then
      thr="$(printf '%s\n' "$line" | awk '{print $2+0}')"
      typ="$(printf '%s\n' "$line" | awk '{print $5}')"
      printf '         %s  threshold %-10s score_type %s\n' "$k" "$thr" "$typ"
      eval "THR_$k=$thr"
    else
      bad "$k not in ko_list -- falling back to the recorded value"
    fi
  done
else
  bad "no ko_list found; using the values recorded in script 53"
fi
THR_K00399="${THR_K00399:-775.53}"
THR_K00401="${THR_K00401:-504.30}"
THR_K00402="${THR_K00402:-360.93}"
printf '  using: K00399 %s  K00401 %s  K00402 %s\n' "$THR_K00399" "$THR_K00401" "$THR_K00402"

# hits: 45_kofam_full uses BIN|PRODIGAL_ID, which maps straight onto the
# CheckM2 proteome. That is the only id space in this workspace that carries
# coordinates, so this script requires it.
mapfile -t TBLS < <(find "$WS2/45_kofam_full" -type f -name '*.tbl' 2>/dev/null | sort)
if [ "${#TBLS[@]}" -eq 0 ]; then
  bad "no .tbl under 45_kofam_full -- cannot map hits to coordinates."
  echo "STEP_STATUS: FAILED"; exit 1
fi
ok "${#TBLS[@]} KOfam table(s) under 45_kofam_full"

# CheckM2 proteomes
PROTDIR="$(find "$WS2" -type d -name protein_files 2>/dev/null | head -1)"
if [ -z "$PROTDIR" ]; then
  bad "no protein_files directory found under $WS2"
  echo "STEP_STATUS: FAILED"; exit 1
fi
ok "proteomes: ${PROTDIR#$WS2/}  ($(find "$PROTDIR" -name '*.faa' | wc -l) faa)"

# ================================================================ B
say "B. COLLECT EVERY mcr HIT WITH ITS COORDINATES"

HITS="$OUT/hits_with_coords.tsv"
: > "$HITS"
printf 'mag\tbin\tko\tprotein\tscore\tevalue\tcontig\tgene_idx\tstart\tend\tstrand\taa_len\n' > "$HITS"

for pair in $BINS; do
  BIN="${pair%%:*}"; MAG="${pair##*:}"

  FAA="$(find "$PROTDIR" -name "${BIN}*.faa" 2>/dev/null | head -1)"
  if [ -z "$FAA" ]; then
    bad "$MAG ($BIN): no proteome file -- skipped"
    continue
  fi

  # protein -> contig, idx, start, end, strand
  grep '^>' "$FAA" | awk -F' # ' '
    { id=substr($1,2); s=$2+0; e=$3+0; st=$4+0
      n=split(id,a,"_"); idx=a[n]+0
      ctg=id; sub(/_[0-9]+$/,"",ctg)
      printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\n", id, ctg, idx, s, e, st, int((e-s+1)/3)-1 }' \
    > "$OUT/.coords.$BIN"

  for T in "${TBLS[@]}"; do
    awk -v bin="$BIN" -v mag="$MAG" -v cf="$OUT/.coords.$BIN" '
      BEGIN{ while((getline l < cf)>0){ split(l,c,"\t"); C[c[1]]=l } }
      /^#/{next}
      { ko=$3
        if(ko!="K00399" && ko!="K00401" && ko!="K00402") next
        t=$1
        p=index(t,"|")
        if(p==0) next                       # not the prefixed id space
        b=substr(t,1,p-1); prot=substr(t,p+1)
        if(b!=bin) next
        if(!(prot in C)) { printf "%s\t%s\t%s\t%s\t%.1f\t%s\tNOCOORD\t-\t-\t-\t-\t-\n", mag,bin,ko,prot,$6+0,$5; next }
        split(C[prot],c,"\t")
        printf "%s\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               mag,bin,ko,prot,$6+0,$5,c[2],c[3],c[4],c[5],(c[6]>0?"+":"-"),c[7] }' "$T"
  done
done >> "$HITS"

rm -f "$OUT"/.coords.*

N=$(($(wc -l < "$HITS") - 1))
if [ "$N" -le 0 ]; then
  bad "no mcr hits matched any of the eight bins. Check the id space in 45_kofam_full."
  echo "STEP_STATUS: FAILED"; exit 1
fi
ok "$N mcr hits across the eight archaea, with coordinates"
echo "  -> ${HITS#$WS2/}"

# ================================================================ C
say "C. REFERENCE LENGTH -- what does an intact, above-threshold mcr ORF look like?"

awk -F'\t' -v t399="$THR_K00399" -v t401="$THR_K00401" -v t402="$THR_K00402" '
  NR==1{next}
  $12=="-"{next}
  { thr = ($3=="K00399"? t399+0 : ($3=="K00401"? t401+0 : t402+0))
    if($5+0 < thr) next
    n[$3]++; s[$3]+=$12+0
    if(mn[$3]==""||$12+0<mn[$3]) mn[$3]=$12+0
    if($12+0>mx[$3]) mx[$3]=$12+0 }
  END{ printf "  %-8s %5s  %8s  %8s  %8s\n","KO","n","mean_aa","min_aa","max_aa"
       for(k in n) printf "  %-8s %5d  %8.1f  %8d  %8d\n", k, n[k], s[k]/n[k], mn[k], mx[k] }' \
  "$HITS" | sort

# ================================================================ D
say "D. PER GENOME, PER KO -- intact, split, or absent"

awk -F'\t' -v t399="$THR_K00399" -v t401="$THR_K00401" -v t402="$THR_K00402" '
  NR==1{next}
  { key=$1"\t"$3
    n[key]++
    prot[key,n[key]]=$4; sc[key,n[key]]=$5+0
    ctg[key,n[key]]=$7;  idx[key,n[key]]=$8+0
    st[key,n[key]]=$9+0; en[key,n[key]]=$10+0
    sd[key,n[key]]=$11;  ln[key,n[key]]=$12+0
    mag[key]=$1; ko[key]=$3
    if(!(($1) in seen)){ order[++nm]=$1; seen[$1]=1 } }
  END{
    for(m=1;m<=nm;m++){
      M=order[m]; printed=0
      split("K00399 K00401 K00402",K," ")
      for(ki=1;ki<=3;ki++){
        k=K[ki]; key=M"\t"k
        if(!(key in n)) continue
        if(!printed){ printf "\n%s\n", M; printed=1 }
        thr = (k=="K00399"? t399+0 : (k=="K00401"? t401+0 : t402+0))

        # sort the hits for this KO by contig then gene index (insertion sort)
        c=n[key]
        for(i=1;i<=c;i++) o[i]=i
        for(i=2;i<=c;i++){ v=o[i]; j=i-1
          while(j>=1 && (ctg[key,o[j]] > ctg[key,v] || (ctg[key,o[j]]==ctg[key,v] && idx[key,o[j]]>idx[key,v]))){ o[j+1]=o[j]; j-- }
          o[j+1]=v }

        best=0; nab=0
        for(i=1;i<=c;i++){ if(sc[key,i]>best) best=sc[key,i]; if(sc[key,i]>=thr) nab++ }

        printf "  %-8s threshold %-8.2f  hits %d  best %.1f  above %d\n", k, thr, c, best, nab
        for(i=1;i<=c;i++){ j=o[i]
          printf "      %-26s %7.1f  %-18s idx %-6d %9d..%-9d %s  %4d aa\n",
                 prot[key,j], sc[key,j], ctg[key,j], idx[key,j], st[key,j], en[key,j], sd[key,j], ln[key,j] }

        # adjacency scan
        ngrp=0
        for(i=1;i<c;i++){ a=o[i]; b=o[i+1]
          if(ctg[key,a]!=ctg[key,b]) continue
          if(idx[key,b]-idx[key,a]!=1) continue
          if(sd[key,a]!=sd[key,b]) continue
          gap = st[key,b] - en[key,a] - 1
          ngrp++
          printf "      >> ADJACENT: %s + %s  same contig, same strand %s, gap %d bp, sum score %.1f, sum %d aa\n",
                 prot[key,a], prot[key,b], sd[key,a], gap, sc[key,a]+sc[key,b], ln[key,a]+ln[key,b]
          if(sc[key,a]+sc[key,b] >= thr) printf "         (summed score clears the threshold)\n"
          else printf "         (summed score %.1f still short of %.2f by %.1f)\n", sc[key,a]+sc[key,b], thr, thr-(sc[key,a]+sc[key,b])
        }

        if(nab>0)            printf "      VERDICT: INTACT (%d ORF(s) above threshold)\n", nab
        else if(ngrp>0)      printf "      VERDICT: SPLIT  (%d adjacent same-strand pair(s), none intact)\n", ngrp
        else                 printf "      VERDICT: SUB-THRESHOLD, not adjacent -- probably not mcr at all\n"
      }
    }
  }' "$HITS"

# ================================================================ E
say "E. WHAT THIS MEANS"
cat <<'EOT'
  Read section D one genome at a time.

  INTACT on K00399 means the genome has a full-length mcrA and belongs in the
  tree. Seven genomes should read that way; they are the control.

  SPLIT on K00399 for ETO_LIQ_MAG_017 means the gene is there but broken into
  two ORFs. Check three things in its ">> ADJACENT" lines: the gap should be
  small (a frameshift leaves the fragments nearly touching, often overlapping,
  so a negative or near-zero gap is the expected signal); the summed aa length
  should approach the reference mean in section C; and the summed score should
  approach the threshold. All three together, and the report's "fragmented, not
  missing" wording is evidence rather than inference.

  If instead the gap is thousands of bp, or the summed length is far short of
  the reference, the pairs are two separate partial genes and the honest
  statement is weaker: MAG_017 carries mcr-like sequence that no single ORF
  supports as a full-length gene.

  EITHER WAY the nomination decision in step 15 stands and the tree does not
  need rebuilding. This only decides what sentence goes in the chapter.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
