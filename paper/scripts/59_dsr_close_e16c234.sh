#!/usr/bin/env bash
# 59_dsr_close.sh -- the last three questions on the dsrAB contigs.
#
#   bash scripts/run.sh scripts/59_dsr_close_XXXXXXX.sh
#
# WHAT 58 SETTLED. Both contigs carry REDUCTIVE-type dsrAB -- no query went
# anywhere near the Allochromatium / Thiobacillus / Chlorobium references, so
# these are sulfate reducers and not sulfur oxidisers. That had never been
# tested and it is the chapter's central premise.
#
#   LIQ_contig_12880  DsrB 89.5% and DsrA 88.5% to Desulfosporosinus meridiei
#                     DSM 13257, full length, Desulfosporosinus alone in the
#                     smallest clade. Genus-level certain.
#   LIQ_contig_19299  DsrB 78.0% to three Desulforamulus species, DsrA 82.9%
#                     to Desulfofalx alkaliphila. Desulfotomaculales.
#
# WHAT IS STILL OPEN, and it is one question: is contig 19299 the missing
# operon of ETO_LIQ_MAG_022 (LIQ_das_metabat_022, GTDB f__Desulfohalotomaculaceae,
# g__JBAIYS01)? Evidence so far -- composition compatible and nearest of 29;
# read taxonomy identical to MAG_022 contig 2833; no protein redundancy; depth
# 11.41 against that genome IQR 11.85-14.68, so just outside and 12% low, on a
# test where five unrelated genomes sit at 11-13x and which therefore has
# almost no power either way.
#
# A. WHICH GTDB FAMILY do Desulforamulus, Desulfofalx, Desulfosporosinus and
#    Desulfotomaculum belong to? NCBI still calls them all Peptococcaceae;
#    GTDB split that family, and MAG_022 is Desulfohalotomaculaceae. If the
#    genera contig 19299 places among ARE Desulfohalotomaculaceae in GTDB, the
#    phylogeny and the bin agree at family level. If they are a different GTDB
#    family, contig 19299 is NOT MAG_022 and the attribution dies here.
#    Answered offline from GTDB-Tk own reference taxonomy. No naming guesswork.
#
# B. Is dsrA on contig 19299 SPLIT, the way MAG_017 mcrA is? 58 showed proteins
#    _14 (213 aa) and _15 (183 aa) hitting THE SAME 400 aa reference over 211
#    and 163 residues. If those two alignments cover different halves of that
#    reference, dsrA here is one gene broken in two -- which also explains the
#    331.7 alpha score Volume XII recorded as anomalously low.
#
# C. MAG_022 own sulfur genes: which genera do THEY hit, and at what identity?
#    If MAG_022 aprA also lands on Desulforamulus/Desulfofalx at roughly the
#    same 78-83%, the binned genome and the unbinned operon sit at the same
#    evolutionary distance from the same relatives. That is corroboration from
#    genes that are not dsrAB at all.
#
# Read-only except 59_dsr_close/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/59_dsr_close"; mkdir -p "$OUT"
export TMPDIR="${TMPDIR:-$WS2/tmp}"; mkdir -p "$TMPDIR"
exec > >(tee "$OUT/59.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }

echo "59 -- dsrAB, closing the attribution"
echo "date : $(date -Is)"
echo "out  : $OUT"

PREV58="$WS2/58_dsr_finish"
PROT="$WS2/06_binqc/checkm2/protein_files"
CAND=LIQ_das_metabat_022
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CURL="curl -s --max-time 60"

# ================================================================ A
say "A. WHAT GTDB CALLS THE GENERA THE PHYLOGENY POINTED AT"

echo "  MAG_022 as GTDB-Tk classified it:"
awk -F'\t' -v b="$CAND" 'FNR>1 && $1==b{print "    "$2}' \
  "$WS2/08_taxonomy/gtdbtk/classify/gtdbtk.bac120.summary.tsv" 2>/dev/null

sub "locating GTDB-Tk reference taxonomy"
TAXF=""
for c in "${GTDBTK_DATA_PATH:-}/taxonomy/gtdb_taxonomy.tsv" \
         "$HOME/databases/gtdbtk"*/taxonomy/gtdb_taxonomy.tsv \
         "$HOME/gtdbtk"*/taxonomy/gtdb_taxonomy.tsv; do
  [ -s "$c" ] && { TAXF="$c"; break; }
done
if [ -z "$TAXF" ]; then
  TAXF="$(find "$HOME" ${SCRATCH_ROOT} -maxdepth 6 -type f -name 'gtdb_taxonomy.tsv' 2>/dev/null | head -1)"
fi
if [ -z "$TAXF" ]; then
  TAXF="$(find "$HOME" ${SCRATCH_ROOT} -maxdepth 6 -type f -name 'bac120_taxonomy*.tsv' 2>/dev/null | head -1)"
fi

if [ -n "$TAXF" ] && [ -s "$TAXF" ]; then
  ok "reference taxonomy: $TAXF  ($(wc -l < "$TAXF") genomes)"
  sub "GTDB family of each genus the DsrAB hits named"
  printf '  %-22s %-34s %s\n' "genus (NCBI name)" "GTDB family" "genomes"
  for G in Desulforamulus Desulfosporosinus Desulfotomaculum Desulfitobacterium \
           Pelotomaculum Desulfohalotomaculum; do
    awk -F'\t' -v g="g__$G;" '
      index($2, g) {
        n=split($2,a,";")
        for(i=1;i<=n;i++) if(a[i] ~ /^f__/) fam=a[i]
        c[fam]++ }
      END{ if(length(c)==0) printf "  %-22s %-34s %s\n", g, "(not present under this name)", "-"
           else for(f in c) printf "  %-22s %-34s %d\n", g, f, c[f] }' \
      OFS='\t' "$TAXF" | sed "s/g__$G;/$G/"
  done

  sub "and the alkaliphilic species that DsrA hit best (NCBI: Desulfofalx alkaliphila)"
  echo "  GTDB entries whose species epithet is alkaliphil*, within Desulfotomaculales:"
  grep -i 'o__Desulfotomaculales' "$TAXF" 2>/dev/null | grep -i 'alkaliphil' \
    | awk -F'\t' '{n=split($2,a,";"); printf "    %-16s %s ; %s\n", $1, a[5], a[7]}' | head -10
  [ -z "$(grep -i 'o__Desulfotomaculales' "$TAXF" 2>/dev/null | grep -i 'alkaliphil')" ] \
    && echo "    (none -- GTDB may place it under another epithet)"

  sub "what else lives in MAG_022 family"
  awk -F'\t' 'index($2,"f__Desulfohalotomaculaceae"){ n=split($2,a,";"); print a[6] }' "$TAXF" \
    | sort | uniq -c | sort -rn | head -12 | sed 's/^/    /'
else
  bad "GTDB-Tk reference taxonomy not found."
  echo "         Set GTDBTK_DATA_PATH, or find it with:"
  echo "           find / -name gtdb_taxonomy.tsv 2>/dev/null"
  echo "         Without it the family question cannot be answered offline."
fi

# ================================================================ B
say "B. IS dsrA ON CONTIG 19299 SPLIT IN TWO?"

if [ -s "$PREV58/refs.faa" ] && [ -s "$PREV58/queries.faa" ]; then
  makeblastdb -in "$PREV58/refs.faa" -dbtype prot -out "$OUT/refdb" >/dev/null 2>&1
  blastp -query "$PREV58/queries.faa" -db "$OUT/refdb" \
         -outfmt '6 qseqid sseqid pident length qstart qend sstart send slen evalue' \
         -max_target_seqs 8 -max_hsps 1 -num_threads 8 > "$OUT/blastp_coords.tsv" 2>/dev/null

  sub "where each query aligns ON the reference protein"
  awk -F'\t' '
    { split($2,s,"|")
      if($1!=q){ q=$1; printf "\n  %s\n", q
                 printf "      %-42s %6s %8s %14s %6s\n","reference","pid%","q range","s range","slen" }
      printf "      %-42s %6.1f %8s %14s %6d\n", s[2], $3, $5"-"$6, $7"-"$8, $9 }' "$OUT/blastp_coords.tsv"

  sub "the test: do _14 and _15 cover DIFFERENT halves of the same protein?"
  awk -F'\t' '
    $1 ~ /19299_1[45]/ {
      split($2,s,"|")
      key=s[3]
      if($1 ~ /_14/){ a14[key]=$7"-"$8; p14[key]=$3; l[key]=$9; n14[key]=s[2] }
      if($1 ~ /_15/){ a15[key]=$7"-"$8; p15[key]=$3; l[key]=$9; n15[key]=s[2] } }
    END{
      found=0
      for(k in a14) if(k in a15){
        found=1
        split(a14[k],x,"-"); split(a15[k],y,"-")
        s14=x[1]+0; e14=x[2]+0; s15=y[1]+0; e15=y[2]+0
        lo=(s14<s15?s14:s15); hi=(e14>e15?e14:e15)
        ov=(s14>s15?s14:s15); ov2=(e14<e15?e14:e15); overlap=ov2-ov+1; if(overlap<0) overlap=0
        printf "\n    reference %s (%s), %d aa\n", k, n14[k], l[k]
        printf "      _19299_15 covers %4d-%-4d  (%.1f%% id)\n", s15, e15, p15[k]
        printf "      _19299_14 covers %4d-%-4d  (%.1f%% id)\n", s14, e14, p14[k]
        printf "      union %d-%d = %d of %d aa (%.0f%%), overlap %d aa\n",
               lo, hi, hi-lo+1, l[k], 100*(hi-lo+1)/l[k], overlap
        if(overlap<=20 && (hi-lo+1) > 0.75*l[k])
          printf "      >>> TILES: two fragments of ONE dsrA gene\n"
        else if(overlap > 60)
          printf "      >>> OVERLAPPING: two copies of the same region, not a split gene\n"
        else
          printf "      >>> ambiguous on this reference\n" }
      if(!found) printf "\n    (no reference was hit by both _14 and _15)\n" }' "$OUT/blastp_coords.tsv"
else
  no "need 58 refs.faa and queries.faa -- run 58 first"
fi

# ================================================================ C
say "C. MAG_022 OWN SULFUR GENES -- same relatives, same distance?"

sub "finding MAG_022 aprA / aprB / sat / dsrMK proteins"
CP="$(find "$PROT" -name "${CAND}*.faa" | head -1)"
KTBL="$(find "$WS2/45_kofam_full" -name '*.tbl' 2>/dev/null | head -1)"
: > "$OUT/cand_ids.txt"
if [ -n "$KTBL" ]; then
  for T in "$WS2"/45_kofam_full/*.tbl "$WS2"/45_kofam_full/chunks/*.tbl; do
    [ -s "$T" ] || continue
    awk -v b="$CAND" '
      /^#/{next}
      { ko=$3
        if(ko!="K00394" && ko!="K00395" && ko!="K00958" && ko!="K27187" && ko!="K27188") next
        t=$1; p=index(t,"|"); if(p==0) next
        if(substr(t,1,p-1)!=b) next
        print substr(t,p+1) "\t" ko "\t" $6 }' "$T"
  done | sort -u > "$OUT/cand_ids.txt"
fi
NC=$(wc -l < "$OUT/cand_ids.txt")
if [ "$NC" -gt 0 ]; then
  ok "$NC candidate sulfur proteins in $CAND"
  awk -F'\t' '{printf "    %-26s %-8s score %s\n",$1,$2,$3}' "$OUT/cand_ids.txt"
  awk -F'\t' '{print $1}' "$OUT/cand_ids.txt" | sort -u > "$OUT/.want"
  awk 'BEGIN{while((getline l < ARGV[1])>0) w[l]=1; ARGV[1]=""}
       /^>/{ id=substr($1,2); p=(id in w); if(p) print ">CAND|" id; next }
       { if(p) print }' "$OUT/.want" "$CP" > "$OUT/cand.faa"
  ok "$(grep -c '^>' "$OUT/cand.faa" | head -1) sequences extracted"
else
  no "no aprA/aprB/sat/dsrMK hits found for $CAND in 45_kofam_full"
fi

if [ -s "$OUT/cand.faa" ] && $CURL -o /dev/null "$EUT/einfo.fcgi"; then
  fetch(){
    local out="$1" n="$2" term="$3" ids
    ids=$($CURL --get --data-urlencode "term=$term" "$EUT/esearch.fcgi?db=protein&retmax=$n" 2>/dev/null \
          | grep -oP '(?<=<Id>)[0-9]+(?=</Id>)' | tr '\n' ',' | sed 's/,$//')
    [ -z "$ids" ] && { printf '      %-56s 0\n' "$term"; return 1; }
    $CURL "$EUT/efetch.fcgi?db=protein&id=${ids}&rettype=fasta&retmode=text" >> "$out" 2>/dev/null
    printf '      %-56s %s\n' "$term" "$(printf '%s' "$ids" | tr ',' '\n' | wc -l)"
    sleep 0.4
  }
  sub "fetching AprA/AprB references for the same genera"
  if [ ! -s "$OUT/apr_raw.faa" ]; then
    : > "$OUT/apr_raw.faa"
    for SUB in alpha beta; do
      for ORG in Desulforamulus Desulfofalx Desulfosporosinus Desulfotomaculum \
                 Desulfitobacterium Pelotomaculum Desulfovibrio Archaeoglobus; do
        fetch "$OUT/apr_raw.faa" 4 "\"adenylylsulfate reductase subunit $SUB\"[Protein Name] AND \"$ORG\"[Organism]"
      done
    done
  fi
  if [ -s "$OUT/apr_raw.faa" ]; then
    awk '/^>/{ acc=substr($1,2); org=""
               if (match($0,/\[[^]]*\]$/)) org=substr($0,RSTART+1,RLENGTH-2)
               gsub(/[^A-Za-z0-9]+/,"_",org); if(org=="") org="unknown"
               if(acc in seen){skip=1;next} seen[acc]=1; skip=0
               print ">REF|" org "|" acc; next }
         { if(!skip) print }' "$OUT/apr_raw.faa" > "$OUT/apr_refs.faa"
    ok "$(grep -c '^>' "$OUT/apr_refs.faa" | head -1) Apr references"
    makeblastdb -in "$OUT/apr_refs.faa" -dbtype prot -out "$OUT/aprdb" >/dev/null 2>&1
    blastp -query "$OUT/cand.faa" -db "$OUT/aprdb" \
           -outfmt '6 qseqid sseqid pident length qlen slen evalue' \
           -max_target_seqs 5 -max_hsps 1 -num_threads 8 > "$OUT/blastp_apr.tsv" 2>/dev/null
    sub "MAG_022 sulfur proteins vs the same reference genera"
    awk -F'\t' '{ split($2,s,"|")
        if($1!=q){ q=$1; printf "\n  %s\n", q
                   printf "      %-42s %6s %7s %6s\n","nearest named relative","pid%","aln","qlen" }
        printf "      %-42s %6.1f %7d %6d\n", s[2], $3, $4, $5 }' "$OUT/blastp_apr.tsv"
    echo
    echo "  COMPARE WITH contig 19299 DsrB at 78.0% to Desulforamulus and DsrA at"
    echo "  82.9% to Desulfofalx. If MAG_022 own aprA/aprB sit at a similar"
    echo "  identity to the same genera, the binned genome and the unbinned"
    echo "  operon are the same evolutionary distance from the same relatives."
    echo "  If MAG_022 is much closer or much further, they are different"
    echo "  organisms and the attribution fails."
  fi
else
  no "skipped -- no candidate proteins or NCBI unreachable"
fi

say "D. THE SENTENCE THIS EARNS"
cat <<'EOT'
  ALREADY SAFE TO WRITE, whatever A and C say:

    Two distinct spore-forming sulfate reducers are encoded in the culture
    assembly and neither was binned. LIQ_contig_12880 carries a complete
    dsrAB pair 88-90% identical to Desulfosporosinus meridiei; LIQ_contig_19299
    carries a dsrAB pair nearest Desulforamulus and Desulfofalx at 78-83%.
    Both are REDUCTIVE-type dsr -- neither places among the reverse-type
    enzymes of sulfur oxidisers -- so they are sulfate reducers, and the
    competition premise of this chapter rests on measured phylogeny rather
    than on gene names.

  Volume XII called the spore-forming Bacillota reading a hypothesis with one
  gene behind it. It is now the answer, from two independent directions.

  THE ATTRIBUTION of contig 19299 to ETO_LIQ_MAG_022 stands or falls on section
  A. Same GTDB family -> write it as the most parsimonious reading, flagged as
  an attribution. Different family -> drop it, and report MAG_022 as a
  Desulfohalotomaculaceae that genuinely lacks dsrAB while two other sulfate
  reducers went unbinned. That second version is not a worse result; it is a
  more interesting one, because it means at least three sulfate-reducing
  lineages are present and the binning recovered the wrong parts of all three.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
