#!/usr/bin/env bash
# 55_dsr_identity.sh -- who do the two unbinned dsrAB contigs belong to?
#
#   bash scripts/run.sh scripts/55_dsr_identity_XXXXXXX.sh
#
# THE QUESTION. LIQ_contig_12880 (6,395 bp) and LIQ_contig_19299 (20,651 bp)
# each carry a dsrA/dsrB pair, far above curated thresholds. Neither belongs to
# any MAG. Neither carries a ribosomal protein, so gene content cannot place
# them, and GC does not match JBAIYS01. Volume XII recorded the source organism
# as UNDETERMINED with one lead: a spore germination gene (K06313) on the long
# contig, pointing at a spore-forming Bacillota sulfate reducer.
#
# THE ROUTE. Volume XII also recorded that 35_kaiju is keyed on READ ids, not
# contig names, which is why grepping it for genus names returned nothing. That
# is a reason to come at it from the other end:
#
#     reads that map to these two contigs  ->  their ids  ->  Kaiju taxonomy
#
# Kaiju classifies reads. If we know which reads built a contig, we know what
# Kaiju called them. No reference download, no network, no new database.
#
# THE CONTROL, AND WHY IT DECIDES EVERYTHING. The same procedure is run on
# contigs from a MAG whose taxonomy is already known from GTDB-Tk. If the
# control contig does not come back as its known taxon, the method is broken and
# the dsr answer is discarded unread. Rule 95: a search that returns nothing --
# or anything -- for a known positive is not a result.
#
# Read-only except 55_dsr_id/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/55_dsr_id"; mkdir -p "$OUT"
export TMPDIR="${TMPDIR:-$WS2/tmp}"; mkdir -p "$TMPDIR"
exec > >(tee "$OUT/55.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }

echo "55 -- dsrAB contig identity"
echo "date : $(date -Is)"
echo "out  : $OUT"
echo "tmp  : $TMPDIR"

DSR_CONTIGS="LIQ_contig_12880 LIQ_contig_19299"

# ================================================================ A
say "A. THE QUERY CONTIGS"

ASM="$(find "$WS2/04_polish" -maxdepth 2 -type f -name 'LIQ_asm*.fasta' 2>/dev/null | head -1)"
[ -z "$ASM" ] && ASM="$(find "$WS2" -maxdepth 3 -type f -name 'LIQ_asm_final.fasta' 2>/dev/null | head -1)"
if [ -n "$ASM" ]; then
  ok "assembly: ${ASM#$WS2/}"
else
  no "LIQ assembly fasta not found -- contig sequences will not be extracted"
fi

if [ -n "$ASM" ] && have samtools; then
  [ -f "$ASM.fai" ] || samtools faidx "$ASM" 2>/dev/null || true
  if [ -f "$ASM.fai" ]; then
    for C in $DSR_CONTIGS; do
      awk -v c="$C" '$1==c{printf "  %-20s %10d bp\n", $1, $2}' "$ASM.fai"
    done
  fi
fi

sub "the dsr proteins recorded in step 49"
DSRFAA="$(find "$WS2/49_dsr" -type f -name 'dsr_contig_proteins.faa' 2>/dev/null | head -1)"
if [ -n "$DSRFAA" ]; then
  ok "${DSRFAA#$WS2/}  ($(grep -c '^>' "$DSRFAA") proteins)"
  grep '^>' "$DSRFAA" | sed 's/^/    /' | head -30
else
  no "49_dsr/dsr_contig_proteins.faa not found"
fi

# ================================================================ B
say "B. IS THERE A BAM OF READS AGAINST THE LIQ ASSEMBLY?"

BAM=""
while IFS= read -r b; do
  [ -s "$b" ] || continue
  if samtools view -H "$b" 2>/dev/null | grep -q "SN:LIQ_contig_12880"; then
    BAM="$b"; break
  fi
done < <(find "$WS2" -type f \( -name '*.bam' -o -name '*.cram' \) 2>/dev/null | sort)

if [ -n "$BAM" ]; then
  ok "found: ${BAM#$WS2/}"
  echo "         reference contains LIQ_contig_12880, so read names are recoverable"
  [ -f "$BAM.bai" ] || [ -f "${BAM%.bam}.bai" ] || { echo "         indexing..."; samtools index -@ 4 "$BAM" 2>/dev/null || true; }
else
  no "no existing BAM is mapped to the LIQ assembly"
fi

# fallback: map reads ourselves. Target is 27 kb, so this is quick.
if [ -z "$BAM" ]; then
  sub "fallback -- map the culture reads to the two contigs directly"
  READS="$(find "$WS2/00_input" -maxdepth 2 -type f \( -name '*LIQ*.fastq.gz' -o -name '*liq*.fastq.gz' -o -name '*culture*.fastq.gz' \) 2>/dev/null | head -1)"
  if [ -z "$READS" ] || ! have minimap2; then
    bad "cannot build one: reads=${READS:-NOT FOUND}, minimap2=$(have minimap2 && echo yes || echo NO)"
    echo "         Sections D and E cannot run. C, F and G still will."
  elif [ -z "$ASM" ]; then
    bad "no assembly to extract the contigs from"
  else
    ok "reads: ${READS#$WS2/}"
    samtools faidx "$ASM" $DSR_CONTIGS > "$OUT/dsr_contigs.fasta" 2>/dev/null
    if [ -s "$OUT/dsr_contigs.fasta" ]; then
      ok "extracted $(grep -c '^>' "$OUT/dsr_contigs.fasta") contigs"
      echo "  running: minimap2 -x map-ont -a --secondary=no -t 8"
      minimap2 -x map-ont -a --secondary=no -t 8 "$OUT/dsr_contigs.fasta" "$READS" 2>"$OUT/minimap2.err" \
        | samtools sort -@ 4 -o "$OUT/dsr_reads.bam" - 2>/dev/null
      samtools index "$OUT/dsr_reads.bam" 2>/dev/null
      if [ -s "$OUT/dsr_reads.bam" ]; then
        BAM="$OUT/dsr_reads.bam"; ok "built ${BAM#$WS2/}"
      else
        bad "minimap2/samtools produced nothing -- see minimap2.err"
      fi
    else
      bad "could not extract the contigs from the assembly"
    fi
  fi
fi

# ================================================================ C
say "C. WHAT IS ACTUALLY IN THE KAIJU FILES?"

mapfile -t KJ < <(find "$WS2/35_kaiju" -type f \( -name '*.kaiju' -o -name '*.out' -o -name '*.tsv' \) 2>/dev/null | sort)
if [ "${#KJ[@]}" -eq 0 ]; then
  mapfile -t KJ < <(find "$WS2" -type f -name '*.kaiju' 2>/dev/null | sort)
fi
if [ "${#KJ[@]}" -eq 0 ]; then
  bad "no kaiju output found"
else
  ok "${#KJ[@]} kaiju file(s)"
  for f in "${KJ[@]}"; do
    n=$(wc -l < "$f")
    c=$(head -1 "$f" | awk -F'\t' '{print NF}')
    printf '  %-44s %10d lines, %2d cols\n' "${f#$WS2/}" "$n" "$c"
    head -2 "$f" | sed 's/^/      /' | cut -c1-160
  done
fi

sub "taxonomy dump for taxid -> name"
NAMES="$(find "$WS2" "$HOME" -maxdepth 6 -type f -name 'names.dmp' 2>/dev/null | head -1)"
if [ -n "$NAMES" ]; then ok "names.dmp: $NAMES"; else no "names.dmp not found -- taxids will be reported as numbers"; fi
have kaiju-addTaxonNames && ok "kaiju-addTaxonNames is on PATH" || no "kaiju-addTaxonNames not on PATH"

# ================================================================ D
say "D. THE POSITIVE CONTROL -- a contig from a genome whose taxonomy is known"

GTDB="$(find "$WS2" -type f -name 'gtdbtk.bac120.summary.tsv' 2>/dev/null | head -1)"
CTRL_BIN=""; CTRL_TAX=""; CTRL_CONTIG=""
if [ -n "$GTDB" ]; then
  ok "GTDB-Tk summary: ${GTDB#$WS2/}"
  # pick the bacterial MAG with the longest classification string as control
  CTRL_BIN="$(awk -F'\t' 'NR>1 && $1 ~ /^LIQ_/ {print $1"\t"$2}' "$GTDB" | head -1 | cut -f1)"
  CTRL_TAX="$(awk -F'\t' -v b="$CTRL_BIN" 'NR>1 && $1==b{print $2}' "$GTDB")"
  [ -n "$CTRL_BIN" ] && echo "  control genome: $CTRL_BIN"
  [ -n "$CTRL_TAX" ] && echo "  known taxonomy: $CTRL_TAX"
else
  no "no GTDB-Tk summary found -- the control cannot name its expected answer"
fi

BINFA="$(find "$WS2" -type f -name "${CTRL_BIN}.fa*" 2>/dev/null | head -1)"
if [ -n "$CTRL_BIN" ] && [ -n "$BINFA" ]; then
  ok "control bin fasta: ${BINFA#$WS2/}"
  CTRL_CONTIG="$(grep '^>' "$BINFA" | head -1 | sed 's/^>//' | awk '{print $1}' | sed 's/.*~//')"
  echo "  control contig: $CTRL_CONTIG"
else
  no "control bin fasta not found -- control will be skipped"
fi

# ================================================================ E
say "E. READS ON EACH CONTIG, AND WHAT KAIJU CALLED THEM"

lookup(){   # $1 = contig, $2 = label
  local C="$1" L="$2"
  local IDS="$OUT/.ids_$C"
  samtools view "$BAM" "$C" 2>/dev/null | awk '{print $1}' | sort -u > "$IDS"
  local N; N=$(wc -l < "$IDS")
  printf '\n  %-22s %-18s reads mapped: %d\n' "$L" "$C" "$N"
  if [ "$N" -eq 0 ]; then echo "      (no reads -- nothing to look up)"; return; fi

  : > "$OUT/.tax_$C"
  for f in "${KJ[@]}"; do
    awk -F'\t' -v ids="$IDS" '
      BEGIN{ while((getline l < ids)>0) keep[l]=1 }
      { r=$2; sub(/[ \t].*$/,"",r)
        if(r in keep) print $0 }' "$f" >> "$OUT/.tax_$C"
  done

  local M; M=$(wc -l < "$OUT/.tax_$C")
  echo "      classified rows recovered from kaiju: $M"
  if [ "$M" -eq 0 ]; then
    echo "      >>> zero. Either these reads are not in the kaiju run, or the id"
    echo "          formats differ. Compare the two samples printed below."
    echo "      bam read id   : $(head -1 "$IDS")"
    echo "      kaiju read id : $(awk -F'\t' 'NR==1{print $2}' "${KJ[0]}" 2>/dev/null)"
    return
  fi

  # classified fraction, then the taxon distribution
  awk -F'\t' '{c[$1]++} END{ for(k in c) printf "      status %s : %d\n", k, c[k] }' "$OUT/.tax_$C"

  if [ -n "$NAMES" ]; then
    awk -F'\t' '$1=="C"{print $3}' "$OUT/.tax_$C" | sort | uniq -c | sort -rn | head -12 \
    | while read -r cnt tid; do
        nm=$(awk -F'\t\\|\t' -v t="$tid" '$1==t && $4 ~ /scientific name/{print $2; exit}' "$NAMES")
        printf '      %6d  %-10s %s\n' "$cnt" "$tid" "${nm:-<taxid not in names.dmp>}"
      done
  else
    awk -F'\t' '$1=="C"{print $3}' "$OUT/.tax_$C" | sort | uniq -c | sort -rn | head -12 \
    | awk '{printf "      %6d  taxid %s\n", $1, $2}'
  fi
}

if [ -z "$BAM" ] || [ "${#KJ[@]}" -eq 0 ]; then
  bad "cannot run: BAM=${BAM:-none}, kaiju files=${#KJ[@]}"
else
  if [ -n "$CTRL_CONTIG" ]; then
    sub "POSITIVE CONTROL -- expected: $CTRL_TAX"
    lookup "$CTRL_CONTIG" "CONTROL"
  else
    no "no control contig; the dsr result below is UNCONTROLLED and must not be quoted"
  fi

  sub "THE QUESTION"
  for C in $DSR_CONTIGS; do lookup "$C" "dsr contig"; done
fi

# ================================================================ F
say "F. rRNA ON THE CONTIGS -- the direct route, if anything is there"

if have barrnap; then
  ok "barrnap present"
  if [ -s "$OUT/dsr_contigs.fasta" ]; then
    for K in bac arc; do
      echo "  --- kingdom $K"
      barrnap --kingdom "$K" --quiet "$OUT/dsr_contigs.fasta" 2>/dev/null \
        | awk '!/^#/{printf "    %-18s %-10s %8s..%-8s %s\n", $1, $3, $4, $5, $9}'
    done
  else
    no "contig fasta was not extracted in section B"
  fi
else
  no "barrnap not on PATH. Volume XII reported no ribosomal PROTEINS; rRNA genes"
  echo "         were never checked, and a 20 kb contig can carry a 16S."
  echo "         If barrnap exists in another env, this is worth one command."
fi

# ================================================================ G
say "G. FALLBACK ROUTE -- what a reference tree would need"

sub "how step 15 obtained its reference sequences"
S15="$(find "$WS2" -maxdepth 2 -type f -name '15*marker*tree*.sh' 2>/dev/null | head -1)"
if [ -n "$S15" ]; then
  ok "${S15#$WS2/}"
  grep -n -iE 'wget|curl|uniprot|kegg|refseq|download|REF_|reference' "$S15" | head -25 | sed 's/^/    /'
else
  no "step 15 script not found at the top level"
fi

sub "alignment and tree tools"
for t in mafft muscle trimal FastTree fasttree iqtree iqtree2 raxml-ng blastp diamond hmmsearch seqkit; do
  if have "$t"; then printf '  [OK]   %-10s %s\n' "$t" "$(command -v $t)"; else printf '  [--]   %-10s not on PATH\n' "$t"; fi
done

sub "local sequence databases"
find "$WS2" "$HOME" -maxdepth 5 -type f \( -name '*.dmnd' -o -name '*.pin' -o -name '*.phr' \) 2>/dev/null \
  | head -15 | sed 's/^/    /' || true

rm -f "$OUT"/.ids_* 2>/dev/null

say "H. HOW TO READ THIS"
cat <<'EOT'
  Section E is the result, and section E is only readable if its CONTROL block
  returned the taxon that GTDB-Tk already assigned to that genome.

  CONTROL CORRECT, dsr contigs return a consistent taxon
      -> that is the answer. Report it as read-level taxonomy of the reads that
         assembled the contig, with the classified fraction, not as a genome
         assignment. Expect a spread: Kaiju classifies individual long reads and
         disagrees with itself at species level. Genus or family agreement
         across most reads is what counts.

  CONTROL CORRECT, dsr contigs return a scatter with no majority
      -> the honest answer stays "undetermined", but it is now undetermined
         WITH A MEASUREMENT behind it instead of an unexplored option.

  CONTROL WRONG OR EMPTY
      -> discard section E entirely. Do not read the dsr numbers. Fix the
         method first: the likely cause is a read-id mismatch between the BAM
         and the kaiju run, which the printed samples will show.

  Section F is the route that would settle it outright. A 16S on
  LIQ_contig_19299 identifies the organism directly, and Volume XII only ever
  ruled out ribosomal PROTEINS.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
