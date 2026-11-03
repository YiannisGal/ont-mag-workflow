#!/usr/bin/env bash
# =============================================================================
# 04_polish.sh  --  Medaka consensus polishing + contig renaming
#
# PURPOSE — this is THE step workspace_2 exists for.
#   workspace_1 MAGs had a median average gene length of 244 aa against an
#   expected 290-320 aa, including in a single closed 4.90 Mb contig where
#   fragmentation cannot explain it. Residual homopolymer/indel errors were
#   frameshifting ORFs, splitting genes into fragments, and every functional
#   result downstream (DRAM pathways, CAZymes, mcrA/dsrA, proteases) inherited
#   the undercount.
#
#   Flye's internal polishing is not sufficient at metagenomic coverage.
#   Medaka re-calls consensus with a neural network trained on the MATCHING
#   basecaller model, which is the only step that reliably clears these errors.
#   It was applied to the culturomics isolates but never to the metagenome.
#
# READS USED: the FULL unfiltered set, not the >=1000 bp assembly set.
#   Assembly needed long reads (short ones cannot form overlaps). Consensus
#   needs COVERAGE, and sub-kilobase reads align perfectly well and contribute.
#   With a sediment median read length of 364 bp, excluding them would discard
#   most of the available consensus support.
#
# ALSO DOES: contig renaming per NAMING_CONVENTION.md.
#   Flye numbers contigs independently per sample, so contig_1 exists in BOTH
#   assemblies and means two different organisms. Every cross-sample step
#   (cross-mapping, ANI, merged MAG sets, shared trees) would silently mix
#   them. Prefixing here, before any mapping, makes every contig ID
#   self-identifying for the rest of the project.
#
# RUNTIME  6-15 h for both samples. RUN IN tmux, LEAVE OVERNIGHT.
#
# USAGE
#   cd ${WS2}
#   bash scripts/04_polish.sh 2>&1 | tee logs/04_polish_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="04_polish"
export TMPDIR="${TMPDIR}"   # short path: AF_UNIX limit
mkdir -p "$TMPDIR" "$WS2/04_polish"

echo "######## 04 — MEDAKA POLISHING + CONTIG RENAMING ########"
echo "  started: $(date -Is)   threads: $THREADS"
echo "  model  : $MEDAKA_MODEL"
record_version "medaka" "${ENV_medaka:-medaka}" medaka --version
record_version "seqkit" "${ENV_seqkit:--}"      seqkit version
runlog "04_polish started (model=$MEDAKA_MODEL, full unfiltered reads)"

# ------------------------------------------------------------------ POLISH
polish() {  # polish <CODE> <assembly.fasta> <reads.fastq.gz>
  local code="$1" asm="$2" reads="$3"
  local out="$WS2/04_polish/${code}_medaka"

  echo
  echo "================================================================"
  echo "  POLISHING: $code"
  echo "  assembly: $asm"
  echo "  reads   : $reads  (full set — coverage matters here)"
  echo "  start   : $(date -Is)"
  echo "================================================================"

  if [ -s "$out/consensus.fasta" ]; then
    echo "  consensus.fasta already present — skipping (delete to force)"
    return 0
  fi

  run_tool medaka_consensus \
      -i "$reads" \
      -d "$asm" \
      -o "$out" \
      -t "$THREADS" \
      -m "$MEDAKA_MODEL" \
      -f
  local rc=$?
  echo "  finished: $(date -Is)  (exit $rc)"
  return $rc
}

polish "SED" "$WS2/03_assembly/sediment/assembly.fasta"        "$SED_READS"; sed_rc=$?
polish "LIQ" "$WS2/03_assembly/liquid_cultures/assembly.fasta" "$LIQ_READS"; liq_rc=$?

# ------------------------------------------------- RENAME CONTIGS (naming spec)
echo
echo "######## CONTIG RENAMING (NAMING_CONVENTION.md §4) ########"
rename_contigs() {  # rename_contigs <CODE>
  local code="$1"
  local src="$WS2/04_polish/${code}_medaka/consensus.fasta"
  local dst="$WS2/04_polish/${code}_asm_final.fasta"
  local map="$WS2/04_polish/${code}_contig_rename.tsv"
  [ -s "$src" ] || { echo "  $code: no consensus.fasta, skipping rename"; return 1; }
  if [ -s "$dst" ]; then echo "  $code: already renamed"; return 0; fi

  sed "s/^>/>${code}_/" "$src" > "$dst"
  # auditable old -> new mapping
  printf 'old_id\tnew_id\n' > "$map"
  grep '^>' "$src" | sed 's/^>//' | awk -v c="$code" '{print $1"\t"c"_"$1}' >> "$map"

  printf '  %s: %s contigs renamed  ->  %s_contig_*\n' "$code" \
    "$(( $(wc -l < "$map") - 1 ))" "$code"
  echo "     $dst"
  echo "     mapping: $map"
}
rename_contigs "SED"
rename_contigs "LIQ"

# ------------------------------------------------------- BEFORE / AFTER STATS
echo
echo "######## ASSEMBLY STATS: UNPOLISHED vs POLISHED ########"
stats() {
  [ -s "$1" ] || { echo "MISSING"; return; }
  awk '/^>/{if(l)print l; l=0; next}{l+=length($0)}END{if(l)print l}' "$1" \
  | sort -rn | awk '{a[NR]=$1; s+=$1} END{
      h=s/2; c=0; for(i=1;i<=NR;i++){c+=a[i]; if(c>=h){n50=a[i]; break}}
      printf "%d contigs | %.1f Mbp | N50 %.1f kb | longest %.2f Mb", NR, s/1e6, n50/1000, a[1]/1e6}'
}
printf '\n  %-22s %s\n' "SED unpolished:" "$(stats "$WS2/03_assembly/sediment/assembly.fasta")"
printf '  %-22s %s\n'   "SED polished:"   "$(stats "$WS2/04_polish/SED_asm_final.fasta")"
printf '\n  %-22s %s\n' "LIQ unpolished:" "$(stats "$WS2/03_assembly/liquid_cultures/assembly.fasta")"
printf '  %-22s %s\n'   "LIQ polished:"   "$(stats "$WS2/04_polish/LIQ_asm_final.fasta")"
echo
echo "  (length should barely change — Medaka corrects bases, it does not"
echo "   restructure the assembly. A large change means something went wrong.)"

# ------------------------------------------------------------- VALIDATION
echo
echo "######## VALIDATION GATE ########"
check "SED consensus produced" $([ -s "$WS2/04_polish/SED_medaka/consensus.fasta" ] && echo 0 || echo 1)
check "LIQ consensus produced" $([ -s "$WS2/04_polish/LIQ_medaka/consensus.fasta" ] && echo 0 || echo 1)
check "SED renamed + final"    $([ -s "$WS2/04_polish/SED_asm_final.fasta" ] && echo 0 || echo 1)
check "LIQ renamed + final"    $([ -s "$WS2/04_polish/LIQ_asm_final.fasta" ] && echo 0 || echo 1)
check "medaka exited cleanly (SED)" "$sed_rc"
check "medaka exited cleanly (LIQ)" "$liq_rc"

# contig count must be preserved by polishing
for c in SED LIQ; do
  src="$WS2/03_assembly/$([ "$c" = SED ] && echo sediment || echo liquid_cultures)/assembly.fasta"
  dst="$WS2/04_polish/${c}_asm_final.fasta"
  if [ -s "$src" ] && [ -s "$dst" ]; then
    a=$(grep -c '^>' "$src"); b=$(grep -c '^>' "$dst")
    check "$c contig count preserved ($a -> $b)" $([ "$a" -eq "$b" ] && echo 0 || echo 1)
  fi
done

# every contig must now carry its sample prefix
for c in SED LIQ; do
  f="$WS2/04_polish/${c}_asm_final.fasta"
  [ -s "$f" ] || continue
  bad=$(grep '^>' "$f" | grep -cv "^>${c}_")
  check "$c all contigs prefixed (${bad} unprefixed)" $([ "$bad" -eq 0 ] && echo 0 || echo 1)
done

runlog "04_polish complete (sed_rc=$sed_rc liq_rc=$liq_rc)"
echo
echo "######## 04 DONE ########"
echo "Outputs for the next step:"
echo "  $WS2/04_polish/SED_asm_final.fasta"
echo "  $WS2/04_polish/LIQ_asm_final.fasta"
echo
echo "Next: 05_map.sh (self + cross mapping), then 06_bin.sh (4 binners)."
