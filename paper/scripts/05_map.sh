#!/usr/bin/env bash
# =============================================================================
# 05_map.sh  --  read mapping (self + cross) and multi-sample contig depth
#
# PURPOSE
#   1. Coverage for binning. Critically, this produces TWO-SAMPLE differential
#      coverage: each assembly is mapped with its own reads AND the other
#      sample's reads. An organism present in sediment but absent from the
#      cultures (or vice versa) then has a distinctive coverage signature that
#      a single-sample depth profile cannot see. workspace_1's sediment
#      consensus took 5/5 bins from a single binner, which is what a
#      weak/under-informative depth profile looks like.
#   2. The headline QC number missing from workspace_1: what fraction of reads
#      map back to their own assembly, i.e. how much of the community was
#      actually assembled. Reviewers ask for this first.
#
# INPUTS   04_polish/{SED,LIQ}_asm_final.fasta   (polished, contigs prefixed)
#          00_input/{sediment,liquid}_reads.fastq.gz  (full unfiltered set)
# OUTPUTS  02_mapping/{SED,LIQ}_map_{self,cross}.sorted.bam (+ .bai)
#          02_mapping/{SED,LIQ}_depth_multi.txt   <- feeds MetaBAT2
#          02_mapping/mapping_summary.tsv
#
# RUNTIME  ~1-3 h. Safe to run in the foreground but tmux is still advisable.
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="05_map"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
M="$WS2/02_mapping"; mkdir -p "$M"

echo "######## 05 — MAPPING (self + cross) ########"
echo "  started: $(date -Is)  threads: $THREADS"
record_version "minimap2" "${ENV_minimap2:--}" minimap2 --version
record_version "samtools" "${ENV_samtools:--}" samtools --version
record_version "jgi_summarize_bam_contig_depths" "${ENV_jgi_summarize_bam_contig_depths:--}" jgi_summarize_bam_contig_depths --help
runlog "05_map started"

SED_ASM_P="$WS2/04_polish/SED_asm_final.fasta"
LIQ_ASM_P="$WS2/04_polish/LIQ_asm_final.fasta"
for f in "$SED_ASM_P" "$LIQ_ASM_P"; do
  [ -s "$f" ] || { echo "  !! missing $f — run 04_polish.sh first"; exit 1; }
done

# ------------------------------------------------------------------- MAP
map() {  # map <label> <ref.fasta> <reads> <outprefix>
  local lab="$1" ref="$2" reads="$3" out="$4"
  if [ -s "${out}.sorted.bam.bai" ]; then echo "  $lab: already mapped, skipping"; return 0; fi
  echo
  echo "--- $lab"
  echo "    ref   : $(basename "$ref")"
  echo "    reads : $(basename "$reads")"
  echo "    start : $(date -Is)"
  run_tool minimap2 -ax map-ont -t "$THREADS" --secondary=no "$ref" "$reads" \
    | run_tool samtools sort -@ 6 -m 2G -o "${out}.sorted.bam" -
  run_tool samtools index -@ 6 "${out}.sorted.bam"
  echo "    done  : $(date -Is)"
}

map "SED reads -> SED assembly (self)"  "$SED_ASM_P" "$SED_READS" "$M/SED_map_self"
map "LIQ reads -> SED assembly (cross)" "$SED_ASM_P" "$LIQ_READS" "$M/SED_map_cross"
map "LIQ reads -> LIQ assembly (self)"  "$LIQ_ASM_P" "$LIQ_READS" "$M/LIQ_map_self"
map "SED reads -> LIQ assembly (cross)" "$LIQ_ASM_P" "$SED_READS" "$M/LIQ_map_cross"

# ------------------------------------------------------- % READS MAPPED
echo
echo "######## % OF READS RECRUITED TO EACH ASSEMBLY ########"
printf 'bam\ttotal_reads\tmapped_reads\tpct_mapped\n' > "$M/mapping_summary.tsv"
for b in SED_map_self SED_map_cross LIQ_map_self LIQ_map_cross; do
  [ -s "$M/$b.sorted.bam" ] || continue
  fs=$(run_tool samtools flagstat -@ 6 "$M/$b.sorted.bam")
  echo "$fs" > "$M/${b}.flagstat.txt"
  tot=$(echo "$fs" | awk 'NR==1{print $1}')
  prim=$(echo "$fs" | grep 'primary mapped' | awk '{print $1}')
  [ -z "$prim" ] && prim=$(echo "$fs" | grep -m1 ' mapped (' | awk '{print $1}')
  pct=$(awk -v a="$prim" -v b="$tot" 'BEGIN{printf "%.1f", (b>0? 100*a/b : 0)}')
  printf '%s\t%s\t%s\t%s\n' "$b" "$tot" "$prim" "$pct" >> "$M/mapping_summary.tsv"
  printf '  %-16s %12s reads | %12s mapped | %5s%%\n' "$b" "$tot" "$prim" "$pct"
done
echo
echo "  The SELF rows are the number to report: the fraction of the metagenome"
echo "  that assembled. The CROSS rows quantify overlap between the two samples."

# ------------------------------------------------- MULTI-SAMPLE DEPTH
echo
echo "######## CONTIG DEPTH (both BAMs -> differential coverage) ########"
for s in SED LIQ; do
  out="$M/${s}_depth_multi.txt"
  if [ -s "$out" ]; then echo "  $s: depth already computed"; continue; fi
  echo "  computing $s depth from self + cross ..."
  run_tool jgi_summarize_bam_contig_depths \
    --outputDepth "$out" \
    "$M/${s}_map_self.sorted.bam" "$M/${s}_map_cross.sorted.bam"
  printf '    %s: %s contigs, %s columns\n' "$s" "$(( $(wc -l < "$out") - 1 ))" \
    "$(head -1 "$out" | awk '{print NF}')"
done

# ------------------------------------------------------------ VALIDATION
echo
echo "######## VALIDATION GATE ########"
for s in SED LIQ; do
  check "$s self BAM indexed"  $([ -s "$M/${s}_map_self.sorted.bam.bai" ] && echo 0 || echo 1)
  check "$s cross BAM indexed" $([ -s "$M/${s}_map_cross.sorted.bam.bai" ] && echo 0 || echo 1)
  check "$s multi-sample depth built" $([ -s "$M/${s}_depth_multi.txt" ] && echo 0 || echo 1)
  # depth table must have 2 sample columns (contigName len totalAvgDepth + 2x(depth,var))
  if [ -s "$M/${s}_depth_multi.txt" ]; then
    nf=$(head -1 "$M/${s}_depth_multi.txt" | awk '{print NF}')
    check "$s depth has both samples (${nf} cols, expect 7)" $([ "$nf" -ge 7 ] && echo 0 || echo 1)
  fi
done
for s in SED LIQ; do
  p=$(awk -F'\t' -v k="${s}_map_self" '$1==k{print $4}' "$M/mapping_summary.tsv")
  [ -n "$p" ] && { awk -v p="$p" 'BEGIN{exit !(p>=70)}' \
      && check "$s self-mapping >=70% (${p}%)" 0 \
      || check "$s self-mapping >=70% (${p}%) — report honestly, do not hide" 1; }
done

runlog "05_map complete"
echo
echo "######## 05 DONE — next: scripts/06_bin.sh ########"
