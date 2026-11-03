#!/usr/bin/env bash
# =============================================================================
# 06_bin.sh  --  three independent binners, named per NAMING_CONVENTION.md
#
# PURPOSE
#   Give DAS_Tool genuinely independent inputs. workspace_1's sediment
#   consensus took 5 of 5 bins from a single binner, which is the signature of
#   an under-diverse input set. Three algorithmically distinct binners, all fed
#   TWO-SAMPLE differential coverage from step 05:
#
#     MetaBAT2  — tetranucleotide frequency + coverage, deterministic
#     VAMB      — variational autoencoder over composition + abundance
#     SemiBin2  — self-supervised deep learning; consistently the strongest on
#                 long-read, low-coverage assemblies, and already in use in
#                 your culturomics work
#
#   TaxVAMB is deliberately NOT run here: it requires a per-contig taxonomy
#   annotation (an mmseqs2/metabuli pre-step measured in hours) and SemiBin2
#   covers the same "informed binning" niche. If sediment recovery is poor
#   after DAS_Tool, adding TaxVAMB is the first escalation — see 06b notes.
#
# NAMING (NAMING_CONVENTION.md §5)
#   Every bin is renamed on output to  <SAMPLE>_<binner>_<NNN>.fa
#     SED_metabat_007.fa   LIQ_semibin_003.fa   SED_vamb_012.fa
#   Directory and filename agree; neither alone is trusted.
#
# INPUTS   04_polish/{SED,LIQ}_asm_final.fasta
#          02_mapping/{SED,LIQ}_depth_multi.txt
#          02_mapping/{SED,LIQ}_map_{self,cross}.sorted.bam
# OUTPUTS  03_binning/{metabat,vamb,semibin}/  + bin_counts.tsv
#
# RUNTIME  ~2-5 h. VAMB/SemiBin2 will use the GPU if torch sees it.
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="06_bin"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
B="$WS2/03_binning"; M="$WS2/02_mapping"
mkdir -p "$B"/{metabat,vamb,semibin}

echo "######## 06 — BINNING (MetaBAT2 / VAMB / SemiBin2) ########"
echo "  started: $(date -Is)  threads: $THREADS"
record_version "metabat2" "${ENV_metabat2:--}"   metabat2 --help
record_version "vamb"     "${ENV_vamb:-vamb_env}" vamb --version
record_version "SemiBin2" "${ENV_SemiBin2:-SemiBin}" SemiBin2 --version
runlog "06_bin started (metabat2 + vamb + semibin2, two-sample depth)"

MINCONTIG=2000     # long-read assemblies tolerate a higher floor than 1500
MINBIN=200000      # 200 kb minimum bin size

# ---------------------------------------------------------------- helpers
# rename a binner's output into <SAMPLE>_<binner>_<NNN>.fa
collect() {  # collect <sample> <binner> <srcdir> <glob>
  local s="$1" bn="$2" src="$3" pat="$4"
  local dst="$B/$bn"; mkdir -p "$dst"
  local i=0
  shopt -s nullglob
  for f in $src/$pat; do
    [ -s "$f" ] || continue
    # keep the binner's own id where it is numeric, else sequential
    local id; id=$(basename "$f" | grep -oE '[0-9]+' | head -1)
    [ -z "$id" ] && { i=$((i+1)); id=$i; }
    printf -v pad '%03d' "$id"
    cp "$f" "$dst/${s}_${bn}_${pad}.fa"
  done
  shopt -u nullglob
  local n; n=$(ls -1 "$dst"/${s}_${bn}_*.fa 2>/dev/null | wc -l)
  printf '    -> %s bins as %s_%s_NNN.fa\n' "$n" "$s" "$bn"
}

run_sample() {  # run_sample <SED|LIQ> <assembly>
  local s="$1" asm="$2"
  local depth="$M/${s}_depth_multi.txt"
  local bam_self="$M/${s}_map_self.sorted.bam" bam_cross="$M/${s}_map_cross.sorted.bam"
  echo
  echo "================================================================"
  echo "  SAMPLE $s"
  echo "================================================================"

  # ---- MetaBAT2 (two-sample depth) ----
  echo "--- MetaBAT2  ($(date +%H:%M))"
  if ls "$B/metabat/${s}_metabat_"*.fa >/dev/null 2>&1; then echo "    already done"; else
    local o="$B/metabat/raw_$s"; mkdir -p "$o"
    run_tool metabat2 -i "$asm" -a "$depth" -o "$o/bin" \
        -m "$MINCONTIG" -s "$MINBIN" -t "$THREADS" --seed 1 --unbinned
    collect "$s" metabat "$o" "bin.*.fa"
  fi

  # ---- VAMB ----
  echo "--- VAMB  ($(date +%H:%M))"
  if ls "$B/vamb/${s}_vamb_"*.fa >/dev/null 2>&1; then echo "    already done"; else
    local o="$B/vamb/raw_$s"; rm -rf "$o"
    run_tool vamb bin default --outdir "$o" --fasta "$asm" \
        --bamfiles "$bam_self" "$bam_cross" \
        -m "$MINCONTIG" --minfasta "$MINBIN" -p "$THREADS" \
      || echo "    !! vamb failed — see log; continuing with the other binners"
    collect "$s" vamb "$o/bins" "*.fna"
  fi

  # ---- SemiBin2 (self-supervised; no pretrained model assumed) ----
  echo "--- SemiBin2  ($(date +%H:%M))"
  if ls "$B/semibin/${s}_semibin_"*.fa >/dev/null 2>&1; then echo "    already done"; else
    local o="$B/semibin/raw_$s"; rm -rf "$o"
    run_tool SemiBin2 single_easy_bin \
        -i "$asm" -b "$bam_self" "$bam_cross" -o "$o" \
        --self-supervised --threads "$THREADS" \
        --minfasta-kbs $((MINBIN/1000)) --sequencing-type=long_read \
      || echo "    !! SemiBin2 failed — see log; continuing"
    collect "$s" semibin "$o/output_bins" "*.fa*"
  fi
}

run_sample SED "$WS2/04_polish/SED_asm_final.fasta"
run_sample LIQ "$WS2/04_polish/LIQ_asm_final.fasta"

# ------------------------------------------------------------ BIN COUNTS
echo
echo "######## BINS PER BINNER PER SAMPLE ########"
printf 'sample\tbinner\tn_bins\ttotal_Mbp\n' > "$B/bin_counts.tsv"
printf '  %-8s %-10s %8s %12s\n' sample binner bins total_Mbp
for s in SED LIQ; do for bn in metabat vamb semibin; do
  n=$(ls -1 "$B/$bn/${s}_${bn}_"*.fa 2>/dev/null | wc -l)
  mb=$(cat "$B/$bn/${s}_${bn}_"*.fa 2>/dev/null | grep -v '^>' | tr -d '\n' | wc -c \
       | awk '{printf "%.1f",$1/1e6}')
  printf '%s\t%s\t%s\t%s\n' "$s" "$bn" "$n" "$mb" >> "$B/bin_counts.tsv"
  printf '  %-8s %-10s %8s %12s\n' "$s" "$bn" "$n" "$mb"
done; done
echo
echo "  workspace_1 for reference:  SED metabat 20 / vamb 21 / taxvamb 7"
echo "                              LIQ metabat ~74 / vamb 55"

# ------------------------------------------------------------ VALIDATION
echo
echo "######## VALIDATION GATE ########"
for s in SED LIQ; do for bn in metabat vamb semibin; do
  n=$(ls -1 "$B/$bn/${s}_${bn}_"*.fa 2>/dev/null | wc -l)
  check "$s $bn produced bins ($n)" $([ "$n" -gt 0 ] && echo 0 || echo 1)
done; done
# naming compliance
bad=$(ls -1 "$B"/*/*.fa 2>/dev/null | grep -vcE '/(SED|LIQ)_(metabat|vamb|semibin)_[0-9]{3}\.fa$' || true)
check "all bin filenames follow the convention (${bad} violations)" $([ "${bad:-0}" -eq 0 ] && echo 0 || echo 1)
# every contig in every bin must carry its sample prefix
for s in SED LIQ; do
  wrong=$(cat "$B"/*/${s}_*_*.fa 2>/dev/null | grep '^>' | grep -cv "^>${s}_" || true)
  check "$s bins contain only ${s}_ contigs (${wrong} foreign)" $([ "${wrong:-0}" -eq 0 ] && echo 0 || echo 1)
done

runlog "06_bin complete"
echo
echo "######## 06 DONE — next: 07_consensus.sh (DAS_Tool) ########"
echo
echo "If SED bin counts are low across all three binners, that is a DATA"
echo "limit (median sediment read length 364 bp), not a method failure —"
echo "escalation options are TaxVAMB and a lowered DAS_Tool score threshold."
