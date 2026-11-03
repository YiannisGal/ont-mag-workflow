#!/usr/bin/env bash
# =============================================================================
# 07b_score_sweep.sh  --  DAS_Tool score_threshold sensitivity analysis
#
# PURPOSE
#   The primary run at the DAS_Tool default (0.5) returned 6 sediment consensus
#   bins from 32-41 candidates per binner — barely above workspace_1's 5. The
#   question is whether that reflects the sediment (genuinely few resolvable
#   genomes) or the threshold (many bins scoring just under 0.5).
#
#   This is a decision that should be MEASURED, not argued. The sweep re-runs
#   DAS_Tool at 0.3 / 0.4 / 0.5 and reports how many bins each threshold admits
#   and how much sequence they carry. Every admitted bin is still subject to
#   CheckM2 + GUNC + MIMAG at step 08, so a lower threshold here does not lower
#   the quality bar — it only decides what gets *assessed*.
#
#   Reporting plan: the primary MAG set stays whichever threshold is chosen,
#   stated explicitly in Methods, with this sweep given as the justification.
#   That is materially more defensible than silently using a default.
#
# CHEAP BECAUSE
#   Prodigal gene calling dominates DAS_Tool's runtime and is threshold-
#   independent. The proteins written by the primary run are reused via
#   --proteins, so each extra threshold costs only the scoring pass.
#
# INPUTS   05_consensus/<SAMPLE>_proteins.faa  (from 07_consensus.sh)
#          05_consensus/contig2bin/*.tsv       (from 07_consensus.sh)
# OUTPUTS  05_consensus/sweep/<SAMPLE>_s<score>_DASTool_bins/
#          05_consensus/score_sweep.tsv
#
# RUNTIME  20-60 min total.
#
# USAGE
#   ws2 && bash scripts/07b_score_sweep.sh 2>&1 | tee logs/07b_sweep_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="07b_sweep"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
C="$WS2/05_consensus"; S="$C/sweep"; mkdir -p "$S"

SCORES="0.3 0.4 0.5"

echo "######## 07b — DAS_TOOL SCORE THRESHOLD SWEEP ########"
echo "  started: $(date -Is)  thresholds: $SCORES"
runlog "07b_score_sweep started ($SCORES)"

# locate the proteins written by the primary run — name varies by DAS_Tool build
find_proteins() {  # find_proteins <SAMPLE>
  local s="$1" p
  for p in "$C/${s}_proteins.faa" "$C/${s}.faa" "$C/${s}_proteins.faa.all"; do
    [ -s "$p" ] && { echo "$p"; return 0; }
  done
  p=$(ls -1 "$C/${s}"*proteins*.faa 2>/dev/null | head -1)
  [ -n "$p" ] && { echo "$p"; return 0; }
  return 1
}

sweep() {  # sweep <SAMPLE> <assembly>
  local s="$1" asm="$2"
  echo
  echo "================================================================"
  echo "  SAMPLE $s"
  echo "================================================================"

  local il ll tables=() labels=()
  for bn in metabat vamb semibin; do
    local t="$C/contig2bin/${s}_${bn}.contigs2bin.tsv"
    [ -s "$t" ] && { tables+=("$t"); labels+=("$bn"); }
  done
  il=$(IFS=,; echo "${tables[*]}"); ll=$(IFS=,; echo "${labels[*]}")

  local prot pflag=""
  if prot=$(find_proteins "$s"); then
    pflag="--proteins $prot"
    echo "  reusing proteins: $(basename "$prot")  ($(grep -c '^>' "$prot") genes)"
  else
    echo "  !! no cached proteins found — Prodigal will re-run per threshold (slower)"
  fi

  local WB="--write_bins"
  run_tool DAS_Tool --help 2>&1 | grep -qE '\-\-write_bins[= ]+(1|<)' && WB="--write_bins 1"

  for sc in $SCORES; do
    local tag="${s}_s${sc}"
    if ls "$S/${tag}_DASTool_bins"/*.fa >/dev/null 2>&1; then
      echo "  score $sc: already done"; continue
    fi
    echo "  --- score_threshold $sc  ($(date +%H:%M))"
    run_tool DAS_Tool -i "$il" -l "$ll" -c "$asm" -o "$S/$tag" \
        --search_engine diamond --score_threshold "$sc" \
        $pflag $WB -t "$THREADS" >"$S/${tag}.log" 2>&1 \
      || echo "      !! failed at $sc — see $S/${tag}.log"
  done
}

sweep SED "$WS2/04_polish/SED_asm_final.fasta"
sweep LIQ "$WS2/04_polish/LIQ_asm_final.fasta"

# ------------------------------------------------------------------ REPORT
echo
echo "######## SCORE THRESHOLD SWEEP ########"
printf 'sample\tscore\tn_bins\ttotal_Mbp\tlargest_Mb\tmedian_Mb\tn_metabat\tn_vamb\tn_semibin\n' \
  > "$C/score_sweep.tsv"
printf '  %-5s %6s %6s %10s %11s %10s   %s\n' sample score bins total_Mbp largest_Mb median_Mb "binner mix"
for s in SED LIQ; do
  for sc in $SCORES; do
    d="$S/${s}_s${sc}_DASTool_bins"
    if ! ls "$d"/*.fa >/dev/null 2>&1; then
      printf '  %-5s %6s %6s %10s %11s %10s\n' "$s" "$sc" 0 "-" "-" "-"; continue
    fi
    read -r n tot big med < <(
      for f in "$d"/*.fa; do grep -v '^>' "$f" | tr -d '\n' | wc -c; done \
      | sort -n | awk '{a[NR]=$1; s+=$1} END{
          printf "%d %.1f %.2f %.2f\n", NR, s/1e6, a[NR]/1e6, a[int((NR+1)/2)]/1e6}')
    nm=$(ls -1 "$d" | grep -c '_metabat_' || true)
    nv=$(ls -1 "$d" | grep -c '_vamb_'    || true)
    ns=$(ls -1 "$d" | grep -c '_semibin_' || true)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$s" "$sc" "$n" "$tot" "$big" "$med" "$nm" "$nv" "$ns" >> "$C/score_sweep.tsv"
    printf '  %-5s %6s %6s %10s %11s %10s   mb=%s vamb=%s semi=%s\n' \
      "$s" "$sc" "$n" "$tot" "$big" "$med" "$nm" "$nv" "$ns"
  done
done

echo
echo "######## HOW TO READ THIS ########"
echo "  If 0.4 and 0.3 add many bins, the 0.5 result was threshold-limited and"
echo "  those bins deserve to be assessed by CheckM2 rather than discarded"
echo "  unseen. If the counts barely move, 6 sediment bins is the biology and"
echo "  the default threshold stands — which is the stronger claim to publish."
echo
echo "  NOTE: bin counts alone do not justify a threshold. The decision is made"
echo "  at step 08 on how many ADDITIONAL bins meet MIMAG medium/high quality."

runlog "07b_score_sweep complete"
echo
echo "######## 07b DONE ########"
