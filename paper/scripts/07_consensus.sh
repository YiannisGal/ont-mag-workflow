#!/usr/bin/env bash
# =============================================================================
# 07_consensus.sh  --  DAS_Tool consensus binning
#
# PURPOSE
#   Three binners disagree, and that disagreement is information. DAS_Tool
#   scores every candidate bin on single-copy gene completeness/redundancy and
#   then assigns each CONTIG to at most one winning bin. Two consequences that
#   matter here:
#
#     - The LIQ catch-all bins (metabat 44.5 Mb, vamb 19.8 Mb) score near zero
#       on SCG redundancy and lose their contigs to better-resolved competitors,
#       most likely SemiBin2's, whose largest LIQ bin is only 5.43 Mb. This is
#       the designed remedy for the 06b oversize FAIL, not a workaround.
#     - Output bins are contig-disjoint, so no contig is double-counted in any
#       downstream abundance, annotation or phylogenomic analysis.
#
#   workspace_1's sediment consensus drew 5 of 5 bins from a single binner,
#   which is the signature of an under-diverse candidate set. With three
#   genuinely independent binners on two-sample differential coverage, the
#   per-binner contribution below is a real measurement. If it again comes out
#   as one binner sweeping the board, that is worth knowing before annotation.
#
# SCORE THRESHOLD
#   DAS_Tool's default 0.5 is used and is what you report. It is deliberately
#   NOT lowered to inflate sediment recovery: CheckM2 + GUNC at step 08 do the
#   real quality filtering against MIMAG, and admitting weak bins here only
#   moves the rejection downstream while making the method harder to defend.
#   Override for a sensitivity check with:  DAS_SCORE=0.4 bash scripts/07_...
#
# INPUTS   03_binning/{metabat,vamb,semibin}/<SAMPLE>_<binner>_NNN.fa
#          04_polish/<SAMPLE>_asm_final.fasta
# OUTPUTS  05_consensus/<SAMPLE>_DASTool_bins/   (raw DAS_Tool names)
#          05_consensus/bins/<SAMPLE>_das_<binner>_NNN.fa   (naming spec)
#          05_consensus/<SAMPLE>_das_summary.tsv
#          05_consensus/DAS_CROSSWALK.tsv
#
# RUNTIME  1-3 h, dominated by Prodigal + DIAMOND on the LIQ assembly.
#
# USAGE
#   cd ${WS2}
#   bash scripts/07_consensus.sh 2>&1 | tee logs/07_consensus_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="07_consensus"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
B="$WS2/03_binning"; C="$WS2/05_consensus"
mkdir -p "$C/contig2bin" "$C/bins"

SCORE="${DAS_SCORE:-0.5}"

echo "######## 07 — DAS_TOOL CONSENSUS ########"
echo "  started: $(date -Is)  threads: $THREADS  score_threshold: $SCORE"
record_version "DAS_Tool" "${ENV_DAS_Tool:-das_tool}" DAS_Tool --version
record_version "diamond"  "${ENV_diamond:--}"          diamond --version
record_version "prodigal" "${ENV_prodigal:--}"         prodigal -v
runlog "07_consensus started (score_threshold=$SCORE, 3 binners)"

# --------------------------------------------------------------- PREFLIGHT
# DAS_Tool 1.1.x shells out to prodigal/diamond and (older builds) R.
# Fail here with a clear message rather than 40 minutes in.
for dep in prodigal diamond; do
  run_tool "$dep" --version >/dev/null 2>&1 || run_tool "$dep" -v >/dev/null 2>&1 \
    || echo "  !! WARNING: '$dep' not callable via run_tool — DAS_Tool may fail"
done

# --write_bins is a store_true flag in DAS_Tool >=1.1.3 but takes '1' in older
# releases. Detect rather than guess.
WB="--write_bins"
if run_tool DAS_Tool --help 2>&1 | grep -qE '\-\-write_bins[= ]+(1|<)'; then
  WB="--write_bins 1"
fi
echo "  write_bins syntax: $WB"

# ------------------------------------------------- CONTIG -> BIN TABLES
# Built directly from the collected FASTAs rather than via
# Fasta_to_Contig2Bin.sh, so this does not depend on the helper being on PATH
# and the bin label is guaranteed to be the naming-spec filename stem.
make_c2b() {  # make_c2b <SAMPLE> <binner>
  local s="$1" bn="$2"
  local out="$C/contig2bin/${s}_${bn}.contigs2bin.tsv"
  : > "$out"
  shopt -s nullglob
  local f n=0
  for f in "$B/$bn/${s}_${bn}_"*.fa; do
    local lab; lab=$(basename "$f" .fa)
    grep '^>' "$f" | sed 's/^>//' | awk -v b="$lab" '{print $1"\t"b}' >> "$out"
    n=$((n+1))
  done
  shopt -u nullglob
  printf '    %-8s %3d bins -> %7d contigs\n' "$bn" "$n" "$(wc -l < "$out")"
  [ -s "$out" ] || return 1
}

# ------------------------------------------------------------- RUN DAS_TOOL
run_das() {  # run_das <SAMPLE> <assembly>
  local s="$1" asm="$2"
  echo
  echo "================================================================"
  echo "  SAMPLE $s"
  echo "================================================================"
  [ -s "$asm" ] || { echo "  !! missing assembly $asm"; return 1; }

  if [ -d "$C/${s}_DASTool_bins" ] && ls "$C/${s}_DASTool_bins"/*.fa >/dev/null 2>&1; then
    echo "  already done — delete $C/${s}_DASTool_bins to force"
    return 0
  fi

  echo "  building contig2bin tables:"
  local labels=() tables=()
  for bn in metabat vamb semibin; do
    if make_c2b "$s" "$bn"; then
      labels+=("$bn"); tables+=("$C/contig2bin/${s}_${bn}.contigs2bin.tsv")
    else
      echo "    !! $bn produced no table for $s — excluded from consensus"
    fi
  done
  if [ "${#tables[@]}" -lt 2 ]; then
    echo "  !! fewer than 2 binners available for $s — consensus is meaningless"
    return 1
  fi

  local il ll
  il=$(IFS=,; echo "${tables[*]}")
  ll=$(IFS=,; echo "${labels[*]}")
  echo "  running DAS_Tool on: $ll"
  echo "  start: $(date -Is)"

  run_tool DAS_Tool \
      -i "$il" \
      -l "$ll" \
      -c "$asm" \
      -o "$C/$s" \
      --search_engine diamond \
      --score_threshold "$SCORE" \
      --write_bin_evals \
      $WB \
      -t "$THREADS"
  local rc=$?
  echo "  finished: $(date -Is)  (exit $rc)"
  return $rc
}

run_das SED "$WS2/04_polish/SED_asm_final.fasta"; sed_rc=$?
run_das LIQ "$WS2/04_polish/LIQ_asm_final.fasta"; liq_rc=$?

# ------------------------------------------------- RENAME TO NAMING SPEC
# DAS_Tool keeps the winning bin's original label, e.g. SED_metabat_007.
# NAMING_CONVENTION.md wants the consensus stage visible: SED_das_metabat_007.
echo
echo "######## RENAMING CONSENSUS BINS (NAMING_CONVENTION.md) ########"
printf 'das_bin\tsource_binner\tsource_bin\tn_contigs\tbp\tsample\n' > "$C/DAS_CROSSWALK.tsv"
for s in SED LIQ; do
  src="$C/${s}_DASTool_bins"
  [ -d "$src" ] || { echo "  $s: no DASTool_bins directory"; continue; }
  rm -f "$C/bins/${s}_das_"*.fa
  n=0
  shopt -s nullglob
  for f in "$src"/*.fa "$src"/*.fasta; do
    [ -s "$f" ] || continue
    stem=$(basename "$f"); stem="${stem%.fa}"; stem="${stem%.fasta}"
    # stem is <SAMPLE>_<binner>_<NNN>; insert das_ after the sample code
    if [[ "$stem" =~ ^(SED|LIQ)_(metabat|vamb|semibin)_([0-9]+)$ ]]; then
      new="${BASH_REMATCH[1]}_das_${BASH_REMATCH[2]}_${BASH_REMATCH[3]}"
      binner="${BASH_REMATCH[2]}"
    else
      n=$((n+1)); printf -v pad '%03d' "$n"
      new="${s}_das_unknown_${pad}"; binner="unknown"
      echo "    ?? unexpected DAS_Tool bin name: $stem"
    fi
    cp "$f" "$C/bins/${new}.fa"
    nc=$(grep -c '^>' "$C/bins/${new}.fa")
    bp=$(grep -v '^>' "$C/bins/${new}.fa" | tr -d '\n' | wc -c)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$new" "$binner" "$stem" "$nc" "$bp" "$s" \
      >> "$C/DAS_CROSSWALK.tsv"
  done
  shopt -u nullglob
  echo "  $s: $(ls -1 "$C/bins/${s}_das_"*.fa 2>/dev/null | wc -l) consensus bins renamed"
done

# ------------------------------------------------------------- SUMMARY
echo
echo "######## CONSENSUS BINS PER SAMPLE ########"
printf '  %-6s %6s %11s %11s %10s\n' sample bins total_Mbp largest_Mb median_Mb
for s in SED LIQ; do
  awk -F'\t' -v s="$s" 'NR>1 && $6==s {a[++k]=$5; t+=$5} END{
      if(k==0){printf "  %-6s %6d %11s %11s %10s\n", s,0,"-","-","-"; exit}
      n=asort(a);
      printf "  %-6s %6d %11.1f %11.2f %10.2f\n", s,k,t/1e6,a[k]/1e6,a[int((k+1)/2)]/1e6
    }' "$C/DAS_CROSSWALK.tsv"
done

echo
echo "######## WHICH BINNER WON — THE workspace_1 QUESTION ########"
echo "  (workspace_1 sediment took 5 of 5 from a single binner. Watch for that.)"
printf '  %-6s %-9s %6s %8s\n' sample binner bins share
for s in SED LIQ; do
  tot=$(awk -F'\t' -v s="$s" 'NR>1 && $6==s' "$C/DAS_CROSSWALK.tsv" | wc -l)
  for bn in metabat vamb semibin unknown; do
    n=$(awk -F'\t' -v s="$s" -v b="$bn" 'NR>1 && $6==s && $2==b' "$C/DAS_CROSSWALK.tsv" | wc -l)
    [ "$n" -eq 0 ] && continue
    pct=$(awk -v a="$n" -v b="$tot" 'BEGIN{printf "%.0f%%",(b>0?100*a/b:0)}')
    printf '  %-6s %-9s %6s %8s\n' "$s" "$bn" "$n" "$pct"
  done
done
awk -F'\t' 'NR>1{print $6"\t"$2}' "$C/DAS_CROSSWALK.tsv" | sort | uniq -c \
  | awk '{printf "%s\t%s\t%s\n",$2,$3,$1}' > "$C/binner_contribution.tsv"

# ------------------------------------------------------------ VALIDATION
echo
echo "######## VALIDATION GATE ########"
check "DAS_Tool exited cleanly (SED)" "$sed_rc"
check "DAS_Tool exited cleanly (LIQ)" "$liq_rc"

for s in SED LIQ; do
  n=$(ls -1 "$C/bins/${s}_das_"*.fa 2>/dev/null | wc -l)
  check "$s consensus bins produced ($n)" $([ "$n" -gt 0 ] && echo 0 || echo 1)
done

# naming compliance
bad=$(ls -1 "$C/bins"/*.fa 2>/dev/null \
  | grep -vcE '/(SED|LIQ)_das_(metabat|vamb|semibin)_[0-9]{3}\.fa$' || true)
check "consensus bin names follow the convention (${bad} violations)" \
  $([ "${bad:-0}" -eq 0 ] && echo 0 || echo 1)

# THE POINT OF DAS_TOOL: contigs must be disjoint across bins, per sample
for s in SED LIQ; do
  dup=$(cat "$C/bins/${s}_das_"*.fa 2>/dev/null | grep '^>' | sort | uniq -d | wc -l)
  check "$s consensus bins are contig-disjoint (${dup} shared contigs)" \
    $([ "$dup" -eq 0 ] && echo 0 || echo 1)
done

# the 06b oversize problem should be gone
over=$(awk -F'\t' 'NR>1 && $5>15000000{n++} END{print n+0}' "$C/DAS_CROSSWALK.tsv")
check "no consensus bin exceeds 15 Mb (${over} oversize)" \
  $([ "$over" -eq 0 ] && echo 0 || echo 1)

# sample prefixes preserved
for s in SED LIQ; do
  wrong=$(cat "$C/bins/${s}_das_"*.fa 2>/dev/null | grep '^>' | grep -cv "^>${s}_" || true)
  check "$s consensus bins contain only ${s}_ contigs (${wrong} foreign)" \
    $([ "${wrong:-0}" -eq 0 ] && echo 0 || echo 1)
done

# more than one binner should contribute, or the candidate set was under-diverse
for s in SED LIQ; do
  nb=$(awk -F'\t' -v s="$s" 'NR>1 && $6==s{print $2}' "$C/DAS_CROSSWALK.tsv" | sort -u | wc -l)
  check "$s consensus draws on >1 binner (${nb} binners)" \
    $([ "$nb" -gt 1 ] && echo 0 || echo 1)
done

runlog "07_consensus complete (score=$SCORE, sed_rc=$sed_rc liq_rc=$liq_rc)"
echo
echo "######## 07 DONE — next: 08_binqc.sh (CheckM2 + GUNC + MIMAG) ########"
echo
echo "Step 08 carries the number this whole re-run exists for: CheckM2's"
echo "Average_Gene_Length. workspace_1 median was 244 aa; 290-320 aa is the"
echo "expected prokaryotic range. That single column tests the --nano-raw +"
echo "no-Medaka diagnosis."
