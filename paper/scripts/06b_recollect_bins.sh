#!/usr/bin/env bash
# =============================================================================
# 06b_recollect_bins.sh  --  repair the bin collection from 06_bin.sh
#
# WHY THIS EXISTS
#   06_bin.sh ran the binners correctly but collected their output badly:
#
#   1. MetaBAT2 was run with --unbinned, which writes bin.unbinned.fa holding
#      every contig that failed to bin. The glob bin.*.fa collected it as if it
#      were a genome. SED metabat totalled 141.5 Mbp against a 144.6 Mbp
#      assembly; LIQ 231.3 against 236.8. Passed to DAS_Tool, that single file
#      competes for every contig in the assembly and corrupts the consensus.
#
#   2. SemiBin2 writes output_bins/*.fa.gz. The glob *.fa* copied the GZIP
#      files to a .fa extension, so every semibin "bin" is a compressed stream
#      with a FASTA name. This is why total_Mbp read 0.0 and grep reported
#      "binary file matches" — and why the contig-prefix validation passed
#      vacuously for those bins.
#
#   3. Bin IDs were taken from the binner's own filename. VAMB and SemiBin2
#      cluster indices run past 999, so %03d emitted 4-digit IDs and broke
#      NAMING_CONVENTION.md. 34 files violated the spec.
#
# WHAT IT DOES
#   Re-collects from 03_binning/<binner>/raw_<SAMPLE>/ — the binners are NOT
#   re-run, nothing is recomputed. Excludes MetaBAT2's non-bin outputs,
#   decompresses anything gzipped, and renumbers sequentially so every ID is
#   exactly three digits. The binner's original filename is preserved in a
#   per-binner provenance TSV, so the crosswalk back to raw output is auditable
#   (NAMING_CONVENTION.md §7) even though the IDs no longer match the binner's.
#
# RUNTIME  seconds to a couple of minutes (decompression only).
#
# USAGE
#   cd ${WS2}
#   bash scripts/06b_recollect_bins.sh 2>&1 | tee logs/06b_recollect_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="06b_recollect"
B="$WS2/03_binning"

echo "######## 06b — RE-COLLECT BINS (repair of 06_bin.sh collection) ########"
echo "  started: $(date -Is)"
runlog "06b_recollect started"

# MetaBAT2 emits these alongside real bins; none of them is a genome.
EXCLUDE_RE='unbinned|tooShort|lowDepth|tooshort|lowdepth|noclass'

# ---------------------------------------------------------------- diagnostic
echo
echo "######## WHAT THE BINNERS ACTUALLY WROTE ########"
for s in SED LIQ; do
  for bn in metabat vamb semibin; do
    for d in "$B/$bn/raw_$s" "$B/$bn/raw_$s/bins" "$B/$bn/raw_$s/output_bins"; do
      [ -d "$d" ] || continue
      n=$(find "$d" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fna' -o -name '*.fa.gz' -o -name '*.fna.gz' \) 2>/dev/null | wc -l)
      [ "$n" -eq 0 ] && continue
      ex=$(find "$d" -maxdepth 1 -type f 2>/dev/null | head -3 | xargs -r -n1 basename | tr '\n' ' ')
      printf '  %-4s %-8s %-14s %4s files   e.g. %s\n' "$s" "$bn" "$(basename "$d")" "$n" "$ex"
    done
  done
done

# ---------------------------------------------------------------- recollect
recollect() {  # recollect <SED|LIQ> <binner> <srcdir> <glob> [glob...]
  local s="$1" bn="$2" src="$3"; shift 3
  local dst="$B/$bn"
  echo
  echo "--- $s / $bn"
  if [ ! -d "$src" ]; then echo "    !! source dir missing: $src"; return 1; fi
  mkdir -p "$dst"

  # wipe only THIS sample+binner's previously collected files
  rm -f "$dst/${s}_${bn}_"*.fa

  local map="$B/${s}_${bn}_bin_provenance.tsv"
  printf 'new_name\toriginal_file\tn_contigs\tbp\n' > "$map"

  shopt -s nullglob
  local files=() pat
  for pat in "$@"; do files+=( "$src"/$pat ); done
  shopt -u nullglob
  if [ "${#files[@]}" -eq 0 ]; then echo "    !! nothing matched in $src"; return 1; fi

  # stable, human-sensible order so IDs are reproducible across re-runs
  mapfile -t files < <(printf '%s\n' "${files[@]}" | sort -V)

  local i=0 skipped=0 f base out nc bp
  for f in "${files[@]}"; do
    [ -s "$f" ] || continue
    base=$(basename "$f")
    if [[ "$base" =~ $EXCLUDE_RE ]]; then
      printf '    SKIP (not a genome bin): %s\n' "$base"
      skipped=$((skipped+1)); continue
    fi
    i=$((i+1))
    if [ "$i" -gt 999 ]; then
      echo "    !! more than 999 bins — naming spec allows 3 digits only. Stopping."
      return 1
    fi
    printf -v pad '%03d' "$i"
    out="$dst/${s}_${bn}_${pad}.fa"

    if [[ "$f" == *.gz ]]; then
      zcat "$f" > "$out" || { echo "    !! zcat failed on $base"; rm -f "$out"; i=$((i-1)); continue; }
    else
      cp "$f" "$out"
    fi

    # refuse to keep anything that is not plain-text FASTA
    if [ "$(head -c1 "$out")" != ">" ]; then
      echo "    !! $base did not decompress to FASTA — dropped"
      rm -f "$out"; i=$((i-1)); continue
    fi

    nc=$(grep -c '^>' "$out")
    bp=$(grep -v '^>' "$out" | tr -d '\n' | wc -c)
    printf '%s\t%s\t%s\t%s\n' "$(basename "$out")" "$base" "$nc" "$bp" >> "$map"
  done

  printf '    -> %s bins collected, %s non-bin files excluded\n' "$i" "$skipped"
  printf '       provenance: %s\n' "$map"
}

recollect SED metabat "$B/metabat/raw_SED"            "bin.*.fa"
recollect LIQ metabat "$B/metabat/raw_LIQ"            "bin.*.fa"
recollect SED vamb    "$B/vamb/raw_SED/bins"          "*.fna" "*.fa" "*.fna.gz" "*.fa.gz"
recollect LIQ vamb    "$B/vamb/raw_LIQ/bins"          "*.fna" "*.fa" "*.fna.gz" "*.fa.gz"
recollect SED semibin "$B/semibin/raw_SED/output_bins" "*.fa.gz" "*.fna.gz" "*.fa" "*.fna"
recollect LIQ semibin "$B/semibin/raw_LIQ/output_bins" "*.fa.gz" "*.fna.gz" "*.fa" "*.fna"

# ------------------------------------------------------------ BIN COUNTS
echo
echo "######## BINS PER BINNER PER SAMPLE (corrected) ########"
printf 'sample\tbinner\tn_bins\ttotal_Mbp\tlargest_Mb\tmedian_Mb\n' > "$B/bin_counts.tsv"
printf '  %-6s %-9s %6s %11s %11s %10s\n' sample binner bins total_Mbp largest_Mb median_Mb
for s in SED LIQ; do for bn in metabat vamb semibin; do
  map="$B/${s}_${bn}_bin_provenance.tsv"
  if [ -s "$map" ]; then
    read -r n tot big med < <(awk -F'\t' 'NR>1{a[++k]=$4; s+=$4} END{
        if(k==0){print "0 0 0 0"; exit}
        n=asort(a); printf "%d %.1f %.2f %.2f\n", k, s/1e6, a[k]/1e6, a[int((k+1)/2)]/1e6}' "$map")
  else n=0; tot=0; big=0; med=0; fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$bn" "$n" "$tot" "$big" "$med" >> "$B/bin_counts.tsv"
  printf '  %-6s %-9s %6s %11s %11s %10s\n' "$s" "$bn" "$n" "$tot" "$big" "$med"
done; done
echo
echo "  workspace_1 for reference:  SED metabat 20 / vamb 21 / taxvamb 7"
echo "                              LIQ metabat ~74 / vamb 55"
echo
echo "  largest_Mb is the catch-all detector. A prokaryotic bin above ~15 Mb is"
echo "  not a genome. Before this repair, SED metabat's largest was ~100+ Mb."

# ------------------------------------------------------------ VALIDATION
echo
echo "######## VALIDATION GATE ########"

for s in SED LIQ; do for bn in metabat vamb semibin; do
  n=$(ls -1 "$B/$bn/${s}_${bn}_"*.fa 2>/dev/null | wc -l)
  check "$s $bn produced bins ($n)" $([ "$n" -gt 0 ] && echo 0 || echo 1)
done; done

# 1. filename compliance
bad=$(ls -1 "$B"/*/*.fa 2>/dev/null | grep -vcE '/(SED|LIQ)_(metabat|vamb|semibin)_[0-9]{3}\.fa$' || true)
check "all bin filenames follow the convention (${bad} violations)" \
  $([ "${bad:-0}" -eq 0 ] && echo 0 || echo 1)

# 2. every collected bin must be plain-text FASTA (catches the gzip bug)
nontext=0
for f in "$B"/*/*.fa; do
  [ -s "$f" ] || continue
  [ "$(head -c1 "$f")" = ">" ] || nontext=$((nontext+1))
done
check "all bins are plain-text FASTA (${nontext} binary/compressed)" \
  $([ "$nontext" -eq 0 ] && echo 0 || echo 1)

# 3. no catch-all survived: no single bin above 15 Mb
oversize=$(awk -F'\t' 'NR>1 && $4>15000000{n++} END{print n+0}' "$B"/*_bin_provenance.tsv)
check "no bin exceeds 15 Mb (${oversize} oversize)" \
  $([ "$oversize" -eq 0 ] && echo 0 || echo 1)

# 4. binned fraction must be well below 100% of the assembly
for s in SED LIQ; do
  asm="$WS2/04_polish/${s}_asm_final.fasta"
  [ -s "$asm" ] || continue
  asm_bp=$(grep -v '^>' "$asm" | tr -d '\n' | wc -c)
  for bn in metabat vamb semibin; do
    map="$B/${s}_${bn}_bin_provenance.tsv"; [ -s "$map" ] || continue
    bn_bp=$(awk -F'\t' 'NR>1{s+=$4} END{print s+0}' "$map")
    pct=$(awk -v a="$bn_bp" -v b="$asm_bp" 'BEGIN{printf "%.1f",(b>0?100*a/b:0)}')
    ok=$(awk -v p="$pct" 'BEGIN{print (p>0 && p<90)?0:1}')
    check "$s $bn binned fraction plausible (${pct}% of assembly)" "$ok"
  done
done

# 5. contig prefixes — now meaningful, because the files are readable
for s in SED LIQ; do
  wrong=$(cat "$B"/*/${s}_*_*.fa 2>/dev/null | grep '^>' | grep -cv "^>${s}_" || true)
  check "$s bins contain only ${s}_ contigs (${wrong} foreign)" \
    $([ "${wrong:-0}" -eq 0 ] && echo 0 || echo 1)
done

# 6. no duplicate contig within a single bin
dups=0
for f in "$B"/*/*.fa; do
  [ -s "$f" ] || continue
  d=$(grep '^>' "$f" | sort | uniq -d | wc -l)
  [ "$d" -gt 0 ] && dups=$((dups+1))
done
check "no bin contains duplicate contig IDs (${dups} bins affected)" \
  $([ "$dups" -eq 0 ] && echo 0 || echo 1)

runlog "06b_recollect complete"
echo
echo "######## 06b DONE — next: 07_consensus.sh (DAS_Tool) ########"
