#!/usr/bin/env bash
# =============================================================================
# 03_assemble.sh  (v2)  --  read length filtering + re-assembly with --nano-hq
#
# CHANGES FROM v1 (both from the 2026-08-17 16:38 failed run)
#
#   1. TMPDIR OVERRIDE.  v1 died at the consensus stage in both samples with
#         OSError: AF_UNIX path too long
#      Flye's consensus uses multiprocessing.Manager(), which opens a Unix
#      domain socket under $TMPDIR. Unix socket paths are limited to 108 chars.
#      00_config.sh set TMPDIR="$WS2/tmp" = 78 chars, and Python appends
#      ~32 more for /pymp-XXXXXXXX/listener-XXXXXXXX. Overflow.
#      Fixed by forcing a short TMPDIR here (and in 00_config.sh).
#
#   2. READ LENGTH FILTER.  The v1 run reported:
#         sediment  8.44 Gb   read N50 3707   read N90 359
#         liquid    4.13 Gb   read N50 2657   read N90 731
#      With --min-overlap 1000, every read <1000 bp is incapable of forming an
#      overlap. It cannot contribute to the assembly but still costs memory,
#      index space and time. Sediment's N90 of 359 bp means a large fraction of
#      the dataset is dead weight.
#
#      This is also the likely reason the sediment assembly is poor
#      (N50 5.3 kb, longest 0.29 Mb) despite Flye estimating 90x overlap
#      coverage -- it is read-LENGTH limited, not depth limited.
#
#      Reads are filtered to >=1000 bp for ASSEMBLY ONLY. The unfiltered set is
#      retained and used later for polishing and mapping, where short reads
#      still contribute coverage.
#
# RUNTIME  filtering ~20-40 min; assembly ~8-12 h. RUN IN tmux.
#
# USAGE
#   cd ${WS2}
#   bash scripts/03_assemble.sh 2>&1 | tee logs/03_assemble_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="03_assemble"

# --- FIX 1: short TMPDIR for Unix domain sockets ----------------------------
export TMPDIR="${TMPDIR}"      # 33 chars, well under 108
mkdir -p "$TMPDIR"
echo "  TMPDIR=$TMPDIR (${#TMPDIR} chars — must stay well under 108)"

echo "######## 03 — READ FILTERING + RE-ASSEMBLY (--nano-hq) ########"
echo "  started: $(date -Is)   threads: $THREADS"
record_version "flye"   "${ENV_flye:-NB_Genome}" flye --version
record_version "seqkit" "${ENV_seqkit:--}"       seqkit version
runlog "03_assemble v2 started (length filter >=1000 bp, --nano-hq --meta)"

mkdir -p "$WS2/01_qc" "$WS2/03_assembly"

# ------------------------------------------------------- FIX 2: LENGTH FILTER
MINLEN=1000
filter_reads() {  # filter_reads <label> <in.fastq.gz> <out.fastq.gz>
  local label="$1" in="$2" out="$3"
  if [ -s "$out" ]; then echo "  $label: filtered reads already present, skipping"; return 0; fi
  echo
  echo "--- filtering $label to >= ${MINLEN} bp"
  echo "    before:"
  run_tool seqkit stats -a -T "$in" | sed 's/^/      /'
  run_tool seqkit seq -m "$MINLEN" -j "$THREADS_LOW" "$in" -o "$out"
  echo "    after:"
  run_tool seqkit stats -a -T "$out" | sed 's/^/      /'
  # keep the stats for the methods section
  { run_tool seqkit stats -a -T "$in"; run_tool seqkit stats -a -T "$out" | tail -1; } \
    > "$WS2/01_qc/${label}_readstats_prepost.tsv"
}

SED_FILT="$WS2/01_qc/SED_reads_min${MINLEN}.fastq.gz"
LIQ_FILT="$WS2/01_qc/LIQ_reads_min${MINLEN}.fastq.gz"
filter_reads "SED" "$SED_READS" "$SED_FILT"
filter_reads "LIQ" "$LIQ_READS" "$LIQ_FILT"

# ------------------------------------------------------------------ ASSEMBLE
assemble() {  # assemble <label> <reads> <outdir>
  local label="$1" reads="$2" out="$3"
  echo
  echo "================================================================"
  echo "  ASSEMBLING: $label"
  echo "  reads : $reads"
  echo "  start : $(date -Is)"
  echo "================================================================"
  if [ -s "$out/assembly.fasta" ]; then
    echo "  assembly.fasta already present — skipping (delete to force)"; return 0
  fi
  local resume=""
  [ -d "$out/00-assembly" ] && { resume="--resume"; echo "  previous run found -> --resume"; }

  run_tool flye \
      --nano-hq "$reads" \
      --out-dir "$out" \
      --meta \
      --threads "$THREADS" \
      --min-overlap "$MINLEN" \
      $resume
  local rc=$?
  echo "  finished: $(date -Is)  (exit $rc)"
  return $rc
}

# NOTE: the previous failed run left 00-assembly/ in both output dirs. --resume
# picks up from the disjointig stage, so the first ~20 min is not repeated.
assemble "sediment"        "$SED_FILT" "$WS2/03_assembly/sediment";        sed_rc=$?
assemble "liquid_cultures" "$LIQ_FILT" "$WS2/03_assembly/liquid_cultures"; liq_rc=$?

# ------------------------------------------------------- COMPARE OLD vs NEW
echo
echo "######## OLD (--nano-raw, unfiltered) vs NEW (--nano-hq, >=${MINLEN} bp) ########"
stats() {
  [ -s "$1" ] || { echo "MISSING"; return; }
  awk '/^>/{if(l)print l; l=0; next}{l+=length($0)}END{if(l)print l}' "$1" \
  | sort -rn | awk '{a[NR]=$1; s+=$1} END{
      h=s/2; c=0; for(i=1;i<=NR;i++){c+=a[i]; if(c>=h){n50=a[i]; break}}
      printf "%d contigs | %.1f Mbp | N50 %.1f kb | longest %.2f Mb", NR, s/1e6, n50/1000, a[1]/1e6}'
}
printf '\n  %-16s %s\n'  "sediment OLD:" "$(stats "$WS2/00_input/sediment_assembly.fasta")"
printf '  %-16s %s\n'    "sediment NEW:" "$(stats "$WS2/03_assembly/sediment/assembly.fasta")"
printf '\n  %-16s %s\n'  "liquid OLD:"   "$(stats "$WS2/00_input/liquid_assembly.fasta")"
printf '  %-16s %s\n'    "liquid NEW:"   "$(stats "$WS2/03_assembly/liquid_cultures/assembly.fasta")"

echo
echo "  circular contigs  (OLD: sediment 118, liquid 240)"
for s in sediment liquid_cultures; do
  f="$WS2/03_assembly/$s/assembly_info.txt"
  [ -s "$f" ] && printf '    %-18s NEW %s circular of %s\n' "$s" \
    "$(awk 'NR>1 && $4=="Y"' "$f" | wc -l)" "$(( $(wc -l < "$f") - 1 ))"
done

for s in sediment liquid_cultures; do
  cp -f "$WS2/03_assembly/$s/assembly_info.txt" "$WS2/03_assembly/${s}_assembly_info.txt" 2>/dev/null
  grep -m1 "Cmd:" "$WS2/03_assembly/$s/flye.log" 2>/dev/null | sed "s|^|  $s cmd: |"
done

# ------------------------------------------------------------- VALIDATION
echo
echo "######## VALIDATION GATE ########"
check "sediment assembly produced"     $([ -s "$WS2/03_assembly/sediment/assembly.fasta" ] && echo 0 || echo 1)
check "liquid assembly produced"       $([ -s "$WS2/03_assembly/liquid_cultures/assembly.fasta" ] && echo 0 || echo 1)
check "flye exited cleanly (sediment)" "$sed_rc"
check "flye exited cleanly (liquid)"   "$liq_rc"

runlog "03_assemble v2 complete (sed_rc=$sed_rc liq_rc=$liq_rc)"
echo
echo "######## 03 DONE ########"
