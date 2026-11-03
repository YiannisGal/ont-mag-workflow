#!/usr/bin/env bash
# =============================================================================
# 02_env_smoketest.sh
#
# PURPOSE
#   Several conda envs on this machine are broken. A binary sitting in
#   <env>/bin/ proves nothing -- Python tools fail at import, R-dependent tools
#   fail on missing packages, and some envs have mismatched shared libraries.
#
#   This script EMPIRICALLY tests every tool in every candidate env and writes
#   env_resolved.sh containing only the combinations that actually ran. Nothing
#   downstream is allowed to guess.
#
#   Two tiers:
#     Tier 1 (seconds)  -- invoke the tool; catches import/library failures,
#                          which is how a broken conda env almost always fails.
#     Tier 2 (minutes)  -- run the tool on real data extracted from your own
#                          assembly; catches database and runtime failures that
#                          --help would never reveal.
#
# USAGE
#   cd ~/Documents/Metagenome_analysis_PhD/workspace_2
#   bash scripts/02_env_smoketest.sh 2>&1 | tee logs/02_smoketest_$(date +%Y%m%d_%H%M).log
#
# RUNTIME  ~10-20 min (Tier 2 dominates; checkm2 testrun is the slowest)
# OUTPUT   workspace_2/env_resolved.sh   <- sourced by 00_config.sh thereafter
#          workspace_2/00_smoketest/smoketest_report.tsv
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh" 2>/dev/null || source ./00_config.sh
export STEP="02_smoketest"
RESOLVED="$WS2/env_resolved.sh"
REPORT="$SMOKE/smoketest_report.tsv"
printf 'tool\tenv\ttier\tresult\tdetail\n' > "$REPORT"
: > "$RESOLVED"
echo "# env_resolved.sh -- generated $(date -Is) by 02_env_smoketest.sh" >> "$RESOLVED"
echo "# Only combinations that ACTUALLY RAN appear here." >> "$RESOLVED"

pass=0; fail=0; unresolved=""

# ---------------------------------------------------------------- TIER 1
# invoke <env> <tool> <args...> ; returns 0 only if the tool actually executed
try_run() {
  local env="$1"; shift
  if [ "$env" = "-" ]; then timeout 120 "$@" >/dev/null 2>&1
  else timeout 180 conda run -n "$env" "$@" >/dev/null 2>&1; fi
}

resolve() {  # resolve <tool> <invocation-args...>
  local tool="$1"; shift
  local var="CAND_${tool//-/_}"
  local cands="${!var:-}"
  local winner=""
  for env in $cands; do
    if try_run "$env" "$@"; then
      winner="$env"
      printf '%s\t%s\tT1\tPASS\t\n' "$tool" "$env" >> "$REPORT"
      break
    else
      printf '%s\t%s\tT1\tfail\t\n' "$tool" "$env" >> "$REPORT"
    fi
  done
  if [ -n "$winner" ]; then
    printf '  %-28s -> %s\n' "$tool" "$winner"
    echo "export ENV_${tool//-/_}=\"$winner\"" >> "$RESOLVED"
    pass=$((pass+1))
  else
    printf '  %-28s -> ** NO WORKING ENV ** (tried: %s)\n' "$tool" "$cands"
    unresolved="$unresolved $tool"
    fail=$((fail+1))
  fi
}

echo "######## TIER 1 — does the tool execute at all? ########"
resolve minimap2    minimap2 --version
resolve samtools    samtools --version
resolve seqkit      seqkit version
resolve flye        flye --version
resolve medaka      medaka --version
resolve metabat2    metabat2 --help
resolve jgi_summarize_bam_contig_depths jgi_summarize_bam_contig_depths --help
resolve vamb        vamb --help
resolve SemiBin2    SemiBin2 --help
resolve DAS_Tool    DAS_Tool --version
resolve checkm2     checkm2 --help
resolve gtdbtk      gtdbtk --version
resolve bakta       bakta --help
resolve prokka      prokka --version
resolve coverm      coverm --version
resolve fastANI     fastANI --version
resolve iqtree3     iqtree3 --version
resolve trimal      trimal --version
resolve mafft       mafft --version
resolve quast       quast --version
resolve barrnap     barrnap --version
resolve tRNAscan_SE tRNAscan-SE -h
resolve hmmsearch   hmmsearch -h
resolve diamond     diamond --version
resolve kraken2     kraken2 --version
resolve mmseqs      mmseqs version
resolve nextflow    nextflow -version
resolve gunc        gunc --version

# ---------------------------------------------------------------- TIER 2
echo
echo "######## TIER 2 — does it work on real data / with its database? ########"
source "$RESOLVED"
mkdir -p "$SMOKE"
T="$SMOKE"

# a small real test genome: the largest contig from the liquid assembly
if [ ! -s "$T/test_genome.fna" ]; then
  echo "--- extracting a test contig from your own assembly"
  run_tool seqkit sort -l -r "$LIQ_ASM" 2>/dev/null | head -40000 > "$T/test_genome.fna" \
    || awk '/^>/{n++} n<=1' "$LIQ_ASM" > "$T/test_genome.fna"
  printf '  test genome: %s bp\n' "$(grep -v '^>' "$T/test_genome.fna" | tr -d '\n' | wc -c)"
fi
# a small read subset
if [ ! -s "$T/test_reads.fastq" ]; then
  zcat "$LIQ_READS" | head -80000 > "$T/test_reads.fastq"
  printf '  test reads: %s\n' "$(( $(wc -l < "$T/test_reads.fastq") / 4 ))"
fi

t2() {  # t2 <name> <command...>
  local name="$1"; shift
  printf '  %-30s ' "$name"
  if timeout 3600 "$@" > "$T/${name}.log" 2>&1; then
    echo "PASS"; printf '%s\t-\tT2\tPASS\t\n' "$name" >> "$REPORT"
  else
    echo "FAIL  (see 00_smoketest/${name}.log)"; printf '%s\t-\tT2\tFAIL\t\n' "$name" >> "$REPORT"
    unresolved="$unresolved ${name}(T2)"
  fi
}

echo "--- mapping chain (minimap2 -> samtools -> depth)"
t2 minimap2_map bash -c "run_tool minimap2 -ax map-ont -t 4 '$T/test_genome.fna' '$T/test_reads.fastq' > '$T/test.sam'"
t2 samtools_sort bash -c "run_tool samtools sort -@ 4 -o '$T/test.bam' '$T/test.sam' && run_tool samtools index '$T/test.bam'"
t2 jgi_depth bash -c "run_tool jgi_summarize_bam_contig_depths --outputDepth '$T/test_depth.txt' '$T/test.bam'"

echo "--- CheckM2 (built-in end-to-end test; exercises diamond + the model files)"
t2 checkm2_testrun bash -c "run_tool checkm2 testrun --threads 4"
echo "--- CheckM2 against YOUR database on the test genome"
t2 checkm2_realdb bash -c "run_tool checkm2 predict --threads 4 --input '$T/test_genome.fna' --output-directory '$T/checkm2_out' --database_path '$CHECKM2DB' --force"

echo "--- GTDB-Tk install check (this is the R232-vs-R220 question)"
t2 gtdbtk_check_install bash -c "GTDBTK_DATA_PATH='$GTDBTK_DATA_PATH' run_tool gtdbtk check_install"
if ! grep -q "gtdbtk_check_install.*PASS" "$REPORT"; then
  echo "    retrying against release220 fallback..."
  t2 gtdbtk_check_install_R220 bash -c "GTDBTK_DATA_PATH='$GTDB_FALLBACK' run_tool gtdbtk check_install"
fi

echo "--- Medaka on the test data (confirms the model loads and torch works)"
t2 medaka_consensus bash -c "run_tool medaka_consensus -i '$T/test_reads.fastq' -d '$T/test_genome.fna' -o '$T/medaka_out' -t 4 -m '$MEDAKA_MODEL' -f"

echo "--- Bakta against your database (confirms the db type and that it loads)"
t2 bakta_annotate bash -c "run_tool bakta --db '$BAKTA_DB' --output '$T/bakta_out' --prefix test --threads 4 --skip-plot --force '$T/test_genome.fna'"

echo "--- DAS_Tool (R + pullseq dependencies; the most fragile tool here)"
t2 dastool_help bash -c "run_tool DAS_Tool --help"

echo "--- binners"
t2 metabat2_run bash -c "run_tool metabat2 -i '$T/test_genome.fna' -a '$T/test_depth.txt' -o '$T/metabat_test/bin' -m 1500 --seed 1"
t2 semibin2_help bash -c "run_tool SemiBin2 single_easy_bin --help"
t2 vamb_help bash -c "run_tool vamb bin default --help"

echo "--- others"
t2 coverm_genome bash -c "run_tool coverm genome --single '$T/test_reads.fastq' --genome-fasta-files '$T/test_genome.fna' -p minimap2-ont -t 4 -o '$T/coverm_test.tsv'"
t2 barrnap_run bash -c "run_tool barrnap --kingdom bac --threads 4 '$T/test_genome.fna'"
t2 iqtree3_run bash -c "run_tool iqtree3 --version"
t2 quast_run bash -c "run_tool quast -o '$T/quast_out' -t 4 '$T/test_genome.fna'"

# ---------------------------------------------------------------- SUMMARY
echo
echo "######## SUMMARY ########"
echo "  Tier 1 resolved: $pass tools | unresolved: $fail"
[ -n "$unresolved" ] && { echo "  NEEDS ATTENTION:"; for u in $unresolved; do echo "      $u"; done; }
echo
echo "  env_resolved.sh:"
grep '^export' "$RESOLVED" | sed 's/^/      /'
echo
echo "  Full report: $REPORT"
echo
echo "  Bakta database type (from the test annotation):"
grep -iE 'database|db version|version' "$T/bakta_annotate.log" 2>/dev/null | head -5 | sed 's/^/      /'
echo
echo "  GTDB release actually validated:"
grep -iE 'release|version|done' "$T/gtdbtk_check_install.log" 2>/dev/null | tail -5 | sed 's/^/      /'

runlog "02_smoketest: $pass resolved, $fail unresolved"
echo
echo "######## CHECK: smoketest_report.tsv, env_resolved.sh, and the SUMMARY above ########"
echo "Do not start step 03 until every tool in the day-1/2 path resolves:"
echo "  minimap2 samtools jgi_summarize_bam_contig_depths metabat2 vamb SemiBin2 DAS_Tool checkm2 medaka"
