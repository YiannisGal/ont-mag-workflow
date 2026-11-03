#!/usr/bin/env bash
# =============================================================================
# 01_setup_and_checks.sh
#
# PURPOSE
#   (a) build the workspace_2 tree
#   (b) COPY every input file into workspace_2 so the analysis is fully
#       self-contained -- nothing outside this directory is read during the run
#   (c) verify each copy with md5 against the source
#   (d) resolve the three open questions before any compute starts:
#         - is GTDB release232 compatible with gtdbtk 2.7.2?
#         - which of the two bakta_db directories is the real one?
#         - is GUNC installed? (database is present, binary is not)
#   (e) capture corrected version strings for the ~12 tools that reported
#       junk in the first inventory pass
#
# USAGE  (see STEP_BY_STEP_SERVER_GUIDE.md)
#   cd ${WS2}
#   bash scripts/01_setup_and_checks.sh 2>&1 | tee logs/01_setup_$(date +%Y%m%d_%H%M).log
#
# RUNTIME  ~10-15 min, almost all of it copying 11 GB of reads.
#          RUN IT INSIDE tmux -- see the guide.
# WRITES   only inside workspace_2. workspace_1 is never modified.
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="01_setup"

echo "######## 01 — SETUP, INPUT COPY, PRE-FLIGHT CHECKS ########"
echo "  source (read-only): $WS1"
echo "  target            : $WS2"

# ------------------------------------------------------------------ 1. TREE
echo; echo "--- 1/6 building directory tree"
mkdir -p "$WS2"/{00_input,00_smoketest,01_qc,02_mapping,03_binning/{metabat,vamb,semibin,taxvamb,dastool},\
04_polish,05_binqc/{checkm2,gunc,mimag},06_taxonomy/{gtdbtk,ani},07_mags/{unpolished,polished},\
08_annotation/{bakta,dram},09_coverage,10_phylogenomics,11_markers,12_results,logs,scripts,tmp}
echo "  done"

# --------------------------------------------------------- 2. COPY INPUTS
# Everything the analysis touches lives inside workspace_2 from here on.
echo; echo "--- 2/6 copying input data into 00_input/  (~11 GB, several minutes)"
copy_in() {  # copy_in <src> <dest-name>
  local src="$1" dst="$WS2/00_input/$2"
  if [ ! -e "$src" ]; then echo "  !! MISSING SOURCE: $src"; return 1; fi
  if [ -e "$dst" ]; then echo "  already present, skipping: $2"; return 0; fi
  echo "  copying $2  ($(du -h "$src" | cut -f1)) ..."
  rsync -h --progress "$src" "$dst"
}
copy_in "$SED_READS" "sediment_reads.fastq.gz"
copy_in "$LIQ_READS" "liquid_reads.fastq.gz"
copy_in "$SED_ASM"   "sediment_assembly.fasta"
copy_in "$LIQ_ASM"   "liquid_assembly.fasta"

# provenance of the carried-forward Flye run -> Materials & Methods
for s in sediment liquid_cultures; do
  cp -n "$WS1/$s/assembly/params.json"       "$WS2/00_input/${s}_flye_params.json"       2>/dev/null
  cp -n "$WS1/$s/assembly/assembly_info.txt" "$WS2/00_input/${s}_flye_assembly_info.txt" 2>/dev/null
  cp -n "$WS1/$s/qc/NanoStats.txt"           "$WS2/00_input/${s}_NanoStats.txt"          2>/dev/null
  cp -n "$WS1/$s/qc/NanoStats_post_filtering.txt" "$WS2/00_input/${s}_NanoStats_post_filtering.txt" 2>/dev/null
done
echo "  Flye params / assembly_info / NanoStats copied for the methods section"

# ------------------------------------------------------ 3. VERIFY THE COPIES
echo; echo "--- 3/6 verifying copies against source (md5)"
verify() {  # verify <src> <dest-name>
  local src="$1" dst="$WS2/00_input/$2"
  local a b
  a=$(md5sum "$src" | awk '{print $1}')
  b=$(md5sum "$dst" | awk '{print $1}')
  if [ "$a" = "$b" ]; then printf '  OK    %-32s %s\n' "$2" "$a"
  else printf '  FAIL  %-32s src=%s dst=%s\n' "$2" "$a" "$b"; return 1; fi
  printf '%s  %s\n' "$b" "$2" >> "$WS2/00_input/input_checksums.md5"
}
: > "$WS2/00_input/input_checksums.md5"
copies_ok=0
verify "$SED_READS" "sediment_reads.fastq.gz" || copies_ok=1
verify "$LIQ_READS" "liquid_reads.fastq.gz"   || copies_ok=1
verify "$SED_ASM"   "sediment_assembly.fasta" || copies_ok=1
verify "$LIQ_ASM"   "liquid_assembly.fasta"   || copies_ok=1

# rewrite config paths to point INSIDE workspace_2 from now on
cat >> "$WS2/00_config.sh" <<'EOF'

# --- appended by 01_setup: all inputs now live inside workspace_2 -----------
export SED_READS="$WS2/00_input/sediment_reads.fastq.gz"
export LIQ_READS="$WS2/00_input/liquid_reads.fastq.gz"
export SED_ASM="$WS2/00_input/sediment_assembly.fasta"
export LIQ_ASM="$WS2/00_input/liquid_assembly.fasta"
EOF
echo "  00_config.sh updated: inputs now resolve inside workspace_2"

# ------------------------------------------------- 4. BASELINE INPUT STATS
echo; echo "--- 4/6 baseline statistics (pre-analysis, for the methods section)"
for s in sediment liquid; do
  a="$WS2/00_input/${s}_assembly.fasta"
  [ -e "$a" ] && printf '  %-10s assembly: %s contigs, %s bp\n' "$s" \
    "$(grep -c '^>' "$a")" "$(grep -v '^>' "$a" | tr -d '\n' | wc -c)"
done
echo "  Flye parameters carried forward:"
for s in sediment liquid_cultures; do
  f="$WS2/00_input/${s}_flye_params.json"
  [ -f "$f" ] && { echo "    [$s]"; grep -oE '"(meta|read_type|min_overlap|threads|asm_coverage|genome_size|keep_haplotypes)"[^,]*' "$f" | sed 's/^/        /'; }
done

# ------------------------------------------------- 5. DATABASE PRE-FLIGHT
echo; echo "--- 5/6 database pre-flight"

echo; echo "  [A] GTDB-Tk version vs database release  <-- CRITICAL"
conda run -n "${ENV_gtdbtk:-gtdbtk_env}" gtdbtk --version 2>&1 | head -2 | sed 's/^/      /'
echo "      releases present:"
ls -d ${SHARED_DB_DIR}/GTDB/release* 2>/dev/null | sed 's/^/        /'
echo "      configured: $GTDBTK_DATA_PATH"
echo "      running check_install (the authoritative test) ..."
conda run -n "${ENV_gtdbtk:-gtdbtk_env}" gtdbtk check_install 2>&1 | tail -20 | sed 's/^/        /'
echo "      >> if this FAILS, re-run with GTDBTK_DATA_PATH=$GTDB_FALLBACK and note which passes"

echo; echo "  [B] which bakta_db is real?"
for d in ${HOME_DIR}/bakta_db \
         ${DB_DIR}/bakta_db \
         ${DB_DIR}/bakta_db; do
  [ -d "$d" ] || continue
  printf '      %s  (%s)\n' "$d" "$(du -sh "$d" 2>/dev/null | cut -f1)"
  [ -f "$d/version.json" ] && sed 's/^/          version.json: /' "$d/version.json"
done
echo "      (bakta 'light' ~1-2 GB, 'full' ~30-80 GB — 4-5 GB means light/older schema)"

echo; echo "  [C] CheckM2 database"
ls -la "$(dirname "$CHECKM2DB")" 2>/dev/null | head -4 | sed 's/^/      /'

echo; echo "  [D] GUNC — database present, binary missing"
du -sh "$GUNC_DB" 2>/dev/null | sed 's/^/      /'
if conda run -n gunc_env gunc --version >/dev/null 2>&1; then echo "      gunc IS installed"
else
  echo "      NOT installed. Run this once (2 min, no download):"
  echo "          conda create -y -n gunc_env -c bioconda -c conda-forge gunc"
fi

echo; echo "  [E] Medaka model"
echo "      configured: $MEDAKA_MODEL"
zcat "$WS2/00_input/liquid_reads.fastq.gz" 2>/dev/null | head -1 | tr ' ' '\n' \
  | grep -iE 'model|basecall|RG' | sed 's/^/      read header says: /' \
  || echo "      (no basecaller tag in read headers — using the configured default)"

# ------------------------------------------- 6. CORRECTED VERSION CAPTURE
echo; echo "--- 6/6 corrected version strings -> versions.tsv"
record_version "flye"        "NB_Genome"  flye --version
record_version "medaka"      "medaka"     medaka --version
record_version "minimap2"    "-"          minimap2 --version
record_version "samtools"    "-"          samtools --version
record_version "vamb"        "vamb_env"   vamb --version
record_version "SemiBin2"    "SemiBin"    SemiBin2 --version
record_version "DAS_Tool"    "dastool"    DAS_Tool --version
record_version "checkm2"     "checkm2"    checkm2 --version
record_version "gtdbtk"      "gtdbtk_env" gtdbtk --version
record_version "coverm"      "coverm_env" coverm --version
record_version "fastANI"     "coverm_env" fastANI --version
record_version "iqtree3"     "iqtree3"    iqtree3 --version
record_version "trimal"      "trimal"     trimal --version
record_version "mafft"       "-"          mafft --version
record_version "quast"       "NB_Genome"  quast --version
record_version "barrnap"     "NB_Genome"  barrnap --version
record_version "diamond"     "-"          diamond --version
record_version "mmseqs"      "env_mmseqs" mmseqs version
record_version "kraken2"     "tax_profile" kraken2 --version
# the ones the first pass got wrong — these are the correct flags
record_version "seqkit"      "-"          seqkit version
record_version "prodigal"    "checkm2"    prodigal -v
record_version "blastn"      "NB_Genome"  blastn -version
record_version "bracken"     "tax_profile" bracken -v
record_version "nextflow"    "env_nf"     nextflow -version
echo "  bakta (its dependency warning suppressed):"
conda run -n NB_Genome bakta --version 2>/dev/null | tail -1 | sed 's/^/      /'
echo "  FastTree / hmmsearch print to stderr with no flag:"
conda run -n NB_Genome FastTree 2>&1 | head -1 | sed 's/^/      /'
conda run -n NB_Genome hmmsearch -h 2>&1 | sed -n '2p' | sed 's/^/      /'

# ------------------------------------------------------------ VALIDATION
echo; echo "######## VALIDATION GATE ########"
n=$(ls -1 "$WS2/00_input"/{sediment,liquid}_{reads.fastq.gz,assembly.fasta} 2>/dev/null | wc -l)
check "all 4 inputs copied into workspace_2" $([ "$n" -eq 4 ] && echo 0 || echo 1)
check "all copies md5-match the source"      "$copies_ok"
check "versions.tsv populated"               $([ -s "$VERSIONS" ] && echo 0 || echo 1)
echo "  workspace_2 size so far: $(du -sh "$WS2" 2>/dev/null | cut -f1)"
df -h "$WS2" | tail -1 | awk '{print "  free space remaining: "$4}'

runlog "01_setup complete; inputs copied and verified"
echo; echo "######## 01 DONE — now run scripts/02_env_smoketest.sh ########"
