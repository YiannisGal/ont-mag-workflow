# 00_config.template.sh -- generated from 00_config.sh for publication.
# Every machine-specific path is replaced by a placeholder. Set each one, then
#   source 00_config.template.sh
# Placeholders used inside the archived scripts (set these too):
#   DB_DIR           local databases directory
#   SHARED_DB_DIR    shared database disk (GTDB, eggNOG)
#   SOFTWARE_DIR     locally built software
#   ENVS_DIR         standalone environments (singlem, deepsig)
#   CONDA_ENVS_USER  user conda envs directory
#   CONDA_ENVS_SHARED / CONDA_SHARED_ROOT  shared miniconda install
#   COMBINED_WS      combined-analysis workspace
#   RAW_DATA_DIR     raw nanopore read archive
#   DOCS_DIR / HOME_DIR / SCRATCH_ROOT / SHARED_DISK  disk roots used for scratch space
#   EMAPPER_LOG      eggNOG-mapper log polled by SP09_final/A7_maps/run_maps.sh

#!/usr/bin/env bash
# =============================================================================
# 00_config.sh  --  single source of truth for workspace_2
#
# *** SOURCE THIS, DO NOT EXECUTE IT ***
#     correct:   source 00_config.sh      or   . 00_config.sh
#     wrong:     bash 00_config.sh         (runs in a subshell, sets nothing)
#
# Built from inventory_report.txt, MGeMServer, 2026-08-17.
#
# NOTE ON CONDA ENVIRONMENTS
#   Several envs on this machine are known to be broken. This file therefore
#   does NOT hard-code one env per tool. It lists CANDIDATE envs in preference
#   order; 02_env_smoketest.sh actually runs each one and writes the winners to
#   env_resolved.sh, which is sourced below. Nothing downstream assumes an env
#   works just because the binary exists in its bin/ directory.
# =============================================================================

# ----------------------------------------------------------------- LOCATIONS
export WS1="<WS1>"   # path to workspace 1 (first-pass assemblies and liquid-culture reads)
export WS2="<WS2>"   # path to workspace 2 (this analysis; step folders 01_qc ... 97_*)
export LOGS="$WS2/logs"
export VERSIONS="$WS2/versions.tsv"
export RUNLOG="$WS2/RUNLOG.md"
export SMOKE="$WS2/00_smoketest"

# ---------------------------------------------------------------- RESOURCES
# 28 logical cores (14 physical, i9-7940X); 125 GB RAM, 122 GB available.
export THREADS=24
export THREADS_LOW=8
export MAXMEM_GB=100
export TMPDIR="<TMPDIR>"   # a scratch directory with >100 GB free

# ------------------------------------------------- CANDIDATE ENVS (ordered)
# "-" means "on the base PATH, no conda run needed" and is tried first where listed.
export CAND_flye="NB_Genome -"
export CAND_medaka="medaka"
export CAND_minimap2="- NB_Genome binning_env coverm_env medaka"
export CAND_samtools="- binning_env coverm_env NB_Genome Busco"
export CAND_jgi_summarize_bam_contig_depths="- binning_env"
export CAND_metabat2="- binning_env"
export CAND_vamb="vamb_env"
export CAND_SemiBin2="SemiBin"
export CAND_DAS_Tool="dastool dastool_env"
export CAND_checkm2="checkm2 checkm2_env"
export CAND_gtdbtk="gtdbtk_env NB_Genome -"
export CAND_bakta="NB_Genome"
export CAND_coverm="coverm_env"
export CAND_fastANI="coverm_env gtdbtk_env - NB_NanoClust"
export CAND_iqtree3="iqtree3"
export CAND_trimal="trimal"
export CAND_mafft="- NB_Genome gtdbtk_env prokka_env Busco"
export CAND_quast="NB_Genome"
export CAND_barrnap="NB_Genome prokka_env gtdbtk_env prokka"
export CAND_tRNAscan_SE="NB_Genome"
export CAND_seqkit="- binning_env prokka"
export CAND_hmmsearch="- NB_Genome gtdbtk_env Busco"
export CAND_diamond="- NB_Genome checkm2 Busco"
export CAND_kraken2="tax_profile myenvname"
export CAND_mmseqs="env_mmseqs NB_Genome"
export CAND_nextflow="env_nf - NB_Genome"
export CAND_gunc="gunc_env"     # NOT YET INSTALLED -- see 01
export CAND_prokka="prokka_env prokka NB_Genome -"

# --------------------------------------------------------------- DATABASES
# !! gtdbtk is 2.7.2 but GTDBTK_DATA_PATH pointed at release220.
#    GTDB-Tk 2.7.x needs R226+; release232 is present on this machine.
#    01_setup_and_checks.sh runs `gtdbtk check_install` against this to confirm.
export GTDBTK_DATA_PATH="<GTDBTK_DATA_PATH>"   # GTDB-Tk reference data, release 232
export GTDB_RELEASE="R232"                 # <- cite in Materials & Methods
export GTDB_FALLBACK="<GTDB_FALLBACK>"   # GTDB-Tk reference data, release 226

export CHECKM2DB="<CHECKM2DB>"   # CheckM2 DIAMOND database (uniref100.KO.1.dmnd)
export BAKTA_DB="<BAKTA_DB>"   # Bakta full database directory
export GUNC_DB="<GUNC_DB>"   # GUNC progenomes 2.1 DIAMOND database
export KRAKEN2_DB="<KRAKEN2_DB>"   # Kraken2 standard 16 GB database
export DRAM_DB="<DRAM_DB>"   # DRAM database directory

# ------------------------------------------------------------ MEDAKA MODEL
export MEDAKA_MODEL="r1041_e82_400bps_sup_v5.0.0"

# ------------------------------------------------------------- INPUT FILES
export SED_READS="$WS2/00_input/sediment_reads_clean.fastq.gz"   # RBK_test 16S amplicon run excluded 2026-08-23, see 13f
export LIQ_READS="$WS1/liquid_cultures/liquid_cultures_dedup.fastq.gz"
export SED_ASM="$WS1/sediment/assembly/assembly.fasta"        # 28,572 contigs
export LIQ_ASM="$WS1/liquid_cultures/assembly/assembly.fasta" # 19,290 contigs
export MICROCOSM_DIR="$WS1/liquid_cultures/profile/fastqs/Liquid_Cultures_Metagenome"

# DECISION 2026-08-17 (user): keep ALL sequencing runs pooled, in both samples.
#   Rationale: maximise read depth for assembly, polishing and MAG recovery.
#   Consequence: no cross-sample or per-depth quantitative comparison is made
#   from the metagenome. Abundance/diversity comes from the amplicon chapter.
#   The pooled *_dedup files below already contain every run, so nothing is
#   excluded -- including the PCR-amplified sediment libraries, which is
#   acceptable precisely because no abundance claim rests on them.
#
#   CoverM is still run (09_abundance), but ONLY for within-sample statistics:
#     - mean coverage per MAG        -> a required column in Table 1
#     - genome breadth per MAG       -> confirms a MAG is genuinely present
#     - % of reads recruited to MAGs -> how much of the community the
#                                       catalogue captures (reviewers ask)
#   No between-sample or between-microcosm comparison.

mkdir -p "$LOGS" "$TMPDIR" "$SMOKE"

# ------------------------------------------- RESOLVED ENVS (from smoketest)
# env_resolved.sh sets ENV_<tool>=<working env>. Until the smoketest has run,
# nothing downstream should execute.
if [ -f "$WS2/env_resolved.sh" ]; then
  source "$WS2/env_resolved.sh"
else
  echo "[config] WARNING: env_resolved.sh not found -- run scripts/02_env_smoketest.sh first"
fi

# ------------------------------------------------------------- SHARED FUNCS
# run a tool in whichever env the smoketest proved working
run_tool() {  # run_tool <tool> <args...>
  local tool="$1"; shift
  local var="ENV_${tool//-/_}"
  local env="${!var:-}"
  if [ -z "$env" ]; then echo "  !! no working env resolved for $tool" >&2; return 127; fi
  if [ "$env" = "-" ]; then "$tool" "$@"; else conda run --no-capture-output -n "$env" "$tool" "$@"; fi
}

record_version() {  # record_version <label> <env|-> <command...>
  local label="$1"; shift; local env="$1"; shift
  local out
  if [ "$env" = "-" ]; then out=$("$@" 2>&1 | grep -viE 'warning|^$' | head -2 | tr '\n' ' ')
  else out=$(conda run -n "$env" "$@" 2>&1 | grep -viE 'warning|^$' | head -2 | tr '\n' ' '); fi
  out=$(echo "$out" | sed 's/\t/ /g' | cut -c1-180)
  [ -f "$VERSIONS" ] || printf 'step\ttool\tversion\tenv\tdate\n' > "$VERSIONS"
  printf '%s\t%s\t%s\t%s\t%s\n' "${STEP:-NA}" "$label" "$out" "$env" "$(date -Is)" >> "$VERSIONS"
  echo "  [version] $label = $out"
}

runlog() {
  [ -f "$RUNLOG" ] || printf '# workspace_2 run log\n\n| when | step | event |\n|---|---|---|\n' > "$RUNLOG"
  printf '| %s | %s | %s |\n' "$(date -Is)" "${STEP:-NA}" "$*" >> "$RUNLOG"
  echo "  [runlog] $*"
}

check() {  # check <description> <0|1>
  if [ "$2" -eq 0 ]; then echo "  [PASS] $1"; runlog "PASS: $1"; return 0
  else echo "  [FAIL] $1"; runlog "**FAIL**: $1"; return 1; fi
}

echo "[config] WS2=$WS2  THREADS=$THREADS  GTDB=$GTDB_RELEASE"

# --- appended by 01_setup: all inputs now live inside workspace_2 -----------
export SED_READS="$WS2/00_input/sediment_reads_clean.fastq.gz"   # RBK_test 16S amplicon run excluded 2026-08-23, see 13f
export LIQ_READS="$WS2/00_input/liquid_reads.fastq.gz"
export SED_ASM="$WS2/00_input/sediment_assembly.fasta"
export LIQ_ASM="$WS2/00_input/liquid_assembly.fasta"

# --- functions must be exported so `bash -c` subshells can see them ---------
export -f run_tool record_version runlog check
# =============================================================================
# 00_config_patch.sh  --  APPEND THIS TO 00_config.sh
#
#   cat 00_config_patch.sh >> 00_config.sh
#
# WHY
#   A failed `conda install` into ${CONDA_ENVS_SHARED}/NB_Genome rolled
#   back and left conda's metadata inconsistent, so `conda run -n NB_Genome`
#   now fails with DirectoryNotACondaEnvironmentError. The environment's
#   binaries are intact -- only conda's bookkeeping is broken, and that
#   directory is not writable by this user so it cannot be repaired here.
#
#   These replacement functions stop using `conda run` altogether. They locate
#   the environment prefix directly and invoke the binary from its bin/ with
#   PATH and LD_LIBRARY_PATH set. This is how workspace_1 originally invoked
#   Flye (see flye.log: Cmd: ${CONDA_ENVS_SHARED}/NB_Genome/bin/flye).
#
#   Side benefits: faster (no conda run overhead), immune to any further conda
#   metadata damage, and works for environments owned by other users.
#
#   Bash uses the LAST definition of a function, so appending this file
#   overrides the earlier definitions without editing them.
# =============================================================================

# roots to search for conda environment prefixes, in order
export CONDA_ROOTS="<CONDA_ROOTS>"   # space-separated list of conda envs directories to search

# resolve an env NAME to its prefix directory
env_prefix() {  # env_prefix <envname>
  local name="$1" root
  for root in $CONDA_ROOTS; do
    [ -d "$root/$name" ] && { echo "$root/$name"; return 0; }
  done
  return 1
}

# run a tool in whichever env the smoketest proved working, WITHOUT conda run
run_tool() {  # run_tool <tool> <args...>
  local tool="$1"; shift
  local var="ENV_${tool//-/_}"
  local env="${!var:-}"

  if [ -z "$env" ]; then
    echo "  !! no working env resolved for $tool" >&2; return 127
  fi

  # base PATH
  if [ "$env" = "-" ]; then
    command "$tool" "$@"; return $?
  fi

  local prefix
  prefix=$(env_prefix "$env") || { echo "  !! env prefix not found for '$env'" >&2; return 127; }

  if [ ! -x "$prefix/bin/$tool" ]; then
    echo "  !! $tool not executable at $prefix/bin/$tool" >&2; return 127
  fi

  PATH="$prefix/bin:$PATH" \
  LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" \
  CONDA_PREFIX="$prefix" \
    "$prefix/bin/$tool" "$@"
}

# same approach for version capture
record_version() {  # record_version <label> <env|-> <command...>
  local label="$1"; shift; local env="$1"; shift
  local tool="$1"; shift
  local out prefix

  if [ "$env" = "-" ]; then
    out=$(command "$tool" "$@" 2>&1 | grep -viE 'warning|^$' | head -2 | tr '\n' ' ')
  else
    prefix=$(env_prefix "$env")
    if [ -n "$prefix" ] && [ -x "$prefix/bin/$tool" ]; then
      out=$(PATH="$prefix/bin:$PATH" LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" \
            "$prefix/bin/$tool" "$@" 2>&1 | grep -viE 'warning|^$' | head -2 | tr '\n' ' ')
    else
      out="NOT FOUND in env '$env'"
    fi
  fi

  out=$(echo "$out" | sed 's/\t/ /g' | cut -c1-180)
  [ -f "$VERSIONS" ] || printf 'step\ttool\tversion\tenv\tdate\n' > "$VERSIONS"
  printf '%s\t%s\t%s\t%s\t%s\n' "${STEP:-NA}" "$label" "$out" "$env" "$(date -Is)" >> "$VERSIONS"
  echo "  [version] $label = $out"
}

export -f env_prefix run_tool record_version runlog check 2>/dev/null || true

echo "[config] patched: bypassing 'conda run', invoking binaries by env prefix"
export CHECKM2_DB="${CHECKM2_DB:-$CHECKM2DB}"   # alias: scripts use CHECKM2_DB
export GENOMAD_DB="<GENOMAD_DB>"   # geNomad database directory
export GENOMAD_DB="<GENOMAD_DB>"   # geNomad database directory
