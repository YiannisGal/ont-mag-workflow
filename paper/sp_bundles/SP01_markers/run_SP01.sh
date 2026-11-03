#!/usr/bin/env bash
# run_SP01.sh -- SP01: provenance of the "31 MAGs" + curated diagnostic-marker matrix.
# Spec: 99_thesis_audit/prompts/SP01_marker_matrix_and_31_vs_29.md
#
#   bash run_SP01.sh [stage ...]      stages: prov kofam_iso markers nif hyd urease_ref urease_reads build
#                                     (no argument = all, in order)
#
# DEVIATION (approved by user 2026-09-29): step 45 did NOT use KofamScan exec_annotation and
# has no versions.tsv. It ran HMMER 3.4 hmmsearch over 11_kofam/profiles.tar.gz and applied
# 11_kofam/ko_list thresholds in Python. The isolates are annotated the same way here, so
# MAGs and isolates are scored by one method. Score type (full/domain) is honoured per ko_list.
set -euo pipefail
WS2=${WS2}
B=${DOCS_DIR}
OUT=$WS2/99_thesis_audit/SP01_markers
cd "$OUT"
exec > >(tee -a "$OUT/COMMANDS.log") 2>&1
export OUT WS2 B
PY=$WS2/envs/smetana/bin/python          # pandas 2.3.3; same interpreter as step 45
PIXI=$B/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default/bin   # mafft, FastTree
HMMSEARCH=/usr/bin/hmmsearch             # HMMER 3.4, as step 45 / 34
DIAMOND=/usr/local/bin/diamond           # 2.1.9, as step 34 (PATH resolution, ENV_diamond="-")
CPU=${CPU:-24}
KOLIST=$WS2/11_kofam/ko_list
PROFILES=$WS2/45_kofam_full/profiles_full/profiles
export PY PIXI HMMSEARCH DIAMOND CPU KOLIST PROFILES

run(){  # every command: timestamp, cwd, command line, exit code
  local ts rc; ts=$(date -Is)
  echo "CMD [$ts] cwd=$(pwd) :: $*"
  set +e; "$@"; rc=$?; set -e
  echo "EXIT [$(date -Is)] rc=$rc :: $1"
  return $rc
}
runlog(){ printf '| %s | SP01_%s | %s |\n' "$(date -Is)" "$1" "$2" >> "$WS2/RUNLOG.md"; }

st_prov(){
  run $PY scripts/s1_provenance.py && runlog prov "PASS: provenance_31.tsv + combined_analysis_marker_method.md" \
    || { runlog prov "**FAIL: provenance**"; return 1; }
}
st_kofam_iso(){
  mkdir -p isolates/hits
  run $PY scripts/s2_prep_isolates.py
  for c in "$WS2"/45_kofam_full/chunks/c*.hmm; do
    b=$(basename "$c" .hmm)
    [ -s "isolates/hits/$b.domtbl" ] && { echo "cached $b"; continue; }
    run $HMMSEARCH --cpu "$CPU" -E 1e-5 --tblout "isolates/hits/$b.tbl.part" \
        --domtblout "isolates/hits/$b.domtbl.part" "$c" isolates/isolates.faa > /dev/null
    mv "isolates/hits/$b.tbl.part" "isolates/hits/$b.tbl"; mv "isolates/hits/$b.domtbl.part" "isolates/hits/$b.domtbl"
  done
  run $PY scripts/s2_isolate_detail.py && runlog kofam_iso "PASS: kofam_isolates_detail.tsv" \
    || { runlog kofam_iso "**FAIL: isolate KOfam**"; return 1; }
}
st_markers(){
  run $PY scripts/s3_targets.py              # target HMMs + combined proteome of 34 genomes
  run $HMMSEARCH --cpu "$CPU" --tblout markers/targets.tbl --domtblout markers/targets.domtbl \
      markers/targets.hmm markers/all34.faa > /dev/null
  # step-16b comparison: the same targets at a flat E-value, no curated threshold
  run $PY scripts/s3_markers.py && runlog markers "PASS: marker calls + controls" \
    || { runlog markers "**FAIL: marker calls or QC controls**"; return 1; }
}
st_nif(){
  run bash scripts/s4_nif.sh && runlog nif "PASS: nif_tree_verdicts.tsv" || { runlog nif "**FAIL: nif tree**"; return 1; }
}
st_hyd(){
  run $PY scripts/s5_hyd.py && runlog hyd "PASS: hydrogenase_classes.tsv" || { runlog hyd "**FAIL: hydrogenases**"; return 1; }
}
st_urease_ref(){
  run bash scripts/s6_urease_ref.sh && runlog urease_ref "PASS: urease reference set" || { runlog urease_ref "**FAIL: urease refs**"; return 1; }
}
st_urease_reads(){
  local R=$WS2/00_input/sediment_reads_clean.fastq.gz
  [ -s urease/SED_reads.tsv ] || run $DIAMOND blastx -d urease/urease_ref.dmnd -q "$R" -o urease/SED_reads.tsv.part \
      --frameshift 15 --range-culling --top 10 -e 1e-10 --threads "$CPU" --quiet \
      -f 6 qseqid sseqid pident length evalue bitscore
  [ -s urease/SED_reads.tsv ] || mv urease/SED_reads.tsv.part urease/SED_reads.tsv
  [ -s urease/SED_nreads.txt ] || run bash -c "zcat '$R' | awk 'END{print NR/4}' > urease/SED_nreads.txt"
  run $PY scripts/s6_urease_count.py && runlog urease_reads "PASS: urease_reads.tsv" || { runlog urease_reads "**FAIL: urease reads**"; return 1; }
}
st_build(){
  run $PY scripts/s7_build.py && runlog build "PASS: markers_long/wide + versions.tsv" || { runlog build "**FAIL: build**"; return 1; }
}

STAGES=${*:-prov kofam_iso markers nif hyd urease_ref urease_reads build}
echo "== run_SP01.sh $(date -Is) stages: $STAGES"
for s in $STAGES; do echo "== stage $s"; "st_$s"; done
echo "== done $(date -Is)"
