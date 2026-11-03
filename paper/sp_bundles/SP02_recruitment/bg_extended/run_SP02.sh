#!/usr/bin/env bash
# run_SP02.sh (bg_extended) -- SP02 re-run with GTDB R232 Clostridium_J/Terrisporobacter/Paraclostridium backgrounds: identity-filtered competitive read recruitment (29 MAGs + isolates + background).
# Spec: 99_thesis_audit/prompts/SP02_identity_filtered_recruitment.md
#   Two competitive indices:  A = 29 MAGs + YGCA01_Clean + YGCA01_Contaminant + YGCA10 + YGCC01 + background
#                             B = the same with YGCA07 INSTEAD of YGCA10 (99.93% ANI; never in one index)
#   RESUMABLE: every stage skips outputs that already exist; files are written as .part and renamed on success.
#   Run detached:  setsid nohup bash run_SP02.sh > run.out 2>&1 < /dev/null &
set -euo pipefail
B=${DOCS_DIR}
WS2=$B/Metagenome_analysis_PhD/workspace_2
OUT=$WS2/99_thesis_audit/SP02_recruitment/bg_extended
cd "$OUT"
exec > >(tee -a "$OUT/COMMANDS.log") 2>&1
export B WS2 OUT
PY=$WS2/envs/smetana/bin/python
COVERM=$HOME/.conda/envs/coverm_env/bin/coverm       # 0.7.0, as steps 01/08d/26
MM2=/usr/local/bin/minimap2                          # 2.28
SAM=/usr/bin/samtools                                # 1.19.2
CPU=${CPU:-24}
SED=$WS2/00_input/sediment_reads_clean.fastq.gz      # step 26/34 sediment library
LIQ=$WS2/00_input/liquid_reads.fastq.gz              # step 26 culture library ($LIQ_READS in 00_config.sh)
export PY COVERM MM2 SAM CPU SED LIQ
THR="90 95 97 99"

run(){ local rc; echo "CMD [$(date -Is)] cwd=$(pwd) :: $*"; set +e; "$@"; rc=$?; set -e
       echo "EXIT [$(date -Is)] rc=$rc :: $1"; return $rc; }
runlog(){ printf '| %s | SP02bg_%s | %s |\n' "$(date -Is)" "$1" "$2" >> "$WS2/RUNLOG.md"; }

st_bg(){
  [ -s bg/kept.tsv ] && { echo "bg cached"; return 0; }
  mkdir -p bg; run bash scripts/s0_bg.sh && runlog bg_ext "PASS: GTDB R232 backgrounds fetched + ANI-screened" || { runlog bg_ext "**FAIL: bg fetch**"; return 1; }
}
st_refs(){
  [ -s refs/A.fa ] && [ -s refs/B.fa ] && { echo "refs cached"; return 0; }
  run $PY scripts/s1_refs.py && runlog refs "PASS: indices A and B built" || { runlog refs "**FAIL: refs**"; return 1; }
}
map1(){  # map1 <index> <sample> <reads>
  local ix=$1 s=$2 r=$3 bam=bam/${1}_${2}.bam
  [ -s "$bam.bai" ] && { echo "cached $bam"; return 0; }
  mkdir -p bam
  run bash -c "set -o pipefail; $MM2 -ax map-ont --secondary=no -t $CPU refs/$ix.fa '$r' 2> bam/${ix}_${s}.mm2.log \
      | $SAM sort -@ 8 -m 2G -T bam/${ix}_${s}.tmp -o $bam.part -"
  mv "$bam.part" "$bam"; run $SAM index "$bam"
  run bash -c "$SAM flagstat -@ 8 $bam > bam/${ix}_${s}.flagstat"
}
st_map(){
  for ix in A B; do map1 $ix SED "$SED"; map1 $ix LIQ "$LIQ"; done
  runlog map "PASS: 4 BAMs (A/B x SED/LIQ)"
}
st_coverm(){
  mkdir -p filt cov
  for ix in A B; do for s in SED LIQ; do for t in $THR; do
    local fb=filt/${ix}_${s}_id$t.bam cv=cov/${ix}_${s}_id$t.tsv
    if [ ! -s "$cv" ]; then
      [ -s "$fb" ] || { run $COVERM filter -b bam/${ix}_${s}.bam -o "$fb.part.bam" \
            --min-read-percent-identity $t --min-read-aligned-percent 50 --threads $CPU && mv "$fb.part.bam" "$fb"; }
      run $COVERM genome -b "$fb" -s "~" -m covered_fraction mean trimmed_mean relative_abundance \
          --trim-min 5 --trim-max 95 --min-covered-fraction 0 --threads $CPU -o "$cv.part"
      mv "$cv.part" "$cv"
    fi
  done; done; done
  # unfiltered run (no identity / aligned-fraction filter) = the step-26 method; used for positive control 1
  for ix in A B; do for s in SED LIQ; do
    local cv=cov/${ix}_${s}_id0.tsv
    [ -s "$cv" ] || { run $COVERM genome -b bam/${ix}_${s}.bam -s "~" -m covered_fraction mean trimmed_mean relative_abundance \
        --trim-min 5 --trim-max 95 --min-covered-fraction 0 --threads $CPU -o "$cv.part" && mv "$cv.part" "$cv"; }
  done; done
  runlog coverm "PASS: coverm filter+genome at 90/95/97/99 + unfiltered"
}
st_hist(){
  mkdir -p hist
  for ix in A B; do for s in SED LIQ; do
    [ -s hist/${ix}_${s}.tsv ] || run $PY scripts/s4_hist.py $ix $s
  done; done
  runlog hist "PASS: identity histograms"
}
st_build(){
  run $PY scripts/s5_build.py && runlog build "PASS: recruitment_by_identity.tsv, QC, versions" \
    || { runlog build "**FAIL: build or QC**"; return 1; }
}
STAGES=${*:-bg refs map coverm hist build}
echo "== run_SP02.sh $(date -Is) stages: $STAGES"
for s in $STAGES; do echo "== stage $s"; "st_$s"; done
echo "== SP02_DONE $(date -Is)"
