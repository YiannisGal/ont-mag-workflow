#!/usr/bin/env bash
# run_SP04.sh -- SP04: order-level classification of archaeal V4-V5 OTUs against GTDB SSU R232 (+ SILVA 138.2 labels of ms v2).
# Spec: 99_thesis_audit/prompts/SP04_archaeal_OTU_order_classification.md
# RESUMABLE + DETACHED:  setsid nohup bash run_SP04.sh > run.out 2>&1 < /dev/null &
# Carried over from SP03 (user decisions 2026-09-29): full normalised archaeal table (otutab_norm.txt + otus.fa);
# archaeal sample YGnnAR takes the depth of YGnnBA. OTU sequences are primer-trimmed (SP03 regions/otus_arc_template.fa)
# so that primer bases do not enter identity or coverage.
set -euo pipefail
B=${DOCS_DIR}
WS2=$B/Metagenome_analysis_PhD/workspace_2
OUT=$WS2/99_thesis_audit/SP04_archaeal_orders
cd "$OUT"; mkdir -p db work
exec > >(tee -a "$OUT/COMMANDS.log") 2>&1
export B WS2 OUT
VS=$WS2/99_thesis_audit/SP03_16S_links/env/.pixi/envs/default/bin/vsearch   # 2.32.0
PY=$WS2/envs/smetana/bin/python
CPU=${CPU:-16}
OTU=$WS2/99_thesis_audit/SP03_16S_links/regions/otus_arc_template.fa
export VS PY CPU OTU
run(){ local rc; echo "CMD [$(date -Is)] cwd=$(pwd) :: $*"; set +e; "$@"; rc=$?; set -e
       echo "EXIT [$(date -Is)] rc=$rc :: $1"; return $rc; }
runlog(){ printf '| %s | SP04_%s | %s |\n' "$(date -Is)" "$1" "$2" >> "$WS2/RUNLOG.md"; }

st_db(){
  if [ ! -s db/ssu_all_r232.fna.gz ]; then
    for i in 1 2 3 4 5; do
      run curl -sf -L -o db/ssu_all_r232.fna.gz.part https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_all/ssu_all_r232.fna.gz \
        && mv db/ssu_all_r232.fna.gz.part db/ssu_all_r232.fna.gz && break; sleep 60; done
  fi
  [ -s db/MD5SUM.txt ] || run curl -sf -L -o db/MD5SUM.txt https://data.gtdb.ecogenomic.org/releases/release232/232.0/MD5SUM.txt || true
  [ -s db/gtdb_sintax.fa ] || run $PY scripts/s1_db.py
  [ -s db/gtdb_ssu.nsq ] || [ -s db/gtdb_ssu.00.nsq ] || run makeblastdb -in db/gtdb_plain.fa -dbtype nucl -out db/gtdb_ssu -parse_seqids
  runlog db "PASS: GTDB SSU R232 downloaded, sintax + BLAST databases built"
}
st_classify(){
  [ -s work/sintax.tsv ] || { run $VS --sintax "$OTU" --db db/gtdb_sintax.fa --sintax_cutoff 0.8 --strand both \
        --tabbedout work/sintax.tsv.part --threads $CPU --quiet && mv work/sintax.tsv.part work/sintax.tsv; }
  [ -s work/blast.tsv ] || { run blastn -task blastn -query "$OTU" -db db/gtdb_ssu -max_target_seqs 10 -evalue 1e-20 \
        -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore qcovs" -num_threads $CPU -out work/blast.tsv.part \
        && mv work/blast.tsv.part work/blast.tsv; }
  runlog classify "PASS: sintax (cutoff 0.8) + BLAST top hits"
}
st_build(){
  run $PY scripts/s2_build.py && runlog build "PASS: orders, depth profiles, 25-30 cm breakdown, QC" \
    || { runlog build "**FAIL: build or QC**"; return 1; }
}
STAGES=${*:-db classify build}
echo "== run_SP04.sh $(date -Is) stages: $STAGES"
for s in $STAGES; do echo "== stage $s"; "st_$s"; done
echo "== SP04_DONE $(date -Is)"
