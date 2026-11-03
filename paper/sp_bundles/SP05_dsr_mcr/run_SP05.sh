#!/usr/bin/env bash
# run_SP05.sh -- SP05: (H2) identity of the unbinned dsrAB contigs LIQ_contig_12880 / LIQ_contig_19299;
#                      (H3) why ETO_LIQ_MAG_017 (Methanogenium) is absent from mcrA_verdicts.tsv.
# Spec: 99_thesis_audit/prompts/SP05_dsrAB_identity_and_Methanogenium_mcrA.md
# RESUMABLE + DETACHED:  setsid nohup bash run_SP05.sh > run.out 2>&1 < /dev/null &
# Deviations (stated in SUMMARY): no step-85 folder exists; the dsrAB reference set is step 58's (58_dsr_finish/refs.faa,
# 82 seqs), the one the step-58/59 placement used. No GTDB/NR protein set is local; LCA uses UniProt SwissProt+TrEMBL
# (MicrobeAnnotator_DB/protein_db) with NCBI taxdump lineages.
set -euo pipefail
B=${DOCS_DIR}
WS2=$B/Metagenome_analysis_PhD/workspace_2
OUT=$WS2/99_thesis_audit/SP05_dsr_mcr
cd "$OUT"; mkdir -p work trees lca
exec > >(tee -a "$OUT/COMMANDS.log") 2>&1
export B WS2 OUT
PY=$WS2/envs/smetana/bin/python
PIXI=$B/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default/bin
MAFFT=$PIXI/mafft; TRIMAL=$HOME/.conda/envs/trimal/bin/trimal; IQ=$HOME/.conda/envs/iqtree3/bin/iqtree3
HMM=/usr/bin/hmmsearch; DMND=/usr/local/bin/diamond
KOLIST=$WS2/11_kofam/ko_list; PROF=$WS2/45_kofam_full/profiles_full/profiles
UP=${SHARED_DB_DIR}/MicrobeAnnotator_DB/protein_db
CPU=${CPU:-20}
export PY MAFFT TRIMAL IQ HMM DMND KOLIST PROF CPU
run(){ local rc; echo "CMD [$(date -Is)] cwd=$(pwd) :: $*"; set +e; "$@"; rc=$?; set -e
       echo "EXIT [$(date -Is)] rc=$rc :: $1"; return $rc; }
runlog(){ printf '| %s | SP05_%s | %s |\n' "$(date -Is)" "$1" "$2" >> "$WS2/RUNLOG.md"; }

st_inputs(){
  [ -s work/contig_proteins.faa ] || run $PY scripts/s1_inputs.py
  runlog inputs "PASS: contigs, contig proteins, dsr refs, MAG_017 joined mcrA"
}
st_kofam(){
  [ -s work/kofam_all.tbl ] || { run $HMM --cpu $CPU -E 1e-5 --tblout work/kofam_all.tbl.part $WS2/45_kofam_full/kofam_all.hmm work/kofam_query.faa > /dev/null \
        && mv work/kofam_all.tbl.part work/kofam_all.tbl; }
  runlog kofam "PASS: KOfam on contig proteins, dsr refs and MAG_017 joined mcrA"
}
st_dsr(){
  run $PY scripts/s2_dsr_prep.py
  for m in dsrA dsrB; do
    [ -s trees/$m.aln ] || run bash -c "$MAFFT --auto --thread $CPU --quiet trees/${m}_in.faa > trees/$m.aln"
    [ -s trees/$m.trim ] || run $TRIMAL -in trees/$m.aln -out trees/$m.trim -gappyout
  done
  run $PY scripts/s2_concat.py
  for m in dsrA dsrB dsrAB; do
    [ -s trees/$m.treefile ] || run $IQ -s trees/$m.trim -m LG+C20+F+G -B 1000 -T $CPU --prefix trees/$m -redo --quiet
  done
  run $PY scripts/s3_dsr_place.py && runlog dsr "PASS: dsrA/dsrB/dsrAB placement + control" \
    || { runlog dsr "**FAIL: dsr placement or control**"; return 1; }
}
st_mcr(){
  run $PY scripts/s6_mcr_prep.py
  [ -s trees/mcrA.aln ] || run bash -c "$MAFFT --auto --thread $CPU --quiet trees/mcrA_in.faa > trees/mcrA.aln"
  [ -s trees/mcrA.trim ] || run $TRIMAL -in trees/mcrA.aln -out trees/mcrA.trim -gappyout
  [ -s trees/mcrA.treefile ] || run $IQ -s trees/mcrA.trim -m MFP -mset LG,WAG,JTT -B 1000 -alrt 1000 -T $CPU --prefix trees/mcrA -redo --quiet
  run $PY scripts/s7_mcr_verdict.py && runlog mcr "PASS: mag017_mcrA_verdict.tsv + 9-protein control" \
    || { runlog mcr "**FAIL: mcrA verdict or control**"; return 1; }
}
st_mag022(){
  run $PY scripts/s5_mag022.py && runlog mag022 "PASS: mag022_compatibility.tsv" || { runlog mag022 "**FAIL: MAG_022 tests**"; return 1; }
}
st_lca(){
  if [ ! -s lca/taxdump/nodes.dmp ]; then
    mkdir -p lca/taxdump
    run curl -sf -L -o lca/taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    run tar -xzf lca/taxdump.tar.gz -C lca/taxdump nodes.dmp names.dmp
  fi
  for db in uniprot_sprot uniprot_trembl; do
    if [ ! -s lca/$db.dmnd ]; then
      run $DMND makedb --in $UP/$db.fasta -d lca/$db.part --threads $CPU --quiet
      mv lca/$db.part.dmnd lca/$db.dmnd
    fi
    [ -s lca/hits_$db.tsv ] || { run $DMND blastp -q work/contig_proteins.faa -d lca/$db.dmnd --more-sensitive --evalue 1e-10 \
        -k 10 --threads $CPU --quiet -f 6 qseqid sseqid pident length evalue bitscore stitle -o lca/hits_$db.tsv.part \
        && mv lca/hits_$db.tsv.part lca/hits_$db.tsv; }
  done
  run $PY scripts/s4_lca.py && runlog lca "PASS: contig_protein_LCA.tsv" || { runlog lca "**FAIL: LCA**"; return 1; }
}
STAGES=${*:-inputs kofam dsr mcr mag022 lca}
echo "== run_SP05.sh $(date -Is) stages: $STAGES"
for s in $STAGES; do echo "== stage $s"; "st_$s" || echo "!! stage $s failed; continuing with the next"; done
echo "== SP05_DONE $(date -Is)"
