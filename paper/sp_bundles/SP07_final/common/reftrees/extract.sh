set -uo pipefail; cd $(dirname $0); EB=../../env/.pixi/envs/default/bin; D=${SHARED_DB_DIR}/MicrobeAnnotator_DB/protein_db
echo "CMD [$(date -Is)] seqkit grep -f candidate_ids.txt uniprot_sprot.fasta uniprot_trembl.fasta (Release 2026_03)" >> ../COMMANDS.log
$EB/seqkit grep -j 16 -f candidate_ids.txt $D/uniprot_sprot.fasta $D/uniprot_trembl.fasta -o candidates.faa
echo "EXIT [$(date -Is)] rc=$? seqkit grep" >> ../COMMANDS.log; echo DONE > extract.done
