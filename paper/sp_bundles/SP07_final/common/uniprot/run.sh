set -uo pipefail
U=$(cd $(dirname $0); pwd); S=$(dirname $(dirname $U)); EB=$S/env/.pixi/envs/default/bin; SP05=$(dirname $S)/SP05_dsr_mcr
cd $U
for db in sprot trembl; do
  [ -s hits_$db.tsv ] && continue
  echo "CMD [$(date -Is)] diamond blastx reads vs SP05 lca/uniprot_$db.dmnd" >> ../COMMANDS.log
  $EB/diamond blastx -q ../reads/reads.fa -d $SP05/lca/uniprot_$db.dmnd --frameshift 15 --more-sensitive --evalue 1e-10 -k 10 \
     --threads 24 -b 6 -c 1 --tmpdir $S/tmp --quiet -f 6 qseqid sseqid pident length evalue bitscore stitle -o hits_$db.tsv.part && mv hits_$db.tsv.part hits_$db.tsv
  echo "EXIT [$(date -Is)] rc=$? :: diamond $db" >> ../COMMANDS.log
done
$EB/python lca.py $SP05/lca/taxdump
