set -uo pipefail; R=$(cd $(dirname $0); pwd); S=$(dirname $(dirname $R)); EB=$S/env/.pixi/envs/default/bin; CE=${CONDA_ENVS_USER}
log() { echo "$(date -Is) $*" >> $S/PROGRESS.log; }
mkdir -p $S/B2_mcrA_reads/v2
tree() {
  local g=$1 o=$R/$1
  $EB/python $R/build_refs.py $g > $o.build.out 2>&1 || { log "v2 $g build_refs FAILED"; return 1; }
  echo "CMD [$(date -Is)] mafft --auto; trimal -gappyout; iqtree3 ($g)" >> $S/common/COMMANDS.log
  $EB/mafft --auto --thread 8 $o/${g}_in.faa > $o/$g.aln 2> $o/mafft.log
  $CE/trimal/bin/trimal -in $o/$g.aln -out $o/$g.trim -gappyout
  if [ $g = mcrA ]; then M="-m MFP -mset LG,WAG,JTT -B 1000 -alrt 1000"; else M="-m LG+C20+F+G -B 1000"; fi
  $CE/iqtree3/bin/iqtree3 -s $o/$g.trim $M -T 8 --prefix $o/$g -seed 42 -redo > $o/iqtree.out 2>&1 || { log "v2 $g iqtree FAILED"; return 1; }
  log "v2 $g tree done ($(grep -c '>' $o/${g}_in.faa) refs)"
}
tree mcrA && REFDIR=$R/mcrA $EB/python $S/common/place.py mcrA $S/B2_mcrA_reads/v2 > $S/B2_mcrA_reads/v2/place_mcrA.out 2>&1 || log "v2 placement mcrA FAILED"
log "v2 placement mcrA done: $(head -1 $S/B2_mcrA_reads/v2/QC_mcrA.txt 2>/dev/null)"
