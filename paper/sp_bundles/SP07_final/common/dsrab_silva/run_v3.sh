#!/bin/bash
# B1 v3: dsrAB@SILVA trees (GTDB r214), LG+G4+F UFBoot 1000 (C20 intractable at 600-800 taxa), then placement + leave-one-out QC.
set -uo pipefail; H=$(cd $(dirname $0); pwd); S=$(dirname $(dirname $H)); EB=$S/env/.pixi/envs/default/bin; CE=${CONDA_ENVS_USER}
log() { echo "$(date -Is) $*" >> $S/PROGRESS.log; }
log "START B1 v3 (dsrAB@SILVA)"
tree() { local g=$1 o=$H/$1
  echo "CMD [$(date -Is)] mafft --auto; trimal -gappyout; iqtree3 -m LG+G4+F -B 1000 ($g, dsrAB@SILVA)" >> $S/common/COMMANDS.log
  [ -s $o/$g.aln ] || $EB/mafft --auto --thread 6 $o/${g}_in.faa > $o/$g.aln 2> $o/mafft.log
  [ -s $o/$g.trim ] || $CE/trimal/bin/trimal -in $o/$g.aln -out $o/$g.trim -gappyout
  [ -s $o/$g.treefile ] || $CE/iqtree3/bin/iqtree3 -s $o/$g.trim -m LG+G4+F -B 1000 -T 6 --prefix $o/$g -seed 42 > $o/iqtree.out 2>&1 || { log "v3 $g iqtree FAILED"; return 1; }
  $EB/python $H/function_check.py $g >> $H/function_check.txt
  log "v3 $g tree done"
  mkdir -p $S/B1_dsr_reads/v3
  REFDIR=$o $EB/python $S/common/place.py $g $S/B1_dsr_reads/v3 > $S/B1_dsr_reads/v3/place_$g.out 2>&1 || log "v3 placement $g FAILED"
  log "v3 placement $g done: $(cat $S/B1_dsr_reads/v3/QC_$g.txt 2>/dev/null | head -2 | tr '\n' ' ')"
}
tree dsrA & tree dsrB & tree dsrAB & wait
log "END B1 v3"
