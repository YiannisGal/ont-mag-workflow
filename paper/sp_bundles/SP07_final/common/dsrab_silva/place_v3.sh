#!/bin/bash
H=$(cd $(dirname $0); pwd); S=$(dirname $(dirname $H)); EB=$S/env/.pixi/envs/default/bin
log() { echo "$(date -Is) $*" >> $S/PROGRESS.log; }
log "START B1 v3 placement (rerun after leaf-name fix)"
p() { local g=$1; REFDIR=$H/$g $EB/python $S/common/place.py $g $S/B1_dsr_reads/v3 > $S/B1_dsr_reads/v3/place_$g.out 2>&1 || log "v3 placement $g FAILED"
  log "v3 placement $g done: $(head -2 $S/B1_dsr_reads/v3/QC_$g.txt 2>/dev/null | tr '\n' ' ')"; }
p dsrAB & p dsrB & p dsrA & wait
log "END B1 v3 placement"
