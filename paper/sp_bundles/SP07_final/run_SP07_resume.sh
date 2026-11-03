#!/bin/bash
S=$(cd $(dirname $0); pwd); cd $S
st() { echo "$(date -Is) $*" >> $S/PROGRESS.log; }
part() { local n=$1; shift; st "START $n"; "$@" > $S/$n/run.out 2>&1; local rc=$?; st "END $n rc=$rc"; }
part B5_mag022_dsr bash B5_mag022_dsr/scripts/b5.sh
part B4_unbinned_phyla bash B4_unbinned_phyla/scripts/b4.sh
env/.pixi/envs/default/bin/python common/summary.py >> PROGRESS.log 2>&1
st "SP07 ALL DONE"
