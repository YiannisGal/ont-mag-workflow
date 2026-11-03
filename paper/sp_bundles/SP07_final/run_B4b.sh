#!/bin/bash
S=$(cd $(dirname $0); pwd); cd $S
echo "$(date -Is) START B4b (restart)" >> PROGRESS.log
bash B4_unbinned_phyla/scripts/b4b.sh > B4_unbinned_phyla/run_b4b.out 2>&1; echo "$(date -Is) END B4b rc=$?" >> PROGRESS.log
env/.pixi/envs/default/bin/python common/summary.py >> PROGRESS.log 2>&1; rm -rf $S/pplacer_scratch; echo "$(date -Is) SP07 ALL DONE" >> PROGRESS.log
