#!/bin/bash
# SP07 master: B6 -> B1 -> B2 -> B5 -> B3 -> B4. A failing part is logged and the next part runs.
S=$(cd $(dirname $0); pwd); cd $S
st() { echo "$(date -Is) $*" | tee -a $S/PROGRESS.log; }
part() { local n=$1; shift; st "START $n"; "$@" > $S/$n/run.out 2>&1; local rc=$?; st "END $n rc=$rc"; }
[ -s B6_arc_otu22/B6_result.tsv ] && st "B6 already complete" || part B6_arc_otu22 bash B6_arc_otu22/scripts/b6.sh
part B1_dsr_reads bash B1_dsr_reads/scripts/run.sh
part B2_mcrA_reads bash B2_mcrA_reads/scripts/run.sh
part B5_mag022_dsr bash B5_mag022_dsr/scripts/b5.sh
part B3_geochem_stats bash -c "cd $S/B3_geochem_stats && echo \"CMD [\$(date -Is)] Rscript scripts/b3.R\" >> COMMANDS.log && $S/env/.pixi/envs/default/bin/Rscript scripts/b3.R"
part B4_unbinned_phyla bash B4_unbinned_phyla/scripts/b4.sh
[ -s common/summary.py ] && env/.pixi/envs/default/bin/python common/summary.py >> PROGRESS.log 2>&1
st "SP07 ALL DONE"
