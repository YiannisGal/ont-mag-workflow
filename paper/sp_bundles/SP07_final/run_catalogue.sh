#!/bin/bash
S=$(cd $(dirname $0); pwd); P=$S/B4_unbinned_phyla
until [ -s $P/existing/QC.txt ]; do sleep 120; done
echo "CMD [$(date -Is)] python scripts/bin_catalogue.py" >> $P/COMMANDS.log
$S/env/.pixi/envs/default/bin/python $P/scripts/bin_catalogue.py > $P/catalogue_run.out 2>&1; echo "$(date -Is) bin catalogue + figure rc=$?" >> $S/PROGRESS.log
