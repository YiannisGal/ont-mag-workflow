#!/usr/bin/env bash
# detached driver for the remaining SP01 stages (survives a dropped connection)
cd "$(dirname "$0")"
while pgrep -f "run_SP01.sh kofam_iso" >/dev/null; do sleep 60; done
[ -s kofam_isolates_detail.tsv ] || bash run_SP01.sh kofam_iso || { echo "STOP: kofam_iso failed"; exit 1; }
bash run_SP01.sh markers || { echo "STOP: markers/QC failed -- matrix not produced"; exit 1; }
bash run_SP01.sh nif hyd || echo "WARN: nif/hyd failed"
while pgrep -f "run_SP01.sh urease_reads" >/dev/null; do sleep 120; done
[ -s urease_reads.tsv ] || bash run_SP01.sh urease_reads
bash run_SP01.sh build
echo "PIPELINE_DONE $(date -Is)"
