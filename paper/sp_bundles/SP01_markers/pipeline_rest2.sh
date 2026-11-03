#!/usr/bin/env bash
cd "$(dirname "$0")"
bash run_SP01.sh markers || { echo "STOP: markers/QC failed -- matrix not produced"; exit 1; }
bash run_SP01.sh nif hyd || echo "WARN: nif/hyd failed"
bash run_SP01.sh build
echo "PIPELINE_DONE $(date -Is)"
