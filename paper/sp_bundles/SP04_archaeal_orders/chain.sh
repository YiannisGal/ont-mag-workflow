#!/usr/bin/env bash
cd "$(dirname "$0")"
while pgrep -f "run_SP04.sh db classify" >/dev/null; do sleep 30; done
grep -q "SP04_DONE" run.out || bash run_SP04.sh db classify > run_retry.out 2>&1
bash run_SP04.sh build > run_build.out 2>&1
