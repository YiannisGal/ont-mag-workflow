# Restart-safe AAI: if the session-managed run dies, resume missing pairs (done pairs are skipped).
export PATH=${CONDA_ENVS_USER}/sp09_aai/bin:$PATH
while true; do
  n=$(ls pairs | wc -l); [ "$n" -ge 4629 ] && { echo "AAI_COMPLETE $n"; break; }
  if ! pgrep -f "ezaai calculate" >/dev/null; then
    sleep 120; pgrep -f "ezaai calculate" >/dev/null && continue
    echo "resuming at $n $(date)"
    cat pairs.txt | xargs -P 24 -n2 sh -c 'o=pairs/$0__$1.tsv; [ -s $o ] || (cd tmp && mkdir -p $0__$1 && cd $0__$1 && ezaai calculate -i ../../db/$0.db -j ../../db/$1.db -o ../../$o -t 1 >/dev/null 2>&1; cd .. && rmdir $0__$1 2>/dev/null)'
  fi
  sleep 300
done
