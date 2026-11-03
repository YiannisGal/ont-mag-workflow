# waits for each Bakta result, runs eggNOG-mapper (COG), then draws the map
until grep -q "emapper-2.1" "${EMAPPER_LOG}" 2>/dev/null; do sleep 30; done
for x in ETO_LIQ_MAG_001:Y ETO_LIQ_MAG_007:Y ETO_SED_MAG_001:N ETO_LIQ_MAG_003:N ETO_LIQ_MAG_009:N ETO_LIQ_MAG_010:N; do
  m=${x%:*}; c=${x#*:}
  until [ -s bakta/$m/$m.gbff ] || grep -q "FAIL2 $m" bakta_run2.log; do sleep 60; done
  [ -s bakta/$m/$m.gbff ] || { echo "MAPFAIL $m (bakta)"; continue; }
  [ -s eggnog/$m.emapper.annotations ] || conda run -n sp09_emapper emapper.py -i bakta/$m/$m.faa --itype proteins -m diamond --data_dir ${SHARED_DB_DIR}/eggNOG --output $m --output_dir eggnog --cpu 8 --override > eggnog/$m.log 2>&1 || echo "EGGFAIL $m"
  conda run -n sp09_bakta python -I scripts/plot_map.py $m bakta/$m/$m.gbff ${m}_circular_map x $c eggnog/$m.emapper.annotations > maplog_$m.txt 2>&1 && echo "MAPDONE $m $(tail -n1 maplog_$m.txt)" || echo "MAPFAIL $m"
done
echo MAPS_ALL_DONE
