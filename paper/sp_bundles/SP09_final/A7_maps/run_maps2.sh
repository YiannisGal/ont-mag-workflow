# one pooled eggNOG-mapper run (DB read once), split per MAG by locus-tag prefix, then draw maps
[ -s eggnog/all6.emapper.annotations ] || conda run -n sp09_emapper emapper.py -i eggnog/all6.faa --itype proteins -m diamond --data_dir ${SHARED_DB_DIR}/eggNOG --output all6 --output_dir eggnog --cpu 12 --override > eggnog/all6.log 2>&1 || echo "EGGFAIL all6"
for x in ETO_LIQ_MAG_001:Y ETO_LIQ_MAG_007:Y ETO_SED_MAG_001:N ETO_LIQ_MAG_003:N ETO_LIQ_MAG_009:N ETO_LIQ_MAG_010:N; do
  m=${x%:*}; c=${x#*:}; p=$(echo ${m#ETO_} | tr -d _)_
  grep "^$p" eggnog/all6.emapper.annotations > eggnog/$m.emapper.annotations
  conda run -n sp09_bakta python -I scripts/plot_map.py $m bakta/$m/$m.gbff ${m}_circular_map x $c eggnog/$m.emapper.annotations > maplog_$m.txt 2>&1 && echo "MAPDONE $m $(tail -n1 maplog_$m.txt)" || echo "MAPFAIL $m"
done
echo MAPS_ALL_DONE
