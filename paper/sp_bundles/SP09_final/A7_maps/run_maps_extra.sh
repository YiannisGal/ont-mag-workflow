# eggNOG + pyCirclize + CGView + Proksee files for the 3 extra candidates (005, 008, 011); runs after the first eggNOG job
S=${WS2}/99_thesis_audit/SP09_final; X=$S/A13_extra_candidates
until [ -s eggnog/all6.emapper.annotations ] || grep -q "EGGFAIL" run_maps.log; do sleep 120; done
until grep -q "A13DONE bakta" $X/run_a13.log; do sleep 120; done
C="ETO_LIQ_MAG_005 ETO_LIQ_MAG_008 ETO_LIQ_MAG_011"
for m in $C; do mkdir -p bakta/$m; cp $X/bakta/$m/$m.* bakta/$m/; done
cat $(for m in $C; do echo bakta/$m/$m.faa; done) > eggnog/extra3.faa
[ -s eggnog/extra3.emapper.annotations ] || conda run -n sp09_emapper emapper.py -i eggnog/extra3.faa --itype proteins -m diamond --data_dir ${SHARED_DB_DIR}/eggNOG --output extra3 --output_dir eggnog --cpu 12 --override > eggnog/extra3.log 2>&1 || { echo "EGGFAIL extra3"; exit 1; }
B=${CONDA_ENVS_USER}/sp09_cgview/bin
for m in $C; do p=$(echo ${m#ETO_} | tr -d _)_; n=$(grep -c ">" bakta/$m/$m.fna)
  grep "^$p" eggnog/extra3.emapper.annotations > eggnog/$m.emapper.annotations
  conda run -n sp09_bakta python -I scripts/plot_map.py $m bakta/$m/$m.gbff ${m}_circular_map x N eggnog/$m.emapper.annotations > maplog_$m.txt 2>&1 && echo "MAPDONE $m" || echo "MAPFAIL $m"
  cp bakta/$m/$m.gbff proksee_upload/
  conda run -n sp09_bakta python -I scripts/cgview_inputs.py $m bakta/$m/$m.gbff eggnog/$m.emapper.annotations cgview proksee_upload > cgview/$m.prep.log 2>&1
  $B/cgview_xml_builder.pl -sequence cgview/$m.fasta -genes cgview/${m}_contigs.tab cgview/${m}_cds_cog.tab cgview/${m}_rna.tab -gc_content T -gc_skew T -size large -legend T -title "$m ($n contigs, concatenated by length)" -output cgview/$m.xml > cgview/$m.xml.log 2>&1
  for f in png svg; do $B/cgview -i cgview/$m.xml -o ${m}_cgview_map.$f -f $f > cgview/$m.$f.log 2>&1; done
  [ -s ${m}_cgview_map.png ] && echo "CGVDONE $m" || echo "CGVFAIL $m"
done
echo EXTRA_MAPS_ALL_DONE
