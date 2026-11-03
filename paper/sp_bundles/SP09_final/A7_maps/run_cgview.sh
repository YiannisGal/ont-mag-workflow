# waits for eggNOG, then builds CGView maps and Proksee upload files for the 6 candidates
until [ -s eggnog/all6.emapper.annotations ] || grep -q "EGGFAIL" run_maps.log; do sleep 60; done
[ -s eggnog/all6.emapper.annotations ] || { echo "CGVFAIL eggnog"; exit 1; }
B=${CONDA_ENVS_USER}/sp09_cgview/bin
for x in ETO_LIQ_MAG_001:Y ETO_LIQ_MAG_007:Y ETO_SED_MAG_001:N ETO_LIQ_MAG_003:N ETO_LIQ_MAG_009:N ETO_LIQ_MAG_010:N; do
  m=${x%:*}; c=${x#*:}; p=$(echo ${m#ETO_} | tr -d _)_
  grep "^$p" eggnog/all6.emapper.annotations > eggnog/$m.cgv.annotations
  cp bakta/$m/$m.gbff proksee_upload/
  conda run -n sp09_bakta python -I scripts/cgview_inputs.py $m bakta/$m/$m.gbff eggnog/$m.cgv.annotations cgview proksee_upload > cgview/$m.prep.log 2>&1 || { echo "CGVFAIL $m prep"; continue; }
  t="$m ($( [ $c = Y ] && echo 'single circular contig' || echo "$(grep -c '>' bakta/$m/$m.fna) contigs, concatenated by length"))"
  $B/cgview_xml_builder.pl -sequence cgview/$m.fasta -genes cgview/${m}_contigs.tab cgview/${m}_cds_cog.tab cgview/${m}_rna.tab -gc_content T -gc_skew T -size large -legend T -title "$t" -output cgview/$m.xml > cgview/$m.xml.log 2>&1 || { echo "CGVFAIL $m xml"; continue; }
  for f in png svg; do $B/cgview -i cgview/$m.xml -o ${m}_cgview_map.$f -f $f > cgview/$m.$f.log 2>&1 || echo "CGVFAIL $m $f"; done
  [ -s ${m}_cgview_map.png ] && echo "CGVDONE $m"
done
echo CGV_ALL_DONE
