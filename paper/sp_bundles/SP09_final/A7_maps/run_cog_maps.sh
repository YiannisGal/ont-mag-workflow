# COGclassifier (Proksee's COG tool) -> pyCirclize + CGView maps + Proksee tables, all 9 candidates
B=${CONDA_ENVS_USER}/sp09_cgview/bin
for x in ETO_LIQ_MAG_001:Y ETO_LIQ_MAG_007:Y ETO_SED_MAG_001:N ETO_LIQ_MAG_003:N ETO_LIQ_MAG_009:N ETO_LIQ_MAG_010:N ETO_LIQ_MAG_005:N ETO_LIQ_MAG_008:N ETO_LIQ_MAG_011:N; do
  m=${x%:*}; c=${x#*:}; n=$(grep -c ">" bakta/$m/$m.fna)
  [ -s cog/$m/cog_classify.tsv ] || conda run -n sp09_cog COGclassifier -i bakta/$m/$m.faa -o cog/$m -d cog/db -t 8 -q > cog/$m.log 2>&1 || { echo "COGFAIL $m"; continue; }
  # 2-column table: locus_tag, COG letter (first letter if several)
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{print $h["QUERY_ID"]"\t"substr($h["COG_LETTER"],1,1)}' cog/$m/cog_classify.tsv > cog/$m.cog2.tsv
  echo "COGDONE $m $(wc -l < cog/$m.cog2.tsv) of $(grep -c '>' bakta/$m/$m.faa)"
  conda run -n sp09_bakta python -I scripts/plot_map.py $m bakta/$m/$m.gbff ${m}_circular_map x $c cog/$m.cog2.tsv > maplog_$m.txt 2>&1 && echo "MAPDONE $m $(tail -n1 maplog_$m.txt)" || echo "MAPFAIL $m"
  cp bakta/$m/$m.gbff proksee_upload/
  conda run -n sp09_bakta python -I scripts/cgview_inputs.py $m bakta/$m/$m.gbff cog/$m.cog2.tsv cgview proksee_upload > cgview/$m.prep.log 2>&1
  t="$m ($( [ $c = Y ] && echo 'single circular contig' || echo "$n contigs, concatenated by length"))"
  $B/cgview_xml_builder.pl -sequence cgview/$m.fasta -genes cgview/${m}_contigs.tab cgview/${m}_cds_cog.tab cgview/${m}_rna.tab -gc_content T -gc_skew T -size large -legend T -title "$t" -output cgview/$m.xml > cgview/$m.xml.log 2>&1
  for f in png svg; do $B/cgview -i cgview/$m.xml -o ${m}_cgview_map.$f -f $f > cgview/$m.$f.log 2>&1; done
  [ -s ${m}_cgview_map.png ] && echo "CGVDONE $m" || echo "CGVFAIL $m"
done
echo COGMAPS_ALL_DONE
