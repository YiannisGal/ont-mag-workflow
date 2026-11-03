S=${WS2}/99_thesis_audit/SP09_final; X=$S/A13_extra_candidates; cd $X
EZ=${CONDA_ENVS_USER}/sp09_aai/bin; export PATH=$EZ:$PATH
C="ETO_LIQ_MAG_005 ETO_LIQ_MAG_008 ETO_LIQ_MAG_011"
# A3 ANI
for m in $C; do tail -n+2 $S/A2_refs/${m}_refs.tsv | cut -f1 | sed "s#^#$S/A2_refs/genomes/#;s#\$#.fna#" > ani/$m.list
  skani dist -q genomes/$m.fa --rl ani/$m.list -t 6 -n 1000 > ani/${m}_skani.tsv 2>/dev/null; fastANI -q genomes/$m.fa --rl ani/$m.list -t 6 -o ani/${m}_fastani.tsv >/dev/null 2>&1; done
skani triangle genomes/ETO_LIQ_MAG_005.fa genomes/ETO_LIQ_MAG_007.fa genomes/ETO_LIQ_MAG_011.fa -E -o ani/methanofollis_etoliko_triangle.tsv 2>/dev/null
echo "A13DONE ani"
# A4 AAI
for g in $C; do [ -s $S/A4_aai/db/$g.db ] || ezaai extract -i genomes/$g.fa -o $S/A4_aai/db/$g.db -l $g >/dev/null 2>&1; done
for a in $(tail -n+2 $S/A2_refs/ETO_LIQ_MAG_008_refs.tsv | cut -f1); do [ -s $S/A4_aai/db/$a.db ] || ezaai extract -i $S/A2_refs/genomes/$a.fna -o $S/A4_aai/db/$a.db -l $a >/dev/null 2>&1; done
{ for m in $C; do tail -n+2 $S/A2_refs/${m}_refs.tsv | cut -f1 | sed "s/^/$m /"; done
  r=($(tail -n+2 $S/A2_refs/ETO_LIQ_MAG_008_refs.tsv | cut -f1)); for i in ${!r[@]}; do for j in ${!r[@]}; do [ $i -lt $j ] && echo "${r[$i]} ${r[$j]}"; done; done
  echo "ETO_LIQ_MAG_005 ETO_LIQ_MAG_007"; echo "ETO_LIQ_MAG_005 ETO_LIQ_MAG_011"; echo "ETO_LIQ_MAG_007 ETO_LIQ_MAG_011"; } > aai/pairs.txt
cat aai/pairs.txt | xargs -P 12 -n2 sh -c 'o=aai/pairs/$0__$1.tsv; [ -s $o ] || [ -s '$S'/A4_aai/pairs/$0__$1.tsv ] && { [ -s $o ] || cp '$S'/A4_aai/pairs/$0__$1.tsv $o; } || (mkdir -p aai/tmp/$0__$1 && cd aai/tmp/$0__$1 && ezaai calculate -i '$S'/A4_aai/db/$0.db -j '$S'/A4_aai/db/$1.db -o ../../../$o -t 1 >/dev/null 2>&1)'
echo "A13DONE aai $(ls aai/pairs | wc -l)/$(wc -l < aai/pairs.txt)"
# A5 POCP (type species / GTDB rep, same rule as A5)
printf "ETO_LIQ_MAG_005\tGCF_013377755.1\nETO_LIQ_MAG_011\tGCF_013377755.1\n" > pocp/pairs.tsv
awk -F'\t' 'NR>1' $S/A2_refs/ETO_LIQ_MAG_008_refs.tsv | sort -t$'\t' -k7,7 | head -1 | awk -F'\t' '{print "ETO_LIQ_MAG_008\t"$1}' >> pocp/pairs.tsv
cut -f1 pocp/pairs.tsv | sort -u | while read g; do prodigal -q -i genomes/$g.fa -a pocp/faa/$g.faa -o /dev/null; done
cut -f2 pocp/pairs.tsv | sort -u | while read g; do [ -s $S/A5_pocp/faa/$g.faa ] && cp $S/A5_pocp/faa/$g.faa pocp/faa/ || prodigal -q -i $S/A2_refs/genomes/$g.fna -a pocp/faa/$g.faa -o /dev/null; done
while IFS=$'\t' read q s; do for p in "$q $s" "$s $q"; do set -- $p; blastp -query pocp/faa/$1.faa -subject pocp/faa/$2.faa -evalue 1e-5 -max_hsps 1 -max_target_seqs 1000 -outfmt "6 qseqid sseqid pident length qstart qend qlen evalue" > pocp/blast/$1__$2.tsv; done; done < pocp/pairs.tsv
echo "A13DONE pocp"
# A7 Bakta
for m in $C; do conda run -n sp09_bakta bakta --db $S/A7_maps/bakta_db_light_copy --output bakta/$m --prefix $m --locus-tag $(echo ${m#ETO_} | tr -d _) --threads 8 --force --skip-plot genomes/$m.fa > bakta/$m.log 2>&1 || echo "A13FAIL bakta $m"; done
echo "A13DONE bakta"
# A6 trees: Methanofollis (005+007+011) and DIGJ01 (008)
export GTDBTK_DATA_PATH=${SHARED_DB_DIR}/GTDB/release232/
for t in methanofollis digj01; do mkdir -p trees/$t/genomes
  if [ $t = methanofollis ]; then for g in ETO_LIQ_MAG_005 ETO_LIQ_MAG_007 ETO_LIQ_MAG_011; do ln -sf $X/genomes/$g.fa trees/$t/genomes/; done; R=$S/A2_refs/ETO_LIQ_MAG_007_refs.tsv; O=$S/A6_trees/ETO_LIQ_MAG_007/outgroup.txt; d=ar53
  else ln -sf $X/genomes/ETO_LIQ_MAG_008.fa trees/$t/genomes/; R=$S/A2_refs/ETO_LIQ_MAG_008_refs.tsv; O=$S/A2_refs/ETO_LIQ_MAG_008_outgroup.txt; d=bac120; fi
  for a in $(tail -n+2 $R | cut -f1) $(cat $O); do ln -sf $S/A2_refs/genomes/$a.fna trees/$t/genomes/ref_$a.fa; done; cp $O trees/$t/outgroup.txt
  /usr/local/bin/gtdbtk identify --genome_dir trees/$t/genomes --extension fa --out_dir trees/$t/align --cpus 8 --tmpdir /tmp/sp09tmp > trees/$t/identify.log 2>&1
  /usr/local/bin/gtdbtk align --identify_dir trees/$t/align --out_dir trees/$t/align --cpus 8 --tmpdir /tmp/sp09tmp > trees/$t/align.log 2>&1
  zcat trees/$t/align/align/gtdbtk.$d.user_msa.fasta.gz > trees/$t/$t.$d.msa.faa
  conda run -n iqtree3 iqtree3 -s trees/$t/$t.$d.msa.faa -m MFP -B 1000 -alrt 1000 -T 6 -o ref_$(head -1 $O) --prefix trees/$t/$t -redo > trees/$t/iqtree.log 2>&1 && echo "A13DONE tree $t" || echo "A13FAIL tree $t"
done
echo A13_ALL_DONE
