export GTDBTK_DATA_PATH=${SHARED_DB_DIR}/GTDB/release232/
for m in $(cut -f1 sets.tsv | tail -n+2); do
  d=bac120; grep -q "ar53" <(awk -v m=$m '$1==m' sets.tsv) && d=ar53
  case $m in ETO_LIQ_MAG_003|ETO_LIQ_MAG_007|ETO_LIQ_MAG_010) d=ar53;; esac
  [ -s $m/align/align/gtdbtk.$d.user_msa.fasta.gz ] || { /usr/local/bin/gtdbtk identify --genome_dir $m/genomes --extension fa --out_dir $m/align --cpus 8 --tmpdir /tmp/sp09tmp > $m/identify.log 2>&1; /usr/local/bin/gtdbtk align --identify_dir $m/align --out_dir $m/align --cpus 8 --tmpdir /tmp/sp09tmp > $m/align.log 2>&1; }
  zcat $m/align/align/gtdbtk.$d.user_msa.fasta.gz > $m/$m.$d.msa.faa
  og=ref_$(head -1 $m/outgroup.txt)
  conda run -n iqtree3 iqtree3 -s $m/$m.$d.msa.faa -m MFP -B 1000 -alrt 1000 -T 6 -o $og --prefix $m/$m -redo > $m/iqtree.log 2>&1 || echo "IQTREE FAIL $m"
  echo "DONE $m $d"
done
echo A6_ALL_DONE
