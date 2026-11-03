export GTDBTK_DATA_PATH=${SHARED_DB_DIR}/GTDB/release232/
[ -s align/align/gtdbtk.bac120.user_msa.fasta.gz ] || { /usr/local/bin/gtdbtk identify --genome_dir genomes --extension fa --out_dir align --cpus 8 --tmpdir /tmp/sp09tmp > identify.log 2>&1; /usr/local/bin/gtdbtk align --identify_dir align --out_dir align --cpus 8 --tmpdir /tmp/sp09tmp > align.log 2>&1; }
for d in bac120 ar53; do
  zcat align/align/gtdbtk.$d.user_msa.fasta.gz > combined.$d.msa.faa || { echo "A10FAIL align $d"; continue; }
  conda run -n iqtree3 iqtree3 -s combined.$d.msa.faa -m LG+F+G4 -B 1000 -alrt 1000 -T 6 --prefix combined_$d -redo > iqtree_$d.log 2>&1 && echo "A10DONE $d" || echo "A10FAIL iqtree $d"
done
echo A10_ALL_DONE
