# ANI screen: 16 unnamed MQ MAGs vs all GTDB r232 genomes (skani index), plus all Etoliko genomes vs each other
skani search -d ${SHARED_DB_DIR}/GTDB/release232/skani/database genomes/*.fa -t 8 -n 5 -o vs_gtdb_r232.tsv > search.log 2>&1 && echo "A12DONE gtdb" || echo "A12FAIL gtdb"
ls ../A10_combined_tree/genomes/ETO_*.fa ../A10_combined_tree/genomes/YG*.fa > etoliko_list.txt
skani triangle -l etoliko_list.txt -t 8 -E -o etoliko_all_vs_all.tsv > triangle.log 2>&1 && echo "A12DONE triangle" || echo "A12FAIL triangle"
echo A12_ALL_DONE
