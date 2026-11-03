#!/bin/bash
PART=B5_mag022_dsr; source $(dirname $0)/../../common/lib.sh
set -uo pipefail; rm -f inputs_md5.tsv; printf 'item\tversion\tnote\n' > versions.tsv; mkdir -p work
ASM=$WS2/04_polish/LIQ_asm_final.fasta; BAM=$WS2/02_mapping/LIQ_map_self.sorted.bam; MAG=$WS2/07_mags/bins/LIQ_das_metabat_022.fa
READS=$WS2/00_input/liquid_reads.fastq.gz
md5in $ASM $MAG $BAM $READS
ver samtools "$($ENVBIN/samtools --version | head -1)" "reads BAM (read-only reuse): minimap2 2.28 -ax map-ont --secondary=no, LIQ reads -> LIQ_asm_final (02_mapping, step 05)"
ver checkm2 "$($CE/checkm2/bin/checkm2 --version 2>&1 | tail -1)" "DB ${DB_DIR}/CheckM2_database/uniref100.KO.1.dmnd"
ver gunc "$($CE/gunc_env/bin/gunc --version 2>&1 | tail -1)" "DB gunc_db_progenomes2.1.dmnd"
ver gtdbtk "$($CE/gtdbtk_env/bin/gtdbtk --version 2>&1 | grep -m1 -i version)" "identify only (bac120 marker counts); GTDBTK_DATA_PATH release232"
grep '>' $MAG | tr -d '>' | cut -d' ' -f1 > work/mag022_contigs.txt
# 1. physical linkage: every read with an alignment (MAPQ>=20) on LIQ_contig_19299 or a MAG_022 contig, with its SA partners
run $ENVBIN/samtools view -@ 8 -q 20 $BAM LIQ_contig_19299 $(cat work/mag022_contigs.txt) > work/region_alns.sam
run $ENVBIN/python scripts/b5_link.py
# 2. forced inclusion
mkdir -p work/genomes; cp $MAG work/genomes/MAG022.fa
$ENVBIN/seqkit grep -p LIQ_contig_19299 $ASM > work/contig_19299.fa
cat $MAG work/contig_19299.fa > work/genomes/MAG022_plus_19299.fa; cp work/genomes/MAG022_plus_19299.fa MAG022_plus_19299.fa
inenv checkm2 checkm2 predict --threads 16 --input work/genomes -x fa --output-directory work/checkm2 --database_path ${DB_DIR}/CheckM2_database/uniref100.KO.1.dmnd --force
mkdir -p work/gunc; inenv gunc_env gunc run --input_dir work/genomes --file_suffix .fa --db_file ${DB_DIR}/gunc_db/gunc_db_progenomes2.1.dmnd --out_dir work/gunc --threads 16 --temp_dir $TMPDIR
GTDBTK_DATA_PATH=${SHARED_DB_DIR}/GTDB/release232/ inenv gtdbtk_env gtdbtk identify --genome_dir work/genomes -x fa --out_dir work/gtdbtk --cpus 16
# 3. compatibility (GC, TNF, coverage) for reference
run $ENVBIN/samtools coverage $BAM -o work/contig_coverage.tsv
run $ENVBIN/python scripts/b5_report.py
