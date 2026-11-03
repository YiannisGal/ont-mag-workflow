#!/bin/bash
# B4: unbinned abundant phyla. Step 1 = multi-sample (SED self + LIQ cross coverage) re-binning of the SED assembly.
# The published workspace_2 binning (scripts/06_bin.sh, 07_consensus.sh, 08d) already used exactly these two coverage
# columns, so step 1 is run with the published parameters and doubles as the reproduction QC.
PART=B4_unbinned_phyla; source $(dirname $0)/../../common/lib.sh
set -uo pipefail; printf 'item\tversion\tnote\n' > versions.tsv; rm -f inputs_md5.tsv
T=24; ASM=$WS2/04_polish/SED_asm_final.fasta; M=$WS2/02_mapping
BS=$M/SED_map_self.sorted.bam; BC=$M/SED_map_cross.sorted.bam; DEPTH=$M/SED_depth_multi.txt
export GTDBTK_DATA_PATH=${SHARED_DB_DIR}/GTDB/release232/
CM2DB=${DB_DIR}/CheckM2_database/uniref100.KO.1.dmnd; GUNCDB=${DB_DIR}/gunc_db/gunc_db_progenomes2.1.dmnd
md5in $ASM $DEPTH $BS $BC $WS2/00_input/sediment_reads_clean.fastq.gz $WS2/00_input/liquid_reads.fastq.gz
ver metabat2 "$($CE/binning_env/bin/metabat2 --help 2>&1 | grep -m1 -i version)" "-m 2000 -s 200000 --seed 1 (as 06_bin.sh)"
ver vamb "$($CE/vamb_env/bin/vamb --version 2>&1 | tail -1)" "bin default -m 2000 --minfasta 200000"
ver SemiBin2 "$($CE/SemiBin/bin/SemiBin2 --version 2>&1 | tail -1)" "single_easy_bin --self-supervised long_read, both BAMs"
ver DAS_Tool "$($CE/dastool/bin/DAS_Tool --version 2>&1 | tail -1)" "diamond, score_threshold 0.3 (adopted at 08d)"
ver checkm2 "$($CE/checkm2/bin/checkm2 --version 2>&1 | tail -1)" "$CM2DB"
ver gunc "$($CE/gunc_env/bin/gunc --version 2>&1 | tail -1)" "$GUNCDB"
ver gtdbtk "$($CE/gtdbtk_env/bin/gtdbtk --version 2>&1 | grep -m1 -i version)" "GTDB R232 ($GTDBTK_DATA_PATH)"
binqc() {  # binqc <dir_with_fa> <outprefix>
  local d=$1 o=$2
  [ -s $o/checkm2/quality_report.tsv ] || inenv checkm2 checkm2 predict --threads $T --input $d -x fa --output-directory $o/checkm2 --database_path $CM2DB --force
  ls $o/gunc/GUNC*maxCSS*.tsv >/dev/null 2>&1 || { mkdir -p $o/gunc; inenv gunc_env gunc run --input_dir $d --file_suffix .fa --db_file $GUNCDB --out_dir $o/gunc --threads $T --temp_dir $TMPDIR; }
  ls $o/gtdbtk/gtdbtk.*.summary.tsv >/dev/null 2>&1 || { mkdir -p $TMPDIR/pp; inenv gtdbtk_env gtdbtk classify_wf --genome_dir $d --extension fa --out_dir $o/gtdbtk --cpus $T --pplacer_cpus 8 --scratch_dir $TMPDIR/pp --tmpdir $TMPDIR --skip_ani_screen; }
}
bin3() {  # bin3 <asm> <depthfile> <bam1> <bam2> <outdir>
  local a=$1 dp=$2 b1=$3 b2=$4 o=$5; mkdir -p $o/c2b
  [ -d $o/metabat ] || inenv binning_env metabat2 -i $a -a $dp -o $o/metabat/bin -m 2000 -s 200000 -t $T --seed 1
  [ -d $o/vamb ] || inenv vamb_env vamb bin default --outdir $o/vamb --fasta $a --bamfiles $b1 $b2 -m 2000 --minfasta 200000 -p $T || echo "vamb failed"
  [ -d $o/semibin ] || inenv SemiBin SemiBin2 single_easy_bin -i $a -b $b1 $b2 -o $o/semibin --self-supervised --threads $T --minfasta-kbs 200 --sequencing-type=long_read || echo "semibin failed"
  local il=() ll=()
  for x in "metabat:$o/metabat:fa" "vamb:$o/vamb/bins:fna" "semibin:$o/semibin/output_bins:fa.gz"; do
    IFS=: read n dd ext <<< "$x"; ls $dd/*.$ext >/dev/null 2>&1 || continue
    [ "$ext" = fa.gz ] && { mkdir -p $o/semibin/unz; for f in $dd/*.fa.gz; do zcat $f > $o/semibin/unz/$(basename $f .gz); done; dd=$o/semibin/unz; ext=fa; }
    inenv dastool Fasta_to_Contig2Bin.sh -i $dd -e $ext > $o/c2b/$n.tsv; il+=($o/c2b/$n.tsv); ll+=($n)
  done
  [ -d $o/das_DASTool_bins ] || inenv dastool DAS_Tool -i $(IFS=,; echo "${il[*]}") -l $(IFS=,; echo "${ll[*]}") -c $a -o $o/das --search_engine diamond --score_threshold 0.3 --write_bin_evals --write_bins -t $T
}
# ---- step 1 (+ reproduction QC): published inputs, published parameters
bin3 $ASM $DEPTH $BS $BC $P/step1
binqc $P/step1/das_DASTool_bins $P/step1
run $ENVBIN/python scripts/b4_eval.py step1 || true
if grep -q "^STEP2_HIT" step1/eval_flag.txt 2>/dev/null; then echo "step 2 found Chloroflexota/Desulfobacterota MIMAG-medium bins; stopping before step 3"; exit 0; fi
# ---- step 3: classify SED contigs (CAT-like: DIAMOND top-10 LCA per protein vs UniProt, majority vote per contig), recruit reads, reassemble
mkdir -p step3; cd step3
[ -s sed.faa ] || inenv prokka_env prodigal -p meta -i $ASM -a sed.faa -o /dev/null -q
SP05=$WS2/99_thesis_audit/SP05_dsr_mcr
[ -s sed_vs_trembl.tsv ] || run $ENVBIN/diamond blastp -q sed.faa -d $SP05/lca/uniprot_trembl.dmnd --fast --evalue 1e-10 -k 10 --threads $T -b 6 -c 1 --tmpdir $TMPDIR --quiet -f 6 qseqid sseqid pident length evalue bitscore stitle -o sed_vs_trembl.tsv
run $ENVBIN/python ../scripts/b4_contig_tax.py $SP05/lca/taxdump
cut -f1 target_contigs.tsv | tail -n +2 > target_contigs.txt
if [ ! -s target_contigs.txt ]; then echo "no Chloroflexota/Desulfobacterota contigs" > ../step3_STOP.txt; exit 0; fi
[ -s recruited_SED.fastq.gz ] || { run $ENVBIN/samtools view -@ 8 -b -F 0x900 -q 20 $BS $(cat target_contigs.txt) -o rec.bam; run $ENVBIN/samtools fastq -@ 8 rec.bam | gzip > recruited_SED.fastq.gz; }
[ -s flye/assembly.fasta ] || run /usr/local/bin/flye --meta --nano-hq recruited_SED.fastq.gz -o flye -t $T
ASM2=$P/step3/flye/assembly.fasta
for r in SED:$WS2/00_input/sediment_reads_clean.fastq.gz LIQ:$WS2/00_input/liquid_reads.fastq.gz; do
  k=${r%%:*}; f=${r#*:}
  [ -s map_$k.bam.bai ] || { run $ENVBIN/minimap2 -ax map-ont -t $T --secondary=no $ASM2 $f | $ENVBIN/samtools sort -@ 8 -m 2G -o map_$k.bam; run $ENVBIN/samtools index map_$k.bam; }
done
[ -s depth.txt ] || inenv binning_env jgi_summarize_bam_contig_depths --outputDepth depth.txt map_SED.bam map_LIQ.bam
cd $P
bin3 $ASM2 $P/step3/depth.txt $P/step3/map_SED.bam $P/step3/map_LIQ.bam $P/step3/bins
[ -d $P/step3/bins/das_DASTool_bins ] && binqc $P/step3/bins/das_DASTool_bins $P/step3/bins
run $ENVBIN/python scripts/b4_eval.py step3/bins || true
