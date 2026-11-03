#!/bin/bash
# B4 (revised 2026-10-11, user decision): instead of re-running the published binning (identical inputs and parameters;
# the first re-run also failed at DAS_Tool on MetaBAT2 header fields), classify the EXISTING workspace_2 bins:
# all raw MetaBAT2 / VAMB / SemiBin2 bins of BOTH libraries (SED and LIQ) + the DAS_Tool union bins (06_binqc/union_bins).
# Then step 3 (contig classification -> read recruitment -> Flye --meta reassembly -> binning) for each library
# that has no Chloroflexota/Desulfobacterota bin at MIMAG medium.
PART=B4_unbinned_phyla; source $(dirname $0)/../../common/lib.sh
set -uo pipefail
T=24; export GTDBTK_DATA_PATH=${SHARED_DB_DIR}/GTDB/release232/
CM2DB=${DB_DIR}/CheckM2_database/uniref100.KO.1.dmnd; GUNCDB=${DB_DIR}/gunc_db/gunc_db_progenomes2.1.dmnd
FLYE=${CONDA_ENVS_SHARED}/NB_Genome/bin/flye   # /usr/local/bin/flye 2.9.4 is broken (import error); used read-only
SP05=$WS2/99_thesis_audit/SP05_dsr_mcr
ver flye "$($FLYE --version 2>&1 | tail -1)" "NB_Genome env (read-only); --meta --nano-hq"
ver gtdbtk_note "classify_wf without ANI-skip flag (2.7.2 has no --skip_ani_screen), as 09_taxonomy.sh"
binqc() {  # binqc <dir_with_fa> <outprefix>
  local d=$1 o=$2
  [ -s $o/checkm2/quality_report.tsv ] || inenv checkm2 checkm2 predict --threads $T --input $d -x fa --output-directory $o/checkm2 --database_path $CM2DB --force
  ls $o/gunc/GUNC*maxCSS*.tsv >/dev/null 2>&1 || { mkdir -p $o/gunc; inenv gunc_env gunc run --input_dir $d --file_suffix .fa --db_file $GUNCDB --out_dir $o/gunc --threads $T --temp_dir $TMPDIR; }
  ls $o/gtdbtk/gtdbtk.*.summary.tsv >/dev/null 2>&1 || { mkdir -p $S/pplacer_scratch; inenv gtdbtk_env gtdbtk classify_wf --genome_dir $d --extension fa --out_dir $o/gtdbtk --cpus $T --pplacer_cpus 8 --scratch_dir $S/pplacer_scratch --tmpdir $TMPDIR; }
}
c2b() {  # c2b <dir> <ext> -> contig<TAB>bin (first header token only; MetaBAT2 adds depth fields to headers)
  for f in $1/*.$2; do b=$(basename $f .$2); grep '>' $f | sed 's/^>//; s/[[:space:]].*//' | awk -v b=$b '{print $1"\t"b}'; done
}
bin3() {  # bin3 <asm> <depthfile> <bam1> <bam2> <outdir>
  local a=$1 dp=$2 b1=$3 b2=$4 o=$5; mkdir -p $o/c2b
  [ -d $o/metabat ] || inenv binning_env metabat2 -i $a -a $dp -o $o/metabat/bin -m 2000 -s 200000 -t $T --seed 1
  [ -d $o/vamb ] || inenv vamb_env vamb bin default --outdir $o/vamb --fasta $a --bamfiles $b1 $b2 -m 2000 --minfasta 200000 -p $T || echo "vamb failed"
  [ -d $o/semibin ] || inenv SemiBin SemiBin2 single_easy_bin -i $a -b $b1 $b2 -o $o/semibin --self-supervised --threads $T --minfasta-kbs 200 --sequencing-type=long_read || echo "semibin failed"
  local il=() ll=()
  ls $o/metabat/*.fa >/dev/null 2>&1 && { c2b $o/metabat fa > $o/c2b/metabat.tsv; il+=($o/c2b/metabat.tsv); ll+=(metabat); }
  ls $o/vamb/bins/*.fna >/dev/null 2>&1 && { c2b $o/vamb/bins fna > $o/c2b/vamb.tsv; il+=($o/c2b/vamb.tsv); ll+=(vamb); }
  if ls $o/semibin/output_bins/*.fa.gz >/dev/null 2>&1; then mkdir -p $o/semibin/unz; for f in $o/semibin/output_bins/*.fa.gz; do zcat $f > $o/semibin/unz/$(basename $f .gz); done
     c2b $o/semibin/unz fa > $o/c2b/semibin.tsv; il+=($o/c2b/semibin.tsv); ll+=(semibin); fi
  [ ${#il[@]} -ge 1 ] || { echo "no binner output for $o"; return 1; }
  [ -d $o/das_DASTool_bins ] || inenv dastool DAS_Tool -i $(IFS=,; echo "${il[*]}") -l $(IFS=,; echo "${ll[*]}") -c $a -o $o/das --search_engine diamond --score_threshold 0.3 --write_bin_evals --write_bins -t $T
}
# ---------------- A. existing bins, both libraries
E=$P/existing; mkdir -p $E/bins
for b in metabat vamb semibin; do for f in $WS2/03_binning/$b/SED_*.fa $WS2/03_binning/$b/LIQ_*.fa; do ln -sf $f $E/bins/$(basename $f); done; done
for f in $WS2/06_binqc/union_bins/*.fa; do ln -sf $f $E/bins/$(basename $f); done
( cd $E/bins && md5sum *.fa > ../bins_md5.tsv ); echo "existing bins: $(ls $E/bins/*.fa | wc -l)" > $E/n_bins.txt
binqc $E/bins $E
run $ENVBIN/python scripts/b4_eval.py existing || true
# ---------------- B. step 3 per library without a qualifying bin
step3() {  # step3 <SED|LIQ>
  local s=$1 d=$P/step3_$1; [ $s = SED ] && d=$P/step3      # SED classification + recruitment already done in step3/
  local ASM BAM; ASM=$WS2/04_polish/${s}_asm_final.fasta; BAM=$WS2/02_mapping/${s}_map_self.sorted.bam
  local RD; [ $s = SED ] && RD=$WS2/00_input/sediment_reads_clean.fastq.gz || RD=$WS2/00_input/liquid_reads.fastq.gz
  mkdir -p $d; cd $d
  [ -s prot.faa ] || { [ -s sed.faa ] && ln -sf sed.faa prot.faa || inenv prokka_env prodigal -p meta -i $ASM -a prot.faa -o /dev/null -q; }
  [ -s vs_trembl.tsv ] || { [ -s sed_vs_trembl.tsv ] && ln -sf sed_vs_trembl.tsv vs_trembl.tsv || run $ENVBIN/diamond blastp -q prot.faa -d $SP05/lca/uniprot_trembl.dmnd --fast --evalue 1e-10 -k 10 --threads $T -b 6 -c 1 --tmpdir $TMPDIR --quiet -f 6 qseqid sseqid pident length evalue bitscore stitle -o vs_trembl.tsv; }
  [ -s target_contigs.tsv ] || run $ENVBIN/python $P/scripts/b4_contig_tax.py $SP05/lca/taxdump prot.faa vs_trembl.tsv
  cut -f1 target_contigs.tsv | tail -n +2 > target_contigs.txt
  [ -s target_contigs.txt ] || { echo "$s: no target contigs" > STOP.txt; cd $P; return 0; }
  [ -s recruited.fastq.gz ] || { [ -s recruited_SED.fastq.gz ] && ln -sf recruited_SED.fastq.gz recruited.fastq.gz || {
       run $ENVBIN/samtools view -@ 8 -b -F 0x900 -q 20 $BAM $(cat target_contigs.txt) -o rec.bam; run $ENVBIN/samtools fastq -@ 8 rec.bam | gzip > recruited.fastq.gz; }; }
  [ -s flye/assembly.fasta ] || run $FLYE --meta --nano-hq recruited.fastq.gz -o flye -t $T
  [ -s flye/assembly.fasta ] || { echo "$s: flye produced no assembly" > STOP.txt; cd $P; return 1; }
  for r in SED:$WS2/00_input/sediment_reads_clean.fastq.gz LIQ:$WS2/00_input/liquid_reads.fastq.gz; do
    k=${r%%:*}; f=${r#*:}
    [ -s map_$k.bam.bai ] || { run $ENVBIN/minimap2 -ax map-ont -t $T --secondary=no flye/assembly.fasta $f | $ENVBIN/samtools sort -@ 8 -m 2G -o map_$k.bam; run $ENVBIN/samtools index map_$k.bam; }
  done
  [ -s depth.txt ] || inenv binning_env jgi_summarize_bam_contig_depths --outputDepth depth.txt map_$s.bam map_$([ $s = SED ] && echo LIQ || echo SED).bam
  local own=map_$s.bam other=map_$([ $s = SED ] && echo LIQ || echo SED).bam
  cd $P; rm -rf $d/bins
  bin3 $d/flye/assembly.fasta $d/depth.txt $d/$own $d/$other $d/bins
  [ -d $d/bins/das_DASTool_bins ] && binqc $d/bins/das_DASTool_bins $d/bins
  run $ENVBIN/python scripts/b4_eval.py $(basename $d)/bins || true
}
for s in SED LIQ; do
  if grep -q "^STEP2_HIT_$s" existing/eval_flag.txt 2>/dev/null; then echo "$s: qualifying bin among existing bins; step 3 not needed" > step3_$s.skip; continue; fi
  step3 $s
done
