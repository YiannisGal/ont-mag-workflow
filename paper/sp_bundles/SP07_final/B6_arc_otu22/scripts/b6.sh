#!/bin/bash
PART=B6_arc_otu22; source $(dirname $0)/../../common/lib.sh
set -uo pipefail; rm -f inputs_md5.tsv versions.tsv; printf 'item\tversion\tnote\n' > versions.tsv
V=$ENVBIN/vsearch; CA=$ENVBIN/cutadapt; SK=$ENVBIN/seqkit
GENE=$WS2/99_thesis_audit/inputs/Arc_Otu22_fulllength_16S.fasta; MS=$S/inputs/Arc_Otu22_Otu674_Illumina_amplicons.fasta; OTU=$AMP/Archaea/out/otus.fa
md5in $GENE $MS $OTU
ver vsearch "$($V --version 2>&1 | head -1)"; ver cutadapt "$($CA --version)"; ver seqkit "$($SK version)"
mkdir -p work; cd work
run $SK grep -r -p '^Otu22(;|$)' $OTU | $SK replace -p '^(\S+).*' -r 'Otu22_server' > server.fa
run $SK grep -p Arc_Otu22_Illumina $MS | $SK replace -p '.*' -r 'Arc_Otu22_ms' > ms.fa
cp $GENE gene.fa
[ "$(grep -c '>' server.fa)" = 1 ] || { echo "QC FAIL: Otu22 not uniquely found in otus.fa"; exit 1; }
run $V --usearch_global server.fa --db ms.fa --id 0.5 --iddef 2 --userout server_vs_ms.tsv --userfields query+target+id+alnlen+mism+gaps+ql+tl --quiet
for a in server ms; do
  # linked adapter: Arc519F at 5', reverse complement of Arc915R at 3'; both required
  run $CA -g 'CAGCMGCCGCGGTAA...AGGAATTGGCGGGGGAGCAC' -e 0.15 --discard-untrimmed -o ${a}_pf.fa $a.fa > cutadapt_$a.log 2>&1
done
# positive control: gene vs itself
run $V --usearch_global gene.fa --db gene.fa --id 0.5 --iddef 2 --userout ctrl.tsv --userfields query+target+id --quiet
for a in server ms; do for d in 0 1 2; do
  run $V --usearch_global ${a}_pf.fa --db gene.fa --id 0.0 --iddef $d --maxaccepts 0 --maxrejects 0 \
     --userout ${a}_iddef$d.tsv --userfields query+target+id+alnlen+ids+mism+gaps+ql+tl+tlo+thi --alnout ${a}_iddef$d.aln --quiet
done; done
cd ..; $ENVBIN/python scripts/b6_report.py
