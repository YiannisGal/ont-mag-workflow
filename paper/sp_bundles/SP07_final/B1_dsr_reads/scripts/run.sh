#!/bin/bash
PART=B1_dsr_reads; source $(dirname $0)/../../common/lib.sh
set -uo pipefail; printf 'item\tversion\tnote\n' > versions.tsv; rm -f inputs_md5.tsv
SP05=$WS2/99_thesis_audit/SP05_dsr_mcr
md5in $S/common/reads/read_ids.tsv $S/common/reads/reads.fa $WS2/34_unbinned/hits/SED_reads.tsv $WS2/31_readlevel/SED_read_hits.tsv $SP05/trees/*.treefile $SP05/trees/*.trim $SP05/trees/*_in.faa
cat $S/common/reads/reads_md5.tsv >> inputs_md5.tsv
for t in diamond mafft epa-ng gappa; do ver $t "$($ENVBIN/$t --version 2>&1 | head -1)"; done
ver UniProt "SP05 lca/uniprot_sprot.dmnd + uniprot_trembl.dmnd (Release 2026_03)" "blastx --frameshift 15 --more-sensitive -e 1e-10 -k 10; top-10 LCA, NCBI taxdump"
ver model "dsr LG+G4+F; mcrA LG+I+G4 (EPA-ng cannot use IQ-TREE's C20 mixture / FreeRate R4; parameters re-estimated on the SP05 tree)"
while [ ! -s $S/common/uniprot/read_LCA.tsv ]; do sleep 60; done
for g in dsrA dsrB; do run $ENVBIN/python $S/common/place.py $g $P > place_$g.out 2>&1 || echo "placement $g failed" >> FAILED.txt; done
