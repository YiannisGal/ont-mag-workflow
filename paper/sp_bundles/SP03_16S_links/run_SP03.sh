#!/usr/bin/env bash
# run_SP03.sh -- SP03: 16S links between genomes (29 MAGs + 5 isolates) and amplicon OTUs, with depth profiles.
# Spec: 99_thesis_audit/prompts/SP03_16S_links_to_amplicon_depth_profiles.md
#   bash run_SP03.sh [stage ...]    stages: collect extract search profiles   (none = all)
# User decisions 2026-09-29:
#   - Illumina tables = FULL normalised (otutab_norm.txt + otus.fa), not the abundance-trimmed qplot subset.
#   - Archaeal samples YGnnAR take the depth of bacterial sample YGnnBA (same YG number; no archaeal metadata file exists).
#   - ONT = Etoliko_138_2/ONT/final_output (cluster_tab_trimmed.txt, rarefied to 1,000; cluster_seq_trimmed.fa).
set -euo pipefail
B=${DOCS_DIR}
WS2=$B/Metagenome_analysis_PhD/workspace_2
OUT=$WS2/99_thesis_audit/SP03_16S_links
cd "$OUT"
exec > >(tee -a "$OUT/COMMANDS.log") 2>&1
export B WS2 OUT
ENVBIN=$OUT/env/.pixi/envs/default/bin          # vsearch 2.32.0, cutadapt 5.2, barrnap 1.10.5
export PATH="$ENVBIN:$PATH"
PY=$WS2/envs/smetana/bin/python
CPU=${CPU:-16}; export CPU

run(){ local rc; echo "CMD [$(date -Is)] cwd=$(pwd) :: $*"; set +e; "$@"; rc=$?; set -e
       echo "EXIT [$(date -Is)] rc=$rc :: $1"; return $rc; }
runlog(){ printf '| %s | SP03_%s | %s |\n' "$(date -Is)" "$1" "$2" >> "$WS2/RUNLOG.md"; }

st_collect(){   # barrnap on every genome; one FASTA of all 16S copies
  mkdir -p rrna
  run $PY scripts/s1_collect.py && runlog collect "PASS: genome 16S copies collected" || { runlog collect "**FAIL: collect**"; return 1; }
}
st_extract(){   # template-only regions: primers removed from BOTH genome copies and OTUs (the OTU FASTAs still carry primers)
  mkdir -p regions
  V34="CCTACGGGNGGCWGCAG;required...GGATTAGATACCCBDGTAGTC;required"    # 341F / rc(805R)
  V45="CAGCMGCCGCGGTAA;required...AGGAATTGGCGGGGGAGCAC;required"      # Arc519F / rc(Arc915R)
  # genome copies and the control: BOTH primers required (cutadapt -g A...B otherwise accepts one), else discarded
  run cutadapt -g "$V34" --revcomp -e 0.15 --discard-untrimmed -o regions/V3V4.fa rrna/genome_16S.fa > regions/V3V4.cutadapt.log
  run cutadapt -g "$V45" --revcomp -e 0.15 --discard-untrimmed -o regions/V4V5.fa rrna/genome_16S.fa > regions/V4V5.cutadapt.log
  run cutadapt -g "$V45" --revcomp -e 0.15 --discard-untrimmed -o regions/control_V4V5.fa rrna/control_Arc_Otu22.fa > regions/control_V4V5.cutadapt.log
  # OTUs: primers removed where present; OTUs without them are kept as they are (not discarded)
  run cutadapt -g "${V34//;required/;optional}" -e 0.15 -o regions/otus_bac_template.fa "$B/Etoliko_138_2/Bacteria/out_3/otus.fa" > regions/otus_bac.cutadapt.log
  run cutadapt -g "${V45//;required/;optional}" -e 0.15 -o regions/otus_arc_template.fa "$B/Etoliko_138_2/Archaea/out/otus.fa" > regions/otus_arc.cutadapt.log
  runlog extract "PASS: primer-free template regions (genomes, control, OTUs)"
}
st_search(){
  mkdir -p hits
  F="query+target+id+alnlen+mism+opens+qcov+tcov+ql+tl"
  BAC=regions/otus_bac_template.fa; ARC=regions/otus_arc_template.fa
  ONT=$B/Etoliko_138_2/ONT/final_output/cluster_seq_trimmed.fa
  S="--id 0.90 --strand both --maxaccepts 10 --output_no_hits --userfields $F --threads $CPU --quiet"
  run vsearch --usearch_global regions/V3V4.fa --db "$BAC" --query_cov 0.9 $S --userout hits/V3V4_vs_bac.tsv
  run vsearch --usearch_global regions/V4V5.fa --db "$ARC" --query_cov 0.9 $S --userout hits/V4V5_vs_arc.tsv
  # full gene vs ONT clusters: query longer than cluster, so no query_cov filter; OTU (target) coverage recorded
  run vsearch --usearch_global rrna/genome_16S.fa --db "$ONT" $S --userout hits/full_vs_ont.tsv
  # positive control: Arc_Otu22 full-length gene (and its primer-free Arc519F/915R region) vs archaeal OTUs
  run vsearch --usearch_global rrna/control_Arc_Otu22.fa --db "$ARC" $S --userout hits/control_full_vs_arc.tsv
  run vsearch --usearch_global regions/control_V4V5.fa --db "$ARC" --query_cov 0.9 $S --userout hits/control_V4V5_vs_arc.tsv
  runlog search "PASS: vsearch searches"
}
st_profiles(){
  run $PY scripts/s3_links_profiles.py && runlog profiles "PASS: links, depth profiles, Fig 6.2 data, QC" \
    || { runlog profiles "**FAIL: links/profiles or QC**"; return 1; }
}
STAGES=${*:-collect extract search profiles}
echo "== run_SP03.sh $(date -Is) stages: $STAGES"
for s in $STAGES; do echo "== stage $s"; "st_$s"; done
echo "== done $(date -Is)"
