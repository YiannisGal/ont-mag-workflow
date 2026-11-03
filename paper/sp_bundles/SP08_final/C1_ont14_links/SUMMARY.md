# C1 — genome 16S → ONT14 links (14 depths)
Outputs: `C1_ont14_links.tsv` (per genome, with SP03 comparison), `C1_all_hits_ont14.tsv` (every copy × hit), `C1_depth_profiles_ont14.tsv`, `C1_Fig3_50_input.tsv` (genome × dataset × depth; Illumina from SP03 unchanged + ONT14), `QC.txt`, `genome_md5_check.tsv` (34/34 genomes match SP03 genomes.tsv), `hits/`.
ONT14 = `Etoliko_new/Results_new/final_output/cluster_seq_trimmed.fa` + `cluster_tab_norm.txt` (NanoCLUST, 882 clusters, **68 libraries, 14 depths**; barcode29 and barcode42 of Metadata_minion are absent). The top-level `Etoliko_new/rep_seqs.fasta`/`rep_table.tsv` (15,467 features, 70 libraries) is a different clustering and was not used.
vsearch 2.32.0 `--usearch_global --id 0.90 --strand both --maxaccepts 10 --output_no_hits`. Classes as SP03; coverage = alignment / shorter sequence. Best link per genome = full-length copy (≥ 1,200 bp; SP03's `|full` header tag is unreliable — 74-bp fragments carry it) → strongest class → identity. Plastid clusters excluded.

## QC — all PASS
See QC.txt: Arc_Otu22 no hit; 156/156 self-hits at 100%; 68 libraries/14 depths; 0–5 cm top clusters are Bacillaceae, no methanogen link.

## Changes vs SP03 (SP03 used 26 libraries / 9 depths; cluster IDs are not comparable between runs)
- Exact links lost: YGCA07, YGCA10 (99.7% → 98.2%, now same-genus; peak 60–65 → 45–50 cm), ETO_LIQ_MAG_012 (99.8% → 97.8%).
- Exact link gained: YGCA01_C (C365, 99.9%; peak 60–65 cm).
- New link: ETO_LIQ_MAG_025 (98.0%, peak 25–30 cm). Lost: ETO_LIQ_MAG_001, YGCA01_Tp (now < 94.5%).
- ETO_LIQ_MAG_004: same-genus → related, peak 10–15 → 45–50 cm. ETO_LIQ_MAG_009 same-genus → related.
- Caveat: SP03's best link was taken at the highest identity of any copy, so some SP03 classes may rest on partial copies; ETO_LIQ_MAG_024 has only a 221-bp copy (200-bp alignment).
A ≥ 99% match supports same genus, possibly species; it does not show the genome population is the amplicon population.
