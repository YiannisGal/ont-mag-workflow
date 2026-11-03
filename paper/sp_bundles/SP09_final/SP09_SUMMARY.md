# SP09 summary: taxonomic novelty of the Etoliko MAGs (run 2026-10-21/22)

**Status: final (updated 2026-10-24).** dDDH added from TYGS (GGDC formula 2 / d4; 12 TYGS result sets submitted by the student, retrieved through the TYGS API). Upload `dDDH_package.tar.gz` and `dDDH_package_2.tar.gz` to GGDC and put the results in `dDDH_results/`.

Order run: A9 → A1 → A8 (first pass). On request, A2–A7 were added, then A10 (combined tree of MAGs + culturomics isolates), A11 (species-level panels in the YGCC01 Figure 5 style), A12 (ANI screen of unnamed medium-quality MAGs) and A13 (3 SeqCode-eligible medium-quality MAGs).

## Status by part (2026-10-22 04:20)
| part | status | output |
|---|---|---|
| A9 NCBI submission | ✅ QC passed; adapters trimmed | `$WS2/NCBI_submission/` |
| A1 GTDB-Tk r232 | ✅ QC passed | `A1_gtdbtk/` |
| A2 reference sets | ✅ | `A2_refs/` |
| A3 ANI (skani + FastANI) | ✅ QC passed (controls ≥95%) | `A3_ani/`, `A13_extra_candidates/ani/` |
| A4 AAI (EzAAI) | ✅ | `A4_aai/`, `A13_extra_candidates/aai/` |
| A5 POCP / POCPu | ✅ | `A5_pocp/pocp.tsv`, `A13_extra_candidates/pocp/pocp.tsv` |
| A6 trees | ✅ 8/8 + *Methanofollis* (005/007/011) + DIGJ01; **QC passed** (both controls inside their named species) | `A6_trees/*_tree.pdf/.svg/.nwk`, `A13_extra_candidates/*_tree.*` |
| A7 maps | ✅ 9 candidates × (pyCirclize + CGView) + Proksee upload files; COG from COGclassifier | `A7_maps/` |
| A8 table | ✅ 9 candidates + 2 controls, incl. dDDH | `A8_table/novel_taxa_summary.{tsv,md}` |
| A10 combined tree (Fig. 6 style) | ✅ 29 MAGs + 5 isolates + 269 GTDB refs; bac120 + ar53 trees (IQ-TREE LG+F+G4, UFBoot/SH-aLRT 1000), midpoint-rooted and joined for display | `A10_combined_tree/Etoliko_combined_tree.{pdf,svg,png}` (full names) and `Etoliko_combined_tree_codes.{pdf,svg,png}` (short codes + key; key also in `tree_code_key.tsv`), `combined_*.treefile` |
| A11 species panels (Fig. 5 style) | ✅ 9 panels (tree + ANI, AAI, dDDH, ΔG+C) | `A11_species_panels/*_species_panel.pdf/.svg` |
| A12 MQ ANI screen | ✅ | `A12_MQ_ani_screen/` |
| A13 extra candidates | ✅ | `A13_extra_candidates/` |
| dDDH (TYGS) | ✅ all 9 candidates < 70% (max 28.4–37.8%) | `dDDH_results/tygs/` (raw API downloads), `dDDH_results/tygs_dDDH_all_MAG_pairs.tsv` |

## Candidate verdicts
Criteria: species, ANI < 95% to every reference (GTDB/Parks); genus, AAI below every genus's internal minimum and < ~65% (Konstantinidis et al. 2017), POCPu < 50% (Qin et al. 2014); RED and the tree as supporting evidence; dDDH (TYGS, formula d4) < 70% for species.

**ETO_SED_MAG_001, *Sedimentibacter* sp. nov. (HQ).** RED 0.989, placed in *Sedimentibacter* by topology. Highest ANI 86.3% (*S.* sp032006665), AAI 89.4% to the genus (internal range 67.9–95.5%), POCPu 54.0% to *S. hydroxybenzoicus* (type species). The tree places it sister to *S.* sp032006665 (97/100). **New species supported; genus confirmed.** Meets the SeqCode type-genome checks (51 contigs, 16S present). Sediment breadth 0.9999, the most widespread candidate in the core.

**ETO_LIQ_MAG_001, Cytobacillaceae, novel genus-level lineage (HQ).** A single circular contig of 4.90 Mb. RED 0.912 lies within the GTDB genus interval (0.844–1.044) and outside the family interval, so GTDB-Tk treats it as genus-level novelty. ANI < 80% to all 150 Cytobacillaceae references. AAI to the 25 genera is 67.5–74.4%: none is below 65%, and it is not below the internal minimum of *Robertmurraya* (74.4 vs 73.7), *Robertmurraya_A* (73.5 vs 72.4) or *Cytobacillus* (71.7 vs 68.3). POCPu is 46.0–58.1% (≥ 50% against 21 of 25 genera). The tree places it sister to the uncultured placeholder genus CALKIU01 (100/100), and that pair sister to *Robertmurraya* (100/100), so it lies outside every named genus. **New species supported. The new genus is supported by RED and phylogeny, but the AAI and POCP threshold criteria are NOT met**; describe it as a novel genus-level lineage, not as a definitive new genus. Its complete circular genome makes it the best available type genome for this lineage.

**ETO_LIQ_MAG_003, MTP4 sp. nov. (Methanosarcinaceae; HQ).** RED 0.917. Highest ANI 88.7%, AAI 89.8% (MTP4 range 88.1–94.8%), POCPu 80.1% (GTDB representative; MTP4 has no type material). Basal to MTP4 in the tree (MTP4 + MAG 100/100; internal split 65.7/73). **New species supported.** MTP4 is an unnamed placeholder genus. One complete *mcr* operon.

**ETO_LIQ_MAG_007, *Methanofollis* sp. nov. (HQ).** A single circular contig of 2.87 Mb. ANI 86.0%, AAI 84.7% (genus range 71.4–94.7%), POCPu 70.8% to *M. tationis*. Its own branch in *Methanofollis*. **New species supported.** Two *mcrBGA* operons, one next to a complete *mtrEDCBAFH* operon.

**ETO_LIQ_MAG_009, *Chakrabartyella* sp. nov. (HQ).** ANI 83.1%, AAI 86.6% (genus range 63.8–100%), POCPu 53.0%. Within *Chakrabartyella* (88.6/85). **New species supported.**

**ETO_LIQ_MAG_010, *Methanosarcina* sp. nov. (HQ).** ANI 89.3% (*M.* sp052057765), AAI 91.3% (genus range 75.3–99.5%), POCPu 69.8% to *M. barkeri*. **New species supported.** One complete *mcr* operon. The medium-quality MAG ETO_LIQ_MAG_018 (*Methanosarcina*) is a separate lineage (93.9% ANI to a GTDB genome).

**ETO_LIQ_MAG_005 and ETO_LIQ_MAG_011, two further *Methanofollis* spp. (MQ).** ANI 88.7% and 89.3% (*M. formosanus*), AAI 88.5% and 89.6%, POCPu 69.7% and 66.9%. They are sister to each other (100/100) but 86.8% ANI apart, so two species; both are < 80% ANI to ETO_LIQ_MAG_007. **New species supported.** Both meet the completeness and contamination limits but **lack a 16S rRNA gene** (011 has no rRNA at all), so they are not ideal type genomes. 005 carries two *mcr* operons.

**ETO_LIQ_MAG_008, DIGJ01 sp. nov. (Lachnospiraceae; MQ).** ANI ≤ 80.9%, AAI 81.3% (genus range 68.4–89.2%), POCPu 49.3% (just below 50%; DIGJ01 is a broad placeholder genus). **New species supported.** Meets the SeqCode type-genome checks (16S present).

**Controls.** ETO_LIQ_MAG_002: ANI 97.65% to *Gaoshiqia sediminis*, sister to the type strain A06ᵀ. ETO_LIQ_MAG_004: ANI 98.92% to *Lysinibacillus sphaericus_D*, inside the *L. sphaericus* complex (96.3/100). **QC passed for A1, A3 and A6.**

**Overall:** 9 new species-level lineages among the candidates, and no new genus established by threshold criteria. Together with A12, about 18 of the 29 MAGs are below 95% ANI to every GTDB r232 genome.


## dDDH (TYGS)
GGDC formula 2 (d4), recommended for incomplete genomes. Species boundary 70%. Source: 12 TYGS result sets (`TYGS_links_SP09.tsv`), downloaded with the TYGS API (`digital_ddh_tsv`, `identification_table`). Comparisons include TYGS-selected type strains ("U vs. T") and the uploaded panel references ("U vs. U").

| MAG | genomes compared | max dDDH d4 (CI) | against |
|---|---|---|---|
| ETO_SED_MAG_001 | 42 | 28.4 (26.1–30.9) | *Tindallia magadiensis* Z-7934ᵀ |
| ETO_LIQ_MAG_001 | 20 | 29.8 (27.4–32.3) | *Neobacillus pocheonensis* KCTC 13943ᵀ |
| ETO_LIQ_MAG_003 | 6 | 34.1 (31.6–36.6) | MTP4 sp9310u |
| ETO_LIQ_MAG_007 | 44 | 29.7 (27.3–32.2) | *Methanofollis* sp028702905 |
| ETO_LIQ_MAG_009 | 158 | 35.2 (32.8–37.7) | *Inediibacterium massiliense* Mt12ᵀ |
| ETO_LIQ_MAG_010 | 44 | 35.7 (33.3–38.3) | *Methanosarcina acetivorans* |
| ETO_LIQ_MAG_005 | 44 | 36.1 (33.6–38.6) | *Methanofollis formosanus* DSM 15483ᵀ |
| ETO_LIQ_MAG_011 | 44 | 37.8 (35.3–40.3) | *Methanofollis formosanus* DSM 15483ᵀ |
| ETO_LIQ_MAG_008 | 68 | 36.8 (34.3–39.3) | *Lacrimispora brassicae* YZC6ᵀ |

- **All 9 candidates are below 70% against every genome compared, so dDDH supports every species call.**
- Between the Etoliko *Methanofollis* MAGs: 005 vs 011 = 32.2%, 007 vs 005 = 20.8%, 007 vs 011 = 20.9%, confirming three species.
- At these low values (< ~35%), formula 2 cannot rank distant relatives, so the "max" genome can come from another genus (e.g. *Tindallia*, *Inediibacterium*). Genus placement rests on AAI and the trees, not on dDDH.
- Coverage gaps: for ETO_SED_MAG_001 and ETO_LIQ_MAG_001 the TYGS jobs compared only type strains, so most uncultured panel references show "n.d." in the panels; ETO_SED_MAG_001 was not compared with its nearest ANI reference (*Sedimentibacter* sp032006665, ANI 86.3%). The controls were not submitted to TYGS.

## A12: ANI screen of the 16 unnamed medium-quality MAGs (skani vs all GTDB r232 genomes)
- No GTDB genome ≥ ~80% ANI: SED_002, SED_003, SED_004, LIQ_008, LIQ_012, LIQ_014, LIQ_019, LIQ_022, LIQ_023.
- 80–90%: LIQ_005, LIQ_011, LIQ_017 (86.8), LIQ_020 (83.9), LIQ_024 (88.7).
- Borderline (close to 95%; may be the named species, low completeness): LIQ_016 → *Methanobacterium paludis* 94.3; LIQ_018 → *Methanosarcina* sp. 2.H.T.1A.6 93.9.
- Shared between Etoliko genomes (≥95% ANI): **SED_MAG_002 = LIQ_MAG_012** (99.75%, *Chakrabartyella*, recovered from both sediment and enrichment); **YGCA07 = YGCA10** (99.94%, same isolate strain).

## Method notes (since the first pass)
- A2: GTDB r232 metadata and trees downloaded (they were not on disk). Cytobacillaceae set capped at 150 (all 96 type strains kept).
- A4: only the needed pairs (candidate vs references, all pairs within each genus) were computed; the full all-vs-all would have taken about 25 h.
- A5: for genera without a type-species genome, a type strain or GTDB representative was used (`A5_pocp/pairs.tsv`, column `basis`). POCPu (unique matches) is used for decisions; raw POCP exceeds 100% where paralogs are present.
- A6: reference IDs prefixed `ref_` (they clash with GTDB-Tk internal IDs). Failed first-run outputs remain in `A6_trees/<set>/` (deletion blocked by a safety check); valid results are in `A6_trees/v2/`. ModelFinder used instead of fixed LG+F+R10 (several sets have fewer than 15 taxa).
- A7: the culturomics isolate (YGCC01) was annotated with Prokka, so the latest Bakta was used. The shared Bakta light DB has an outdated AMRFinderPlus DB, so a local copy was refreshed. The light DB gives no COG. eggNOG-mapper could not be used: the local eggNOG DB is version 7.0.0, but every eggNOG-mapper release (2.1.12, 2.1.15) reads only 5.0.2. COG categories therefore come from COGclassifier 2.0.0 (RPS-BLAST vs NCBI COG; the tool Proksee uses): 74–83% of CDS assigned. Archaeal mcr/mtr labels use COG names (Bakta light leaves them hypothetical); the restriction enzyme 'McrA' (COG3183) was excluded.
- New conda envs (no existing env modified): sp09_aai, sp09_bakta, sp09_cog, sp09_cgview; sp09_emapper and sp09_emapper215 were created for the failed eggNOG attempt and can be removed.

## Input notes
- The stated input folder `99_thesis_audit/SP09_final/inputs/` **does not exist**, and neither does the path cited in the prompt (`04_Shotgun_Metagenome/Tables/gtdbtk.bac120.summary.tsv`).
- Crosswalk: `07_mags/MAG_CROSSWALK.tsv` (md5 0080bd3351654f5c3f2533e1ebbc1b6c), the only copy in workspace_2.
- Reference GTDB-Tk summary: `08_taxonomy/gtdbtk/classify/gtdbtk.{bac120,ar53}.summary.tsv` (GTDB-Tk 2.7.2, r232). A second copy, `SP07_final/B4_unbinned_phyla/existing/gtdbtk/`, gives identical classifications and RED values for all 8 genomes, so the choice does not affect any result.

## A9: NCBI submission copy. **DONE, QC PASSED** (rerun 2026-10-21 after the student chose the inputs)
First pass stopped on ambiguous assemblies and LIQ reads. Student's choice: **polished assemblies** (`04_polish/{SED,LIQ}_asm_final.fasta`) and the **smaller LIQ read file** (`01_qc/LIQ_reads_min1000.fastq.gz`, 1,265,566 reads). SED reads: `00_input/sediment_reads_clean.fastq.gz` (6,827,767 reads, the only exact match).

- Files: 29 in `MAGs_raw/`, 29 in `MAGs_for_GenBank/`. All 29 match the crosswalk n_contigs and bp. ✅
- Removed contigs: **0** (none < 200 bp). Terminal Ns trimmed: 0. Internal runs of ≥10 N: 0. Total bp: raw 92,504,985 = cleaned 92,504,985 before adapter trimming; 92,504,576 after (−409 bp adapters). ✅
- FCS-adaptor 0.5.5 (Docker): **11 adapter hits in 7 MAGs** (ETO_LIQ_MAG_013 ×4, ETO_SED_MAG_002 ×2, ETO_LIQ_MAG_020/022/024, ETO_SED_MAG_003/004), all ACTION_TRIM at contig ends. Listed in `metadata/fcs_report.tsv`, then **trimmed on the student's instruction**: 409 bp removed; 9 end hits cut off, 2 internal hits excised with the contig split in two (no artificial join; contigs renumbered, header map updated; `metadata/adapter_trim_log.tsv`). Re-screen: 0 hits in all 7. No hits in the 8 HQ genomes evaluated in SP09. FCS-GX is not installed; NCBI will run it.
- Metadata: 29/29 rows in `biosample_MAGs.tsv` and `genome_submission.tsv`, plus 2 rows in `biosample_metagenomes.tsv`. `organism` for ETO_LIQ_MAG_014 set to `Christensenellales bacterium` by the student (GTDB f__/g__UBA6094 are placeholders). New umbrella BioProject chosen. **Empty cells:** none. `derived_from` holds the placeholder SAMNxxxxxxxx until the metagenome BioSamples exist.
- Originals unchanged (md5 of all 37 sources identical before and after). ✅
- Added beyond the requested layout: `*_contig_rename.tsv` (needed to link the polished contig names to the Flye `assembly_info.txt`), `metadata/cleaning_stats.tsv`, `removed_short_contigs.tsv`, `internal_N_runs.tsv`, `fcs_raw/`.
- Single-contig circularity (Flye circ. flag): ETO_LIQ_MAG_001 = Y, ETO_LIQ_MAG_007 = Y.

## A1: GTDB-Tk classify_wf. **QC PASSED**
GTDB-Tk 2.7.2, GTDB r232, `--keep_intermediates`. Outputs: `A1_gtdbtk/out/`, `A1_gtdbtk/A1_classification.tsv`, `A1_gtdbtk/A1_RED_medians.tsv`. Two earlier launches failed and are kept as logs: the default data path pointed to r220, and the tmpdir path was too long for a Unix socket. Both were fixed (r232 set explicitly; tmpdir reached through the short link `/tmp/sp09tmp` → `A1_gtdbtk/tmp`).

- **QC 1:** bacterial RED reproduces the existing summary exactly (|Δ| = 0.00000 for all 3 bacterial candidates; archaeal values also identical). **Pass.**
- **QC 2:** controls are assigned to their named species by ANI screen: ETO_LIQ_MAG_002 → *Gaoshiqia sediminis* (GCF_025907915.1, ANI 97.65, AF 0.812); ETO_LIQ_MAG_004 → *Lysinibacillus sphaericus_D* (GCF_038593855.1, ANI 98.92, AF 0.906). **Pass.**

Rank medians in r232 and the ±0.1 intervals used by GTDB-Tk:

| domain | rank | median RED | interval |
|---|---|---|---|
| Bacteria | genus | 0.9438 | 0.8438–1.0438 |
| Bacteria | family | 0.7973 | 0.6973–0.8973 |
| Archaea | genus | 0.9251 | 0.8251–1.0251 |
| Archaea | family | 0.7588 | 0.6588–0.8588 |

## First-pass notes (GTDB-Tk + RED only, superseded by the results above)
- **ETO_LIQ_MAG_001** (f__Cytobacillaceae; g__): RED 0.912. This lies inside the bacterial genus interval (0.844–1.044) and outside the family interval (upper bound 0.897), so RED places it at **genus rank**, consistent with a new genus within Cytobacillaceae. Single circular contig (Flye circ.=Y), 4.90 Mb, meets SeqCode type-genome criteria. **New genus is indicated by RED but not tested**: AAI (A4), POCP (A5) and a de novo tree (A6) were not run.
- **ETO_SED_MAG_001** (g__Sedimentibacter; s__): RED 0.989, placed in the genus by topology; no species within the ANI radius. **Indicates new species, not confirmed** (no skani/FastANI vs the full genus, A3 not run).
- **ETO_LIQ_MAG_003** (g__MTP4; s__): RED 0.917. Classification method: "taxonomic novelty determined using RED", yet a genus was still assigned. **Indicates new species, not confirmed**; the genus assignment via RED should be checked with AAI and a tree.
- **ETO_LIQ_MAG_007** (g__Methanofollis; s__): RED 0.936, topology; single circular contig (Flye circ.=Y). **Indicates new species, not confirmed.**
- **ETO_LIQ_MAG_009** (g__Chakrabartyella; s__): RED 0.978, topology. **Indicates new species, not confirmed.**
- **ETO_LIQ_MAG_010** (g__Methanosarcina; s__): RED 0.964, topology. **Indicates new species, not confirmed.**

All 8 meet the SeqCode type-genome checks (completeness >90%, contamination <5%, 16S present, ≤100 contigs). See `A8_table/novel_taxa_summary.{tsv,md}`. Sediment breadth is from `26_crossmap/cross_abundance.tsv` (the Table 3.10 source; it agrees with `47_master/MASTER.tsv`).

## A8: summary table. **Done**
11 rows (9 candidates + 2 controls), 33 columns (incl. dDDH max, CI, nearest): genome statistics, CheckM2/GUNC, rRNA/tRNA, CDS/COG, RED, ANI/AF, AAI, POCPu, tree placement, sediment breadth, SeqCode check, supported rank. Built by `A8_table/build_a8.py`.

## Medium-quality MAGs (not in the core SP09 set; see A12/A13 for the ANI screen and the 3 SeqCode-eligible ones)
- ETO_SED_MAG_002 (SED_das_vamb_004), medium quality
- ETO_SED_MAG_003 (SED_das_metabat_018), medium quality
- ETO_SED_MAG_004 (SED_das_semibin_003), medium quality
- ETO_LIQ_MAG_005 (LIQ_das_vamb_011), medium quality
- ETO_LIQ_MAG_006 (LIQ_das_semibin_002), medium quality
- ETO_LIQ_MAG_008 (LIQ_das_semibin_003), medium quality
- ETO_LIQ_MAG_011 (LIQ_das_vamb_017), medium quality
- ETO_LIQ_MAG_012 (LIQ_das_semibin_004), medium quality
- ETO_LIQ_MAG_013 (LIQ_das_metabat_075), medium quality
- ETO_LIQ_MAG_014 (LIQ_das_vamb_015), medium quality
- ETO_LIQ_MAG_015 (LIQ_das_semibin_006), medium quality
- ETO_LIQ_MAG_016 (LIQ_das_metabat_047), medium quality
- ETO_LIQ_MAG_017 (LIQ_das_vamb_004), medium quality
- ETO_LIQ_MAG_018 (LIQ_das_metabat_030), medium quality
- ETO_LIQ_MAG_019 (LIQ_das_metabat_009), medium quality
- ETO_LIQ_MAG_020 (LIQ_das_semibin_016), medium quality
- ETO_LIQ_MAG_021 (LIQ_das_metabat_078), medium quality
- ETO_LIQ_MAG_022 (LIQ_das_metabat_022), medium quality
- ETO_LIQ_MAG_023 (LIQ_das_semibin_010), medium quality
- ETO_LIQ_MAG_024 (LIQ_das_metabat_059), medium quality
- ETO_LIQ_MAG_025 (LIQ_dassub_semibin_017), medium quality

Outputs were written only inside `SP09_final/`, plus the temp link `/tmp/sp09tmp`. No workspace_2 source file was modified (all were read-only inputs).
