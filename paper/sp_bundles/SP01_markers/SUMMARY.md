# SP01 summary: the "31 MAGs" and the curated marker matrix

All four QC controls pass (`markers/QC.txt`): mcrA is called in 7/8 archaea (not in ETO_LIQ_MAG_017, which has the known split mcr) and 0/21 bacteria. The YGCC01 ure cluster JFCCLLDG_00871–00877 (7/7 genes) and YGCA01_C arcABC BPCFBLIA_02024–02026 are detected. No dsrAB is called in any methanogen.
Levels: gene presence is L2 and module completeness is L3. Absence is stated only for genomes ≥ 90% complete; every other absence is marked NA.

**Deviations from the prompt (both approved by the user on 2026-09-29):**
- KofamScan was not used. Step 45 used hmmsearch plus ko_list thresholds, so the isolates were scored the same way.
- The "≥60% of profile length" rule failed the YGCC01 urease control on full-length UreA and UreG, because KOfam profile lengths are inflated by fused homologs. It was replaced: a hit is a "fragment" only if its gene is truncated (Prodigal partial≠00, or the CDS touches a contig end).
- HydDB web did not respond, so hydrogenases were classified against the HydDB-labelled Greening-lab reference sets (see `versions.tsv`).
- The prompt's "mtrA K00584" is mtrH; K00577 (mtrA) was added.

## Which "31"?
There is no 31-genome set. `Combined_analysis/00_inputs/metagenomes/mags/` holds 29 FASTAs, and all 29 are byte-identical (md5) to the final workspace_2 MAGs. None is a GUNC-rejected bin, a workspace_1 bin or of unknown origin (`provenance_31.tsv`).

"31" occurs only in a step name and in prose. The "34 genomes" are those 29 MAGs plus 5 isolates, and the file is named `all_36_...`.

The names MAG_Methanofollis_aquaemaris/formosanus, MAG_Syntrophomonas_sp and MAG_Lysinibacillus_sphaericus are hand labels. They are typed in `02_run_read_recruitment.sh` (lines 60–63) for LIQ_das_semibin_016, LIQ_das_vamb_017, LIQ_das_semibin_003 and LIQ_das_vamb_002, i.e. ETO_LIQ_MAG_020, 011, 008 and 002.

The Combined_analysis marker matrix is not a curated call:
- For MAGs it counted raw `hmmsearch -E 1e-5` rows with no threshold.
- For isolates it used a gene-name substring search.
- Several KO labels are wrong: ureA and ureC are swapped, "hydA" is hoxH, and "echA" is MttB (`combined_analysis_marker_method.md`).

A flat-E-value call on the same targets (E<1e-10, as in step 16b) would add 2,785 proteins to 309 curated calls, an inflation of about 900% (`markers/flat_evalue_inflation.tsv`).

## Answers (29 MAGs; tables: `markers_long.tsv`, `markers_wide.tsv`)
- **Urease:** no MAG carries ureA, ureB or ureC above threshold. Among the 34 genomes the ure genes are called only in YGCC01. In the sediment reads, 89 reads have a best hit to ureC: 13.0 per million, 0.013 × rpoB (6,733 rpoB reads; glnA and tuf also detected). This is L2 at community level only and is not attributed to YGCC01 or any MAG (`urease_reads.tsv`).
- **dsrAB:** none in any MAG or isolate. Steps 55–59 place both reductive dsrAB operons on unbinned contigs. The LIQ_contig_19299 link to ETO_LIQ_MAG_022 remains an open hypothesis only.
- **aprAB:** called only in ETO_LIQ_MAG_022 (1 copy each). sat is non-diagnostic and is not interpreted.
- **True Nif:** Nif group III homologs pass the tree check in ETO_LIQ_MAG_002, 003, 005, 007, 010, 018 and 020 (nifH+nifD+nifK), in LIQ_MAG_015 (nifH+nifK) and in LIQ_MAG_016 (nifH+nifD) (nearest-reference placement, `nif_tree_verdicts.tsv`).
  - The nifH homologs in LIQ_MAG_005, 007, 011, 016 and 017 that fall with CfbC are **not nitrogenases**; they are F430 biosynthesis genes.
  - Hits where fewer than 4/5 nearest references agree are reported as AMBIGUOUS and not counted: LIQ_MAG_006, 018 (second copy), 019, 020 (second copy), SED_MAG_003, YGCA07 and YGCA10.
- **Hydrogenases** (`hydrogenase_classes.tsv`):
  - [NiFe] catalytic subunits of groups 1a/1d/1k, 3a/3c (archaea) and 4b/4d–4i.
  - [FeFe] group A2/A3 in LIQ_MAG_019 (a 120-aa hit), LIQ_MAG_021, LIQ_MAG_023 and SED_MAG_003; group C in LIQ_MAG_003 and LIQ_MAG_021.
  - Short hits (under 200 aa) are fragments and are flagged there.
- **FDH:** fdhF/fdoG (K00123) is called in ETO_LIQ_MAG_003, 004, 005, 007, 011, 016, 017 and 018. fdhA (K05299 and K22015) is called in no MAG; K22015 appears only in the isolate YGCA01_Tp.
- **β-oxidation (M00087, L3):** complete (3/3 steps) only in YGCC01. No MAG is complete: at most 2/3 (LIQ_MAG_001, 004, 025). For MAGs under 90% complete, the missing steps are NA (`markers/M00087_completeness.tsv`).

Provenance of each cell: input `.faa` (md5 in `markers/inputs.tsv`, `isolates/inputs.tsv`) → hmmsearch command in `COMMANDS.log` → row in `markers/calls_detail.tsv` → row in `markers_long.tsv` → statement above.
