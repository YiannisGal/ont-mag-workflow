# E4: consistency checks (report only; nothing was changed)

Paths are relative to `workspace_2/` unless they start with `workspace/` (= workspace_1, `Metagenome_analysis_PhD/workspace/`). Line numbers count the header as line 1.

| # | statement | verdict |
|---|---|---|
| 1 | sediment breadth 0.8936 / 0.6912 gap | **confirmed** |
| 2 | SingleM Methanomassiliicoccales ≈ 0.38 coverage, ≈ 0.03% of filled phylum coverage | 0.38 **confirmed**; the % share **cannot be checked** (script/output not on server) |
| 3 | 29,404 / 19,374 contigs; 34.5% / 63.4% self-mapping | contig counts **confirmed**; mapping rates **confirmed as printed, but the denominator includes supplementary alignments** |
| 4 | dbCAN3 931 proteins; MAG_002 421, MAG_013 165, MAG_015 180, MAG_020 118 | per-MAG sums **confirmed**; **931 not found** (file total = 1,774) |
| 5 | LIQ library design | the run files support **5 runs, each from one depth interval, with amendment codes incl. "Sul"**; no run sheet on the server |
| 6 | DNA extraction kit and ONT kit per run | ONT kit only from file names / read headers; **extraction kit not recorded anywhere on the server** |

## 1. Sediment breadth
Source `26_crossmap/cross_abundance.tsv` (column `sed_breadth`, raw file re-read).
- Line 12: `ETO_LIQ_MAG_011  … sed_breadth 0.8936`
- Line 15: `ETO_LIQ_MAG_014  … sed_breadth 0.6912`

All 29, descending (value, line, MAG): 1.0000 (30, SED_004), 0.9999 (27, SED_001), 0.9998 (29, SED_003), 0.9855 (28, SED_002), 0.9780 (21, LIQ_020), 0.9714 (13, LIQ_012), 0.9138 (26, LIQ_025), 0.9056 (9, LIQ_008), 0.9031 (16, LIQ_015), **0.8936 (12, LIQ_011)** ‖ **0.6912 (15, LIQ_014)**, 0.6573 (22, LIQ_021), 0.6170 (4, LIQ_003), 0.4649 (25, LIQ_024), 0.4330 (10, LIQ_009), 0.4287 (7, LIQ_006), 0.3755 (6, LIQ_005), 0.3183 (8, LIQ_007), 0.2607 (23, LIQ_022), 0.2578 (20, LIQ_019), 0.2503 (24, LIQ_023), 0.2239 (5, LIQ_004), 0.1796 (19, LIQ_018), 0.1582 (2, LIQ_001), 0.0825 (11, LIQ_010), 0.0276 (18, LIQ_017), 0.0248 (14, LIQ_013), 0.0223 (3, LIQ_002), 0.0090 (17, LIQ_016).
The gap between the 10th and 11th genome (0.8936 → 0.6912) is the largest step in the list; the 10 genomes above it are the "in-situ ten".

## 2. SingleM Methanomassiliicoccales
- `22_singlem/SED_profile.tsv` line 302: `sediment_reads  0.38  Root; d__Archaea; p__Thermoplasmatota; c__Thermoplasmata; o__Methanomassiliicoccales; f__JAAEEP01; g__JAAEEP01` → coverage **0.38**, confirmed.
- **Share of filled phylum-assigned coverage (≈ 0.03%): not checked.** `singlem_filled_phylum.py` and any output of it are not in workspace_2 or workspace (searched by name), and no notebook defines "filled" coverage. Recomputing it would mean guessing the definition, so the item is skipped. The script is probably on the laptop; run it there on `SED_profile.tsv` (md5 in `inputs_md5.tsv`).

## 3. Contigs and self-mapping
- Contigs: `03_assembly/sediment/assembly.fasta` **29,404** and `03_assembly/liquid_cultures/assembly.fasta` **19,374** (Flye). The polished assemblies used for binning, `04_polish/SED_asm_final.fasta` and `LIQ_asm_final.fasta`, have the same counts (29,404 / 19,374).
- Mapping: `02_mapping/mapping_summary.tsv` line 2 `SED_map_self 7810950 2692947 34.5`, line 4 `LIQ_map_self 3241803 2054488 63.4`.
  - **Which assembly:** the **polished** assemblies `04_polish/{SED,LIQ}_asm_final.fasta` (`scripts/05_map.sh` lines 38–39, 59–61; minimap2 `-ax map-ont --secondary=no`).
  - **Which reads:** `00_input/sediment_reads_clean.fastq.gz` and `00_input/liquid_reads.fastq.gz` (`00_config.sh` lines 150–151), i.e. the deduplicated reads **before** the ≥ 1 kb filter, not the `01_qc/*_min1000` reads used for Flye.
  - **Caveat on the denominator:** `total_reads` is the flagstat "in total" count, which includes supplementary alignments (`02_mapping/SED_map_self.flagstat.txt`: 7,810,950 total = 7,528,854 primary + 282,096 supplementary). 34.5% = primary mapped / all records. Per primary read, the rates are **35.77% (SED)** and **72.97% (LIQ)** (flagstat line 8, "primary mapped"). Counting supplementary alignments as mapped, the rates are 38.09% / 76.53%. The manuscript should state which one it uses; 35.8% / 73.0% is the conventional "% reads mapped".

## 4. dbCAN3
Source `30_dbcan/family_matrix.tsv` (columns `total`, GH, PL, CE, GT, AA, CBM) and `30_dbcan/cazyme_families.tsv` (one row per protein).
- MAG_002 = `ETO_LIQ_MAG_002` line 3: total **421** (= GH+PL+CE+GT+AA+CBM). MAG_013 line 14: **165**. MAG_015 line 16: **180**. MAG_020 line 21: **118**. All confirmed.
- **931 total proteins: not found.** Sum of `total` over the 29 MAGs = **1,774**, equal to the 1,774 data rows (unique locus tags) of `cazyme_families.tsv`. Other sums for comparison: GH 951, GT 539, CE 182, PL 40, AA 30, CBM 30; multidomain 317; Bacteroidota MAGs 653. No subset found in the files equals 931, and "931" does not appear in the notebooks. The manuscript number needs a source or should be replaced with 1,774.

## 5. Enrichment culture design behind the LIQ library
**No run sheet or sample sheet exists on the server** (searched `Metagenome_analysis_PhD/` to depth 6 for `sample_sheet*`, `*samplesheet*`, `final_summary*`, MinKNOW `report_*`; the only hit, `workspace/mags_samplesheet.csv`, is a MAG list). The closest primary evidence is the set of per-run FASTQs that make up the LIQ read set:

`workspace/provenance_harvest/02_qc/per_run_read_counts.tsv` (lines 1–5), files in `workspace/liquid_cultures/profile/fastqs/Liquid_Cultures_Metagenome/`:
| run file | reads | depth interval (cm) and codes in the name |
|---|---|---|
| `metagenome_3540glu.fastq` | 43,658 | 35–40, glu(cose?) |
| `metagenome_35-50_Fum_Sul_But.fastq` | 219,566 | 35–50, Fum / **Sul** / But |
| `metagenome_35-50_No_Pro_Dl.fastq` | 1,624,718 | 35–50, No / Pro / **Dl** |
| `metagenome_50-65-DL-No-But.fastq` | 22,900 | 50–65, **DL** / No / But |
| `metagenome_5-20_No-Pro-But_211023YG.fastq` | 904,634 | 5–20, No / Pro / But |
| total | 2,815,476 | = `liquid_cultures_all.fastq.gz` (line 19) |

- This matches the second record, `labnotes/workspace_2_lab_notebook.md` line 300: *"enrichment-culture libraries comprised five runs, each pooling two to three amendment conditions from one depth interval"*.
- It contradicts `labnotes/MEDIUM_DEFINITION.md` lines 98–100 (*"The nine conditions were pooled before sequencing … pooled for the liquid-culture shotgun library"*), `labnotes/CHAPTER_STRUCTURE.md` lines 24 and 54–57, and `labnotes/CORRECTIONS_INDEX.md` line 55 (*"The enrichment medium contains no sulfate"*).
- **"Sul" in run 2 suggests a sulfate amendment, and "Dl"/"DL" suggests DL-lactate**, but the codes are not defined anywhere on the server. "No" could mean nitrate or "no addition". Only the lab's culture records can settle what the codes mean and whether the 5 runs cover 9 conditions. **Not settled by a run sheet: needs the lab book.**

## 6. Kits per run
There are no MinKNOW reports, sample sheets or `final_summary` files on the server. Of the 15 per-run FASTQs, only two carry Dorado header metadata, and only `protocol_group_id`, not `kit=`. **DNA extraction kit: not recorded in any file on the server, for any run.**

| run file | library kit (evidence) | PCR-based? |
|---|---|---|
| `Sediment_Metagenome/ligation_new.fastq` | Ligation (header `protocol_group_id=GALIATSOS_LIGATION`; LSK114 assumed, kit code not written) | no |
| `Sediment_Metagenome/Ligation_promithion.fastq` | **byte-identical copy of `ligation_new.fastq`** (same size 4,959,637,787 B, same first read, same 3,897,829 reads) | no |
| `Sediment_Metagenome/Rapid_PCR_promithion.fastq` | Rapid PCR Barcoding, RPB (header `protocol_group_id=GALIATSOS_RAPID_PCR_BARCODING`, `barcode04`) | **yes** |
| `Sediment_Metagenome/TSIAMIS_EXP-PBC096.fastq` | PCR Barcoding Expansion EXP-PBC096 (name) | **yes** |
| `Sediment_Metagenome/TSIAMIS_SQK-NBD114-24.fastq` | Native Barcoding SQK-NBD114.24 (name) | no |
| `Sediment_Metagenome/metagenome_LSK_test.fastq` | Ligation, LSK (name) | no |
| `Sediment_Metagenome/metagenome_RBK_test.fastq` | Rapid Barcoding, RBK (name); 16S amplicon library, **excluded** (notebook vol2 lines 355–374) | amplicon |
| `Sediment_Metagenome/metagenome_2.fastq`, `metagenome_3.fastq`, `metagenome_new.fastq` | **not recorded** | unknown |
| 5 × `Liquid_Cultures_Metagenome/*.fastq` | **not recorded** (no kit in name or header) | unknown |

**Consequence for the read counts.** `sediment_all.fastq.gz` (11,426,683 reads) contains the ligation run twice. Deduplication removed exactly one copy: 11,426,683 − 3,897,829 = **7,528,854** = `sediment_dedup.fastq.gz` (`per_run_read_counts.tsv` lines 16–17). The ten sediment "runs" are therefore nine distinct runs plus one duplicate file. Of these, one (RBK_test) was excluded, which leaves **eight distinct sediment runs in the assembly read set**. The SED_LIGATION and Rapid PCR runs share flow cell PBG45329 (2026-03-19/20). MD5 checksums of all 15 files are in `per_run_fastq_md5.txt`.
