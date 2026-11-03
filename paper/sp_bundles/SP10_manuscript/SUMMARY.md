# SP10 manuscript export bundle: summary

Built 2026-10-26 from `99_thesis_audit/prompts/SP10_manuscript_export_bundle.md`. **No analysis was re-run.** Existing outputs were re-read, figures were re-drawn in English and tables exported. Nothing outside `SP10_manuscript/` was changed, and no conda environment was created or modified (the figures used `workspace_2/envs/smetana`, as in the thesis, and `r_viz`). Every input path and md5 is in `inputs_md5.tsv` (210 files); tool versions are in `versions.tsv`.

## Things to fix in the manuscript (found while checking)
1. **dbCAN3 "931 proteins" has no source.** The CAZyme total over the 29 MAGs is **1,774**, and no subset found equals 931. The per-MAG numbers (421 / 165 / 180 / 118) are correct. (E4 #4)
2. **Self-mapping 34.5% / 63.4% count supplementary alignments in the denominator.** Per primary read the rates are **35.8% / 73.0%**. They are against the polished assemblies, using the deduplicated reads before the 1 kb filter. (E4 #3)
3. **The LIQ library is 5 runs, not 9 pooled conditions.** The run files are named by depth interval and amendment, and one carries "Sul" (sulfate?). This contradicts the "nine conditions pooled, no sulfate" text in `MEDIUM_DEFINITION.md` / `CHAPTER_STRUCTURE.md`. The amendment codes are not defined on the server: **check the lab book.** (E4 #5)
4. **The sediment "10 runs" include one file twice** (`ligation_new.fastq` = `Ligation_promithion.fastq`, byte-identical). Deduplication removed exactly that copy, and RBK_test was excluded, so **8 distinct sediment runs** went into the assembly. (E4 #6)
5. **The DNA extraction kit is recorded nowhere on the server**; the library kit is known for 7 of 10 sediment files (from file names or read headers; unknown for `metagenome_2`, `metagenome_3` and `metagenome_new`) and for 0 of 5 culture files. (E4 #6)
6. The SingleM "≈ 0.03% of filled phylum coverage" cannot be checked here: `singlem_filled_phylum.py` is not on the server. The 0.38 coverage is confirmed. (E4 #2)

## Status per item
### E1: figures (`E1_figures/`, PNG 600 dpi + PDF; caption drafts in `captions.md`)
| figure | status | how |
|---|---|---|
| S_MAG_quality | done | thesis script, English labels |
| S_function_heatmap | done | thesis script; it reads `47_master/MASTER.tsv`, not the two panel files named in the prompt (unchanged so numbers do not move) |
| S_exchange_matrix | done | thesis script; matrix = all 610 `pairs_detailed.tsv` rows summed (re-checked: 32 pairs, max 10.78); completeness printed, < 90% marked ° |
| S_compound_classes | done | thesis script |
| S_biogeography | done | thesis script; public runs at their real values (0.000–0.007) |
| S_dsr_binning | done | thesis script (inputs `57_dsr_bin/`) |
| S_ani_isolates | done | thesis script |
| S_mcrA_nomination | done | thesis script (reads `45_kofam_full/all_hits.tbl` by bin; agrees with `90_mag_bin_mapping/mcr_scores_by_mag.tsv`; `54_mcr_fragmentation` not used) |
| S_mcrA_tree | done | **new plot** from the SP05 tree; 2 restriction-enzyme McrA pruned; MCR α / ACR α / ACR other / component A2 coloured; both joined MAG_017 proteins included; ACR α not monophyletic, so no root drawn |
| S_dsr_trees | done | **new plot** from the SP05 DsrA/DsrB consensus trees, one panel each, UFBoot ≥ 70; drawn as written (phangorn/phytools not installed, so midpoint rooting was not used rather than modify an env) |
| S_phylogenomics_archaea / _bacteria | done | thesis script; bacteria tree drawn without a root |
| S_rrna_corrected | done | new plot (no thesis script) |
| S_arc_otu22 | done | copied unchanged from the delivered Figure 7 (already English, caption-free); **16.0 cm** wide |
| S_viruses_defence | done | new plot |
| S_bgc | done | new plot; 141 regions, 33 with a MIBiG hit |

Style changes, made in a local copy of `etoliko_figstyle.py`: Arial, white background, "." decimal mark, Bacteria #1F5FA6 / Archaea #D9622B, ≥ 7 pt, 17 cm. `caption()` no longer draws, and titles were blanked. Two rendering-only fixes in the copy of `fig14_15_mcra.py`, needed for matplotlib 3.11: an Agg canvas for the probe figure, and an early `plt.close(probe)` removed. No PDF contains Greek text (checked with pdftotext). Most figures are 13–16.6 cm wide at their natural size; none is wider than 17 cm. **Font check:** the thesis scripts were authored with 7.5–10.5 pt sizes; the style module sets ≥ 7 pt, but a few hard-coded sizes in the thesis scripts (e.g. small tree tip labels) may be below 7 pt. Look at the PNGs before submission.

### E2: tables (`E2_tables/`, 29 TSV; gaps listed in `E2_tables/E2_notes.tsv`)
All 21 requested tables were written; some are split into several files:
- T10: `T_pyrrolysine` + `_summary`
- T11: `T_readlevel_markers` + `_mtbB_besthit`
- T14: `T_smetana_pairs` + `_viable_pairs`
- T16: `T_mobilome_*` (×4)
- T21: `T_unbinned_fraction` + `_marker_levels`

Text sources are copied to `E2_tables/text/`. **Partly filled (NA, with reason in `E2_notes.tsv`):**
- T2: mean Q / % > Q15 before filtering;
- T6: country;
- T8: KO-list column, and the panels not present in `function_panels.tsv`;
- T9: PUL and T9SS counts;
- T15: thresholds;
- T19: SingleM filled-phylum table and unnamed fractions;
- T21: % reads mapped to MAGs.

T3 uses Flye `assembly_info.txt` / `.fai` because QUAST was never run. Per-run read statistics (T2) were computed for SP10 with `seqkit stats` on the original per-run FASTQs.

### E3: novel-taxa trees (`E3_trees/`)
SP09 originals were copied to `original_SP09/` and checked. All labels are English, but **these had caption text inside the image and were re-drawn** (`redrawn/`, from the same Newick files and tables, E1 style):
- the 10 per-candidate trees (title + subtitle), including the A13 trees (methanofollis, digj01);
- the 9 species panels (caption);
- the 2 combined trees (caption).

The 9 circular maps (`A7_maps`) have no caption and were not re-drawn, but their legend text is 5.5–6 pt. The combined trees stay at their original size (38 × 36 and 46 × 33 cm): 303 tips cannot be legible at 17 cm. The A6 trees are at `SP09_final/A6_trees/`, not `A6_trees/v2/` as the prompt says.

### E4: consistency checks
Done: `E4_checks/E4_report.md`, plus `per_run_fastq_md5.txt` for the 15 per-run FASTQs.

## Folder
```
SP10_manuscript/
  E1_figures/   16 × (PNG+PDF), captions.md, scripts/ (copied+translated thesis scripts, style module, new plots), render logs
  E2_tables/    29 TSV, E2_notes.tsv, text/, scripts/build_e2.py (+ per_run_seqkit_stats.tsv)
  E3_trees/     original_SP09/, redrawn/{trees,species_panels,combined_tree}, scripts/
  E4_checks/    E4_report.md, per_run_fastq_md5.txt
  inputs_md5.tsv  versions.tsv  SUMMARY.md
```
Note: `E1_figures/scripts/translate.py` records the label translation applied to the copied scripts. It has already been applied, so running it again stops at the first string that is no longer present.
