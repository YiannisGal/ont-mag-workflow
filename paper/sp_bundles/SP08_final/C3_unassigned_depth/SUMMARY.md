# C3 — unassigned / uOTU share vs depth (H2)
Output: `C3_unassigned_vs_depth.tsv` (tests), `C3_per_depth_profiles.tsv` (14 values per series), `QC_inputs.txt`. Script: `scripts/c3.R` (+ `scripts/load.R`).
Aggregation: counts summed over the extractions of each depth → % of non-plastid reads (SILVA o__Chloroplast removed). Spearman p = 9,999 permutations (exact p also given); Mann–Kendall = R `Kendall`. BH across all 9 C3 tests.

## Inputs chosen (no authoritative path given in the prompt; choices made from file provenance)
| dataset | table | pass 1 (uOTU threshold) | pass 2 |
|---|---|---|---|
| V3–V4 | Bacteria/out_3/otutab_trim.txt (1,025 OTUs = the classified set) | SILVA 138.2, 95% (qza export) | SILVA 138.2, 75% |
| V4–V5 | Archaea/out/otutab_trim.txt (395) | **pre-138.2 SILVA (RESCRIPt), 95%** — `taxonomy_blast_silva_138_2.qza` was actually run at 75% (provenance), so no 138.2 95% run exists; the exported `taxonomy_138_2/taxonomy.tsv` is the 75% run | SILVA 138.2, 75% |
| ITS2 | Eukaryotes/out/otutab.txt restricted to the 334 OTUs of otus_trimmed.fa (400–519 bp) | **UNITE_new, 95%** (no 97% run exists) | UNITE, 75% |
| ONT14 | Etoliko_new/Results_new/final_output/cluster_tab_norm.txt (882 clusters, 68 libraries) | SILVA 138.2, 95% | **SILVA 138, 90%** (only second pass on file) |
Depth: Bacteria Metadata.txt by YG number (V4–V5, ITS2 via the same YG number, as in SP03); ONT via Metadata_minion.tsv. ITS2 was not rarefied to 27,312 (no sample has ≥ 27,312 reads; shares are relative, so rarefying does not change them materially).

## QC
- ONT14 uOTU share overall = **18.9%** (thesis Table A.5 ≈ 18–19%) → **PASS**.
- V3–V4: **283 non-plastid OTUs (284 incl. plastids) unassigned at 95%; thesis says 250 → MISMATCH (FAIL, not adjusted).** The same count (284) is in the older `taxonomy/` export, so the thesis value probably comes from another OTU set or version.
- V4–V5, ITS2: no thesis reference value in the prompt → not checked. V4–V5 and ITS2 first passes are substitutes (see above), so their uOTU rows are not strictly the thesis definition.

## Result
No series has a depth trend after BH (all p_BH ≥ 0.25). Nominal only: archaeal uOTU (ρ = −0.53, MK p = 0.049) and DPANN share (ρ = −0.59, p = 0.028; MK p = 0.038) decrease with depth. "Unassigned" ≠ novel; trends can reflect reference coverage.
