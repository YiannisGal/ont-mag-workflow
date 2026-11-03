# SP07 summary — final analyses B1–B6

Run order B6 → B1 → B2 → B5 → B3 → B4. Per-part details, commands, versions and md5s are in each part folder.

## B6 Arc_Otu22

| amplicon | md5 | length (primer-free) | identical positions | iddef0 % | iddef1 % | iddef2 % |
|---|---|---|---|---|---|---|
| Arc_Otu22_ms | 8bd6b191f8b4eaa35bd2d79b2126c333 | 415 | 405/415 (97.6%) | 97.6 | 25.9 | 97.1 |
| Arc_Otu22_server | 8bd6b191f8b4eaa35bd2d79b2126c333 | 415 | 405/415 (97.6%) | 97.6 | 25.9 | 97.1 |

| verdict (one line) |
|---|
| Server Otu22 = manuscript amplicon (identical, 415 bp primer-free). vsearch measures 405/415 (97.6%; 97.1% iddef2), not 411/415. The difference is vsearch aligning gene positions 841–850 as mismatches where BLAST uses single-base indels (SP06). |

QC: positive control gene vs itself: 100.0% ; server Otu22 vs ms amplicon (full, with primers): 100.0% (450 cols, 0 mism); primer-free identical: True

## B1 dsrA/dsrB reads (sediment library)

References: dsrAB@SILVA (Jun 2023, GTDB r214; 792 DsrA / 661 DsrB / 597 concatenated DsrAB after 90% clustering). Family- and order-level QC fail in every tree (family 75–81%, order 86–89%); the reductive/oxidative clade passes (98.7–99.4%). B1 is therefore reported at clade level only.

**Clade-level result (QC PASS; this is the thesis-usable B1 result):**

| clade | dsrA reads | dsrB reads | % of assigned (dsrA+dsrB) | concatenated dsrAB reads (27 reads spanning both genes) |
|---|---|---|---|---|
| Reductive | 52 | 45 | 85.1% | 21 |
| Oxidative | 9 | 8 | 14.9% | 4 |

**Family / genus calls: FAIL family- and order-level QC; exploratory only, not for thesis use:**

| clade (reductive/oxidative) | family / genus | dsrA reads | dsrB reads | % of assigned |
|---|---|---|---|---|
| Oxidative | Chlorobiaceae / Chlorobium_A | 1 | 0 | 0.9% |
| Oxidative | Chromatiaceae / unassigned | 1 | 0 | 0.9% |
| Oxidative | Sedimenticolaceae / unassigned | 4 | 1 | 4.4% |
| Oxidative | UBA6429 / Thiohalophilus | 0 | 1 | 0.9% |
| Oxidative | unassigned at family level / unassigned | 3 | 6 | 7.9% |
| Reductive | B3-LCP / unassigned | 1 | 0 | 0.9% |
| Reductive | C00003060 / unassigned | 1 | 0 | 0.9% |
| Reductive | CG2-30-66-41 / CG2-30-66-41 | 0 | 1 | 0.9% |
| Reductive | Desulfarculaceae / unassigned | 0 | 1 | 0.9% |
| Reductive | Desulfatiglandaceae / 4484-190-2 | 1 | 0 | 0.9% |
| Reductive | Desulfatiglandaceae / CG2-30-51-40 | 1 | 0 | 0.9% |
| Reductive | Desulfatiglandaceae / unassigned | 1 | 0 | 0.9% |
| Reductive | Desulfobaccaceae / unassigned | 4 | 0 | 3.5% |
| Reductive | Desulfobacteraceae / Desulfamplus | 2 | 0 | 1.8% |
| Reductive | Desulfobulbaceae / unassigned | 1 | 0 | 0.9% |
| Reductive | Desulfocapsaceae / unassigned | 2 | 0 | 1.8% |
| Reductive | Desulfosarcinaceae / Desulfatitalea | 1 | 0 | 0.9% |
| Reductive | Desulfosarcinaceae / Desulfosarcina | 1 | 0 | 0.9% |
| Reductive | Desulfotomaculaceae / Desulfotomaculum | 1 | 0 | 0.9% |
| Reductive | FEB-10 / Uncultured_family_lineage_9 | 2 | 0 | 1.8% |
| Reductive | PGXB01 / PGXB01 | 0 | 1 | 0.9% |
| Reductive | SG8-13 / SG8-13 | 0 | 1 | 0.9% |
| Reductive | SM23-35 / UBA2194 | 1 | 2 | 2.6% |
| Reductive | Syntrophomonadaceae / UBA4851 | 0 | 1 | 0.9% |
| Reductive | UBA1163 / UBA1163 | 0 | 1 | 0.9% |
| Reductive | UBA2268 / unassigned | 0 | 1 | 0.9% |
| Reductive | UBA3076 / UBA3076 | 1 | 1 | 1.8% |
| Reductive | unassigned at family level / unassigned | 31 | 35 | 57.9% |
| unassigned (no hit, segment < 60 aa, or clade LWR < 0.5) | none | 9 | 15 | none |

| QC line |
|---|
| dsrA: leave-one-out re-placement of 792 references -> function clade recovered 98.9%; family recovered 81.1% of 594 evaluable (taxon present elsewhere in tree). Threshold 95% -> FAIL |
| dsrB: leave-one-out re-placement of 661 references -> function clade recovered 99.4%; family recovered 74.9% of 470 evaluable (taxon present elsewhere in tree). Threshold 95% -> FAIL |
| dsrAB: leave-one-out re-placement of 597 references -> function clade recovered 98.7%; family recovered 81.2% of 464 evaluable (taxon present elsewhere in tree). Threshold 95% -> FAIL |

## B2 mcrA reads

**QC FAIL: see QC line. The table below is not for thesis use.**

| order | reads | nearest reference | % identity range | falls with MAG_011? |
|---|---|---|---|---|
| MCR Methanomicrobiales | 14 | UP|A0A8G1EFZ4|Methanofollis_formosanus | 84.9–99.1 | yes, 5 reads (identity to MAG_011: 95.4, 93.3, 90.2, 90.1, 90.2%) |
| MCR Methanosarcinales | 1 | UP|P07962|Methanosarcina_barkeri | 95.6–95.6 | no |

QC: mcrA: leave-one-out re-placement of 113 SP05 references -> function clade recovered 99.1%; order recovered 93.8% of 113 evaluable (taxon present elsewhere in tree). Threshold 95% -> FAIL placement model used: LG+I+G4; reads 15, blastx hits 15, placed (segment >= 60 aa) 15

References: v2 = SP05 references + step-58 + UniProt top-10 hits of the reads (KOfam-thresholded, 95% clustered); see common/reftrees/*/build_log.txt. Reference choice is read-guided, which favours lineages present in the sample.

## B3 community ~ geochemistry

| dataset | n | variables retained (VIF) | selected (ordiR2step) | adj R² | F | p | p (BH) | partial db-RDA | Condition(depth) adj R² | p |
|---|---|---|---|---|---|---|---|---|---|---|
| Illumina bacteria V3-V4 | 14 | LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate,TP_rank (dropped: Polysaccharides(17.1)) | Phenols,TKN | 0.245 | 1.529 | 0.0129 | 0.0258 | Phenols,TKN + Condition(depth), F = 1.696 | 0.091 | 0.0403 (BH 0.0604) |
| Illumina archaea V4-V5 | 14 | LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate,TP_rank (dropped: Polysaccharides(17.1)) | TKN,Phenols | 0.268 | 1.594 | 0.0027 | 0.0108 | TKN,Phenols + Condition(depth), F = 1.782 | 0.106 | 0.0131 (BH 0.0393) |
| ITS2 | 14 | LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate,TP_rank (dropped: Polysaccharides(17.1)) | none | 0.063 | 1.110 | 0.2795 | 0.3727 | LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate,TP_rank + Condition(depth), F = 1.032 | 0.019 | 0.4260 (BH 0.4260) |
| ONT bacteria | 9 | LOI950,Proteins,Phenols,NH4N,TKN,Sulfate,TP_rank (dropped: LOI550(Inf) > Polysaccharides(30.2)) | none | 0.036 | 1.043 | 0.5536 | 0.5536 | LOI950,Proteins,Phenols,NH4N,TKN,Sulfate,TP_rank + Condition(depth), F = none | none | none (BH none) |

| varpart: depth only / geochem only / shared / residual, per dataset |
|---|
| Illumina bacteria V3-V4: 0.000 / 0.091 / 0.128 / 0.781 (geochem set: Phenols,TKN) |
| Illumina archaea V4-V5: 0.007 / 0.106 / 0.076 / 0.812 (geochem set: TKN,Phenols) |
| ITS2: 0.032 / 0.019 / 0.044 / 0.904 (geochem set: LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate,TP_rank) |
| ONT bacteria: none / none / none / none (geochem set: LOI950,Proteins,Phenols,NH4N,TKN,Sulfate,TP_rank) |

```
Illumina bacteria V3-V4: depth-only Mantel r = 0.632 (p = 0.0001); ms v2 Mantel r ~ 0.60 -> PASS (same direction, similar size)
Illumina archaea V4-V5: depth-only Mantel r = 0.565 (p = 0.0001) [info]
ITS2: depth-only Mantel r = 0.431 (p = 0.0017) [info]
ONT bacteria: depth-only Mantel r = 0.519 (p = 0.0044) [info]

db-RDA residual diagnostics (permutest on the full model):
Illumina bacteria V3-V4: permutest full model F = 1.529, p = 0.0129; residual df = 5; full-model constrained eigenvalues = 0.352,0.199,0.130,0.094,0.075
Illumina archaea V4-V5: permutest full model F = 1.594, p = 0.0027; residual df = 5; full-model constrained eigenvalues = 0.256,0.224,0.118,0.082,0.066
ITS2: permutest full model F = 1.110, p = 0.2795; residual df = 5; full-model constrained eigenvalues = 0.662,0.334,0.253,0.219,0.191
ONT bacteria: permutest full model F = 1.043, p = 0.5536; residual df = 1; full-model constrained eigenvalues = 0.418,0.219,0.083,0.062,0.038
```

## B4 unbinned phyla

| bin | phylum / genus (GTDB) | completeness | contamination | GUNC | SED breadth ≥90% | identity mode | call | dsrAB | 16S link |
|---|---|---|---|---|---|---|---|---|---|
| SED_semibin_014 (SED) | Desulfobacterota / B33-G16 | 7.18 | 0.0 | 0.35/True | none (not run: no MIMAG-medium bin) | none | below MIMAG medium: not a MAG | none | none |
| SED_semibin_041 (SED) | Chloroflexota / none | 57.14 | 17.09 | 0.29/True | none (not run: no MIMAG-medium bin) | none | below MIMAG medium: not a MAG | none | none |
| SED_vamb_035 (SED) | Chloroflexota / RBG-13-51-36 | 10.34 | 0.08 | 0.33/True | none (not run: no MIMAG-medium bin) | none | below MIMAG medium: not a MAG | none | none |
| SED_reassembly_bin_398 (reassembly) | Chloroflexota / none | 34.45 | 6.64 | 0.62/False | none (not run: no MIMAG-medium bin) | none | below MIMAG medium: not a MAG | none | none |
| SED_reassembly_bin_405 (reassembly) | Desulfobacterota / none | 49.03 | 10.69 | 0.77/False | none (not run: no MIMAG-medium bin) | none | below MIMAG medium: not a MAG | none | none |

Bins: all existing workspace_2 bins of both libraries (raw MetaBAT2/VAMB/SemiBin2 + DAS_Tool union), then step-3 reassembly per library. The prompt's reproduction QC (re-run of the published binning) was not completed: the re-run failed at DAS_Tool and was replaced, by decision, with classification of the published bins themselves.

```
existing: 416 bins (SED 119, LIQ 297); Chloroflexota/Desulfobacterota bins: 3; of these MIMAG medium + GUNC pass: 0
step3/bins: 2 bins (SED 0, LIQ 0); Chloroflexota/Desulfobacterota bins: 2; of these MIMAG medium + GUNC pass: 0
none
none
LIQ: only 3 liquid-assembly contigs classified Chloroflexota/Desulfobacterota (1 + 2); 218 recruited reads (0.33 Mb) reassembled into 3 contigs (14.8 kb) -> nothing to bin. No Chloroflexota/Desulfobacterota genome is recoverable from the liquid library.
```

## B5 LIQ_contig_19299 vs MAG_022

| test | result | compatible? |
|---|---|---|
| linking reads LIQ_contig_19299 <-> MAG_022 contig (MAPQ>=20, >=1 kb each side) | 0 reads (19299 linked to any contig: 1 reads; see B5_links_summary.tsv) | no |
| positive control: MAG_022 contig-contig links | 40 contig pairs bridged by 78 reads | PASS |
| CheckM2 completeness without -> with | 68.03 -> 68.23 | info |
| CheckM2 contamination without -> with | 1.06 -> 0.95 | yes |
| GUNC CSS / pass without -> with | 0.16/True -> 0.16/True | info |
| duplicated single-copy markers (GTDB-Tk bac120) | newly duplicated bac120 markers: none; markers added: TIGR00337, TIGR00456 | yes |
| GC | 0.441 vs MAG_022 5-95% 0.381-0.439 | no |
| tetranucleotide distance to MAG_022 centroid | 0.0157 vs own contigs 95th pct 0.0176; nearest MAG by TNF: LIQ_das_metabat_022.fa | yes |
| coverage (LIQ reads, mean depth) | 11.4x vs MAG_022 5-95% 8.9-16.1x | yes |

| linking reads (n, MAPQ ≥ 20) | 0 reads (19299 linked to any contig: 1 reads; see B5_links_summary.tsv) |
|---|---|
| CheckM2 completeness / contamination without → with contig | 68.03 -> 68.23 / 1.06 -> 0.95 |
| duplicated markers | newly duplicated bac120 markers: none; markers added: TIGR00337, TIGR00456 |
| verdict: attributed / not attributed / undecided | sensitivity (B5_sensitivity.tsv): assembly graph has no links for contig_19299; relaxed linkage gives 3 MAG_022 reads against 16 background reads to other contigs -> verdict unchanged: not attributed |

