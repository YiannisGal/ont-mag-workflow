## C3 unassigned vs depth
| dataset | measure | min % (depth) | max % (depth) | Spearman rho | p | Mann–Kendall tau | p | p (BH) |
|---|---|---|---|---|---|---|---|---|
| V3-V4 | phylum-unassigned | 0.05 (0-5) | 32.36 (45-50) | 0.345 | 0.2252 | 0.209 | 0.3244 | 0.4050 |
| V3-V4 | uOTU | 2.96 (65-70) | 43.45 (45-50) | 0.314 | 0.2700 | 0.297 | 0.1546 | 0.4050 |
| V4-V5 | phylum-unassigned | 0.49 (55-60) | 5.51 (35-40) | -0.187 | 0.5192 | -0.055 | 0.8267 | 0.5841 |
| V4-V5 | uOTU | 5.70 (15-20) | 26.30 (0-5) | -0.525 | 0.0562 | -0.407 | 0.0487 | 0.2529 |
| V4-V5 | DPANN share | 20.71 (65-70) | 66.22 (10-15) | -0.587 | 0.0280 | -0.429 | 0.0375 | 0.2520 |
| ITS2 | phylum-unassigned | 6.61 (45-50) | 81.66 (15-20) | -0.411 | 0.1450 | -0.319 | 0.1253 | 0.4050 |
| ITS2 | uOTU | 13.49 (45-50) | 87.38 (0-5) | -0.341 | 0.2368 | -0.297 | 0.1546 | 0.4050 |
| ONT14 | phylum-unassigned | 1.23 (65-70) | 6.50 (35-40) | 0.134 | 0.6443 | 0.121 | 0.5841 | 0.6443 |
| ONT14 | uOTU | 5.48 (65-70) | 25.45 (20-25) | 0.196 | 0.4988 | 0.077 | 0.7426 | 0.5841 |

## C1 ONT links, 14 depths
| genome | 16S copies | ONT14 cluster | identity | aln len | class | rank | peak depth (14-depth) | SP03 cluster | SP03 peak depth | changed? |
|---|---|---|---|---|---|---|---|---|---|---|
| ETO_LIQ_MAG_001 | 10 | C936 | 93.9 | 1487 | no_link | none | none | C315 | 65-70 | yes: class related->no_link; peak 65-70->none |
| ETO_LIQ_MAG_004 | 14 | C779 | 96.5 | 1527 | related | 172 | 45-50 | C127 | 10-15 | yes: class same_genus_level->related; peak 10-15->45-50 |
| ETO_LIQ_MAG_009 | 9 | C804 | 96.3 | 1471 | related | 342 | 60-65 | C284 | 60-65 | yes: class same_genus_level->related |
| ETO_LIQ_MAG_012 | 5 | C804 | 97.8 | 1472 | same_genus_level | 342 | 60-65 | C284 | 60-65 | yes: class exact->same_genus_level |
| ETO_LIQ_MAG_015 | 5 | C854 | 97.8 | 1457 | same_genus_level | 187 | 60-65 | C305 | 60-65 | no |
| ETO_LIQ_MAG_024 | 1 | C805 | 96.5 | 200 | related | 295 | 60-65 | C284 | 60-65 | no |
| ETO_LIQ_MAG_025 | 2 | C775 | 98.0 | 1473 | same_genus_level | 209 | 25-30 | none | none | yes: class no_hit_>=90%->same_genus_level; peak none->25-30 |
| ETO_SED_MAG_002 | 5 | C804 | 97.8 | 1471 | same_genus_level | 342 | 60-65 | C284 | 60-65 | no |
| YGCA01_C | 9 | C365 | 99.9 | 1441 | exact | 392 | 60-65 | C231 | none | yes: class no_link->exact; peak none->60-65 |
| YGCA01_Tp | 11 | C894 | 93.8 | 1443 | no_link | none | none | C99 | 60-65 | yes: class related->no_link; peak 60-65->none |
| YGCA07 | 19 | C894 | 98.2 | 1439 | same_genus_level | 12 | 45-50 | C99 | 60-65 | yes: class exact->same_genus_level; peak 60-65->45-50 |
| YGCA10 | 17 | C894 | 98.2 | 1439 | same_genus_level | 12 | 45-50 | C99 | 60-65 | yes: class exact->same_genus_level; peak 60-65->45-50 |
| YGCC01 | 9 | C798 | 94.9 | 1496 | related | 41 | 5-10 | C141 | 5-10 | no |

| genomes with no link | 21 |
|---|---|
| ETO_LIQ_MAG_002 | no_hit_>=90% |
| ETO_LIQ_MAG_003 | no_hit_>=90% |
| ETO_LIQ_MAG_005 | no_hit_>=90% |
| ETO_LIQ_MAG_006 | no_hit_>=90% |
| ETO_LIQ_MAG_007 | no_hit_>=90% |
| ETO_LIQ_MAG_008 | best hit C678 at 92.9% (below 94.5%) |
| ETO_LIQ_MAG_010 | no_hit_>=90% |
| ETO_LIQ_MAG_011 | unlinkable (no 16S) |
| ETO_LIQ_MAG_013 | best hit C517 at 93.0% (below 94.5%) |
| ETO_LIQ_MAG_014 | no_hit_>=90% |
| ETO_LIQ_MAG_016 | no_hit_>=90% |
| ETO_LIQ_MAG_017 | no_hit_>=90% |
| ETO_LIQ_MAG_018 | no_hit_>=90% |
| ETO_LIQ_MAG_019 | unlinkable (no 16S) |
| ETO_LIQ_MAG_020 | unlinkable (no 16S) |
| ETO_LIQ_MAG_021 | unlinkable (no 16S) |
| ETO_LIQ_MAG_022 | no_hit_>=90% |
| ETO_LIQ_MAG_023 | unlinkable (no 16S) |
| ETO_SED_MAG_001 | no_hit_>=90% |
| ETO_SED_MAG_003 | best hit C24 at 91.5% (below 94.5%) |
| ETO_SED_MAG_004 | best hit C804 at 92.3% (below 94.5%) |

| QC line |
|---|
| [PASS] Arc_Otu22 full-length gene vs ONT14: no hit >=90% (must be no hit, no error) |
| [PASS] self-identity: 156 of 156 genome 16S copies (all isolates + MAGs) hit themselves at 100% |
| [PASS] ONT14: 68 libraries, 14 depths (Results_new/final_output/cluster_tab_norm.txt; 2 of the 70 Metadata_minion barcodes absent: barcode29,barcode42) |
| [PASS] 0-5 cm: no top-5 ONT14 cluster at 0-5 cm links (>=97%) to a methanogen genome (C945=d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Mesobacillus;s__Bacillus_subterraneus / C950=d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Mesobacillus;s__Bacillus_jeotgali / C933=d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Sutcliffiella;s__Bacillus_sp. / C936=d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Niallia;s__Bacillus_oceanisediminis / C951=d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Niallia;s__uncultured_Bacillus); methanogen-linked clusters peaking at 0-5 cm: none |

## C2 community ~ geochemistry
| dataset | n | retained (VIF order) | selected (ordiR2step) | adj R² | F | p | p (BH) | Condition(depth) adj R² | p |
|---|---|---|---|---|---|---|---|---|---|
| V3-V4 | 14 | LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate (dropped: Polysaccharides(11.4)) | Phenols | 0.196 | 1.45 | 0.0808 | 0.1616 | -0.006 (terms: Phenols) | 0.4754 (BH 0.4754) |
| V4-V5 | 14 | LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate (dropped: Polysaccharides(11.4)) | TKN,NH4N | 0.285 | 1.74 | 0.0114 | 0.0456 | 0.221 (terms: TKN,NH4N) | 2.0e-04 (BH 8.0e-04) |
| ITS2 | 14 | LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate (dropped: Polysaccharides(11.4)) | Phenols | 0.104 | 1.22 | 0.1306 | 0.1741 | 0.037 (terms: Phenols) | 0.0584 (BH 0.0779) |
| ONT14 | 14 | LOI550,LOI950,Proteins,Phenols,NH4N,TKN,Sulfate (dropped: Polysaccharides(11.4)) | TKN | 0.092 | 1.19 | 0.2092 | 0.2092 | 0.107 (terms: TKN) | 0.0068 (BH 0.0136) |

| varpart (adj R²) | depth only | chemistry only | shared | residual |
|---|---|---|---|---|
| V3-V4 (chem = Phenols) | 0.007 | -0.006 | 0.157 | 0.842 |
| V4-V5 (chem = TKN,NH4N) | 0.125 | 0.221 | -0.008 | 0.662 |
| ITS2 (chem = Phenols) | 0.007 | 0.037 | 0.066 | 0.890 |
| ONT14 (chem = TKN) | 0.084 | 0.107 | -0.023 | 0.831 |
