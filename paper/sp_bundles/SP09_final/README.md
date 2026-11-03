# SP09_final: novel taxa among the Etoliko MAGs

Start here: `SP09_SUMMARY.md` (status, results, QC).

| folder / file | content |
|---|---|
| `A1_gtdbtk/` | GTDB-Tk 2.7.2 classify_wf (r232) of the 8 HQ genomes; RED medians |
| `A2_refs/` | reference lists per candidate, GTDB r232 metadata/trees, downloaded genomes |
| `A3_ani/` | skani + FastANI raw outputs, `ani_summary.tsv` |
| `A4_aai/` | EzAAI pairs, `aai_summary.tsv`, `all_pairs_aai.tsv` |
| `A5_pocp/` | BLASTP outputs, `pocp.tsv` (POCP + POCPu) |
| `A6_trees/` | trees per candidate set (`*.nwk`, `*_tree.pdf/.svg`); runs in `v2/` |
| `A7_maps/` | Bakta, eggNOG, `*_circular_map.png/.svg` (pyCirclize), `*_cgview_map.png/.svg` (CGView), `proksee_upload/` |
| `A8_table/` | protologue-style summary table |
| `A10_combined_tree/` | tree of all 29 MAGs + 5 culturomics isolates + GTDB context (Figure 6 style) |
| `A11_species_panels/` | per-candidate tree + ANI/AAI/dDDH/ΔG+C panels (YGCC01 Figure 5 style) |
| `A12_MQ_ani_screen/` | skani screen of unnamed medium-quality MAGs vs GTDB r232; Etoliko all-vs-all |
| `A13_extra_candidates/` | ETO_LIQ_MAG_005, 008, 011 (SeqCode-eligible MQ MAGs) |
| `dDDH_package.tar.gz`, `dDDH_package_2.tar.gz` | inputs for GGDC (dDDH); put results in `dDDH_results/` |
| `versions.tsv` | all tools and databases |

The NCBI submission package is outside this folder: `workspace_2/NCBI_submission/` (with its own README).
