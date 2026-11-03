# Caption drafts and legends: SP10 manuscript figures

Every figure is drawn without a title or caption; the text below is meant for the manuscript. Colours used throughout: **Bacteria #1F5FA6 (blue), Archaea #D9622B (orange)**; genomes from this study are in **bold**. MAG IDs are the final `ETO_SED_MAG_###` / `ETO_LIQ_MAG_###` IDs; inside figures the `ETO_` prefix is dropped (`LIQ_MAG_011` = `ETO_LIQ_MAG_011`). All paths are relative to `workspace_2/`.

## E1: re-drawn thesis figures

### S_MAG_quality (thesis fig01; script `fig01_mag_quality.py`)
**Caption.** Quality of the 29 MAGs. Completeness and contamination were estimated with CheckM2 v1.1.0. Dashed lines mark the MIMAG thresholds: high quality is ≥ 90% complete and < 5% contaminated (shaded box); medium quality is ≥ 50% complete and < 10% contaminated (dotted line at 10%).
**Legend.** x = completeness (%); y = contamination (%). Marker area = genome size (key: 1.6 / 3.5 / 5.5 Mb). Colour = domain. A black outline marks MAGs that meet the full MIMAG high-quality standard (23S, 16S and 5S rRNA plus ≥ 18 tRNA types).
**Source.** `47_master/MASTER.tsv` (completeness, contamination, bp, mimag, domain).

### S_function_heatmap (thesis fig03; script `fig03_function_heatmap.py`)
**Caption.** Functional gene sets in the 29 MAGs. Each cell gives the number of genes in a set that were found in a genome (KOfam above the curated thresholds; dbCAN3 for CAZymes). Cell colour is the count scaled to the column maximum, so columns can be compared with each other; it does not show pathway completeness. Empty cells mean zero genes. Genomes are ordered by domain and then by decreasing sediment breadth.
**Legend.** Rows = MAGs, with a domain strip on the left. Columns = gene sets (mcr, methylamine MT, methanol MT, fwd/fmd, hydrogenase, FDH, dsr/apr/sat, nif, ack/pta, B12 synthesis, B12-dependent, B12 transport). Colour bar = gene count relative to the column maximum.
**Source.** `47_master/MASTER.tsv`, the columns the thesis script reads. These columns are derived from `45_kofam_full/tables/function_panels.tsv` and `30_dbcan/family_matrix.tsv`. **Note:** the thesis script reads MASTER.tsv, not the two panel files named in the SP10 prompt; the script was kept unchanged so the numbers do not change.

### S_exchange_matrix (thesis fig06; script `fig06_exchange_matrix.py`)
**Caption.** Predicted metabolite exchange between the in-situ genomes (SMETANA, pairwise simulations on the sediment medium ETO_SED). Each cell is the summed SMETANA score of all compounds transferred from the donor (row) to the receiver (column). Grey cells are true zeros (no predicted transfer); the hatched diagonal is not defined. The bars give the total each genome donates (right) and receives (bottom). Both axes share the same order, by decreasing total received. Genome completeness is printed beside each row. Genomes < 90% complete are marked ° and shown in purple, because gap-filling may have supplied their needs, so an empty row or column for them is not interpretable. The scores describe potential, not activity.
**Legend.** Colour scale = summed SMETANA score per pair (0–10.8). Domain strips run along both axes.
**Source.** `67_network/fig06_matrix.tsv`, which is all 610 rows of `43_smetana/run_43q/pairs/pairs_detailed.tsv` summed per donor–receiver pair. The re-sum was checked: 610 rows → 32 pairs, 11 genomes, maximum 10.78, matching the matrix. Completeness is from `47_master/MASTER.tsv`.

### S_compound_classes (thesis fig07; script `fig07_what_is_exchanged.py`)
**Caption.** Compounds exchanged in the SMETANA predictions. The panel shows the twelve compounds with the highest summed SMETANA score across all pairs. Below them, separately and on the same linear scale, are the four substrates discussed in the text (H₂, formate, acetate, methanol). The axis is linear deliberately: on a log scale the contrast between amino acids and the energy substrates would disappear.
**Legend.** x = summed SMETANA score over all pairs. Bar colour = compound class (amino acids; energy and C1; other: nucleosides, organic acids, vitamins).
**Source.** `67_network/compound_summary.tsv` (derived from `66_network/`, `67_network/`).

### S_biogeography (thesis fig09; script `fig09_biogeography.py`)
**Caption.** None of the 29 MAGs is recovered outside Aitoliko. Breadth of coverage is the fraction of a genome covered by at least one read. Filled circles show the higher of the two Aitoliko libraries (sediment or culture). Open circles show the highest breadth in any of 30 public metagenome runs (15 BioProjects; 2,000,000 spots per run streamed with `fastq-dump -X`, i.e. the first spots, not a random sample). Open circles are plotted at their real values (0.000–0.007), not at zero. The dashed line is the detection threshold (0.30).
**Legend.** x = breadth of coverage (fraction of genome); rows = MAGs with genus.
**Source.** `14_comparative/recruitment_matrix.tsv`; `47_master/MASTER.tsv`.

### S_dsr_binning (thesis fig10; script `fig10_dsr_binning.py`)
**Caption.** Tetranucleotide test for the two unbinned contigs carrying dsrAB. Each contig's tetranucleotide profile is compared with every genome in the catalogue; the ten nearest genomes are shown, and the range of the rest is given at the top right. The grey bar is each genome's own yardstick: the 90th percentile of distances of its own contigs > 5 kb to its profile. A point inside the bar is compatible with membership; a point outside is not. Genomes with too few long contigs get no yardstick (square).
**Legend.** x = Euclidean tetranucleotide distance to the genome profile. Filled = inside (compatible), open = outside (incompatible), square = no yardstick.
**Source.** `57_dsr_bin/bin_paths.tsv`, `57_dsr_bin/dsr_contigs.fasta` (as in the thesis; coverage from `58_dsr_finish/coverage.tsv` and `34_unbinned/` is reported in the text, not drawn).

### S_ani_isolates (thesis fig11; script `fig11_ani_isolates.py`)
**Caption.** Whole-genome identity (skani ANI) between the MAGs, the cultured isolates and public reference genomes, computed all-against-all in one run. The boxes count, for each class of pair, how many pairs were compared and how many returned a value. skani reports a pair only above an internal similarity floor, so a missing value means the pair was computed and fell below it. No MAG–isolate pair returned a value. The lower panel shows every pair that returned a value and acts as a positive control for the method.
**Legend.** Pair classes: MAG × isolate, MAG × reference, isolate × isolate, isolate × reference, reference × reference. Dashed line = species threshold (95% ANI).
**Source.** `40_integration/skani_raw.tsv`, `mag_list.txt`, `iso_list.txt`.

### S_mcrA_nomination (thesis fig14; script `fig14_15_mcra.py`)
**Caption.** Every K00399 (mcrA) hit in each archaeal MAG, scored against the curated KOfam full-sequence threshold (775.53). Filled circles are ORFs above the threshold, which were nominated for the mcrA tree. Open circles are below it. ETO_LIQ_MAG_017 (*Methanogenium*) carries four K00399 fragments (138–628), all below the threshold, because its mcrA genes are split across consecutive ORFs; see S_mcrA_tree for the joined proteins.
**Legend.** x = KOfam score (full sequence); dashed line = curated threshold. Beside each MAG: sediment breadth and the number of its proteins in the tree.
**Source.** `45_kofam_full/all_hits.tbl` (K00399 rows, joined to MAGs by bin name); sediment breadth from `47_master/MASTER.tsv`. The per-MAG scores agree with `90_mag_bin_mapping/mcr_scores_by_mag.tsv` (e.g. MAG_003 993.9, MAG_005 1028.9). `54_mcr_fragmentation/hits_with_coords.tsv`, which has permuted labels, was not used.

### S_phylogenomics_archaea / S_phylogenomics_bacteria (thesis fig12a/b; script `fig12_phylogenomics.py`)
**Caption.** Phylogenomic placement of the MAGs among GTDB reference genomes. These are maximum-likelihood trees (IQ-TREE) of concatenated GTDB-Tk marker proteins (ar53 for Archaea, bac120 for Bacteria). Bold tips are MAGs from this study, with their completeness. Triangles are collapsed reference-only clades, with the number of tips beside them. Ultrafast bootstrap support is shown where ≥ the threshold printed in the panel. Archaea: rooted on the outgroup. Bacteria: the outgroup is not monophyletic in this tree, so **no root is drawn**: the tree is shown as IQ-TREE wrote it and the root position must not be read.
**Legend.** Scale bar = substitutions per site.
**Source.** `12_phylogeny/ref_trees/ar53_ref.contree`, `bac120_ref.contree` (+ `*_labels.tsv`, `*_outgroup.txt`, `*.iqtree`); MASTER.tsv for genome identity.

### S_mcrA_tree (thesis fig15; new plot `plot_sp05_trees.R`)
**Caption.** Maximum-likelihood tree of McrA proteins (SP05 IQ-TREE 3.1.1 run, ModelFinder, UFBoot and SH-aLRT 1,000 each). It includes 11 sequences from this study: 9 KOfam-nominated McrA proteins from seven archaeal MAGs (two each in ETO_LIQ_MAG_005 and ETO_LIQ_MAG_007), plus the two joined McrA proteins of ETO_LIQ_MAG_017. All 11 are classed as MCR (SP05 `QC_mcr.txt`: 9/9 nominated proteins CONFIRMED; both MAG_017 proteins CONFIRMED, nearest MCR WP_460005979.1); none falls within the ACR clades. The two bacterial restriction-endonuclease "McrA" sequences (WP_484489188.1, WP_471601478.1) were pruned from the drawing. ACR α subunits, the other ACR subunits (β/γ) and component A2 are coloured separately. The ACR α sequences are not monophyletic in this tree, so no root is drawn (as written by IQ-TREE). UFBoot is shown where ≥ 70.
**Legend.** Orange bold = this study; brown = MCR α; blue = ACR α; light blue = ACR other subunits; grey = ACR component A2. Scale = substitutions per site.
**Source.** `99_thesis_audit/SP05_dsr_mcr/trees/mcrA.treefile`, `mcrA_in.faa`; subunit categories from the reference descriptions in `15_marker_trees/mcrA/refs/{acr,mcr}.faa` (same rule as `fig14_15_mcra.py`).

### S_dsr_trees (thesis fig13; new plot `plot_sp05_trees.R`)
**Caption.** Maximum-likelihood trees of (**a**) DsrA and (**b**) DsrB (SP05, IQ-TREE 3.1.1, LG+C20+F+G, UFBoot 1,000; consensus trees). The sequences from this study come from the two unbinned contigs LIQ_contig_12880 and LIQ_contig_19299. All of them fall among reductive-type (sulfate-reducing) DsrAB (5/5 nearest references reductive; SP05 `dsr_placement.tsv`). The oxidative-type rDsr references (*Chlorobium*, *Allochromatium*) are separated from the reductive references by one edge with UFBoot 100 in both trees (SP05 `QC_dsr.txt`). The trees are drawn as IQ-TREE wrote them (no root is implied). UFBoot is shown where ≥ 70.
**Legend.** Orange bold = this study; green = reductive type; ochre = oxidative type (rDsr). Scale = substitutions per site.
**Source.** `99_thesis_audit/SP05_dsr_mcr/trees/dsrA.contree`, `dsrB.contree`, `dsrA_in.faa`, `dsrB_in.faa`.

### S_rrna_corrected (step 27; new plot `fig_new_panels.py`)
**Caption.** Effect of rRNA operon copy number on read-based genus shares in (**a**) the sediment and (**b**) the enrichment cultures. Grey bars show raw shares of SSU reads assigned to the five genera with a MAG; coloured bars show the shares after dividing by each genus's rRNA operon count (in parentheses, from its MAG). *Lysinibacillus* (13 operons) is the most inflated in the raw profile.
**Legend.** x = % of genus-assigned SSU reads within this subset; colour of the corrected bar = domain.
**Source.** `27_rrna_corrected/SED_corrected_profile.tsv`, `LIQ_corrected_profile.tsv`.

### S_arc_otu22 (Arc_Otu22 full-length 16S tree; copied from the delivered Figure 7)
**Caption.** (**a**) Maximum-likelihood tree of 145 16S rRNA gene sequences over 1,372 aligned positions (MAFFT v7.505 L-INS-i; trimAl v1.5.rev1 gappyout; IQ-TREE 2.3.6, GTR+F+I+R4). Numbered nodes give SH-aLRT/UFBoot (key in panel) and are shown where support is ≥ 80/≥ 95. The four non-Asgard sequences are not a clade, so the tree is midpoint-rooted. The bracket marks the smallest clade containing Arc_Otu22 and both Thorarchaeota MAG 16S genes (45 tips). Open circles are the 26AF clones placed in Thorarchaeota by Zhang et al. (**b**) The Illumina amplicon aligned to the 1,562-bp gene recovered from a single shotgun read (positions 481–932): 411 of 415 template bases are identical.
**Source.** `98_arc22_longonly/DELIVERY_Arc_Otu22_2026-09-27/1_figures/Figure_7_Arc_Otu22.{pdf,png}` (already English and caption-free; copied unchanged). Full caption: `2_captions/CAPTIONS.md`. Width is 16.0 cm, not 17.

### S_viruses_defence (steps 18/25/28; new plot `fig_new_panels.py`)
**Caption.** Mobile-element and defence content per MAG. Each bar stacks three counts: proviruses detected by geNomad, defence systems (all types), and CRISPR arrays. MAGs are grouped by domain and sorted by total.
**Legend.** Purple = proviruses; green = defence systems; yellow = CRISPR arrays. x labels coloured by domain.
**Source.** `18_viruses/provirus_host_links.tsv` (one row per provirus), `28_crispr/defence_systems.tsv` (`total`), `28_crispr/spacers.fasta` (arrays = distinct array IDs in the headers).

### S_bgc (step 33/33b; new plot `fig_new_panels.py`)
**Caption.** Biosynthetic gene cluster regions per MAG (antiSMASH; 141 regions in total). Regions are grouped into broad classes: terpene; RiPP; NRPS/PKS; hybrid (regions whose products span more than one class); other (e.g. cyclic-lactone autoinducer, deazapurine, β-lactone, ectoine). Filled segments are regions with at least one MIBiG hit (33 of 141); open segments have none.
**Legend.** Colour = class; filled vs open = with vs without a MIBiG hit. x labels coloured by domain.
**Source.** `33_antismash/mibig_hits.tsv` (`products`, `n_mibig_hits`). The class grouping is a display grouping defined in `fig_new_panels.py`; the per-region product strings are in `E2_tables/T_antismash.tsv`.

## E3: supplementary trees (re-drawn without captions; see `../E3_trees/`)
- `redrawn/trees/*_tree`: per-candidate GTDB-Tk marker trees (SP09 A6/A13). Node labels give SH-aLRT/UFBoot, shown where ≥ 80/≥ 95; ᵀ = type strain; bold coloured tips = Etoliko MAGs.
- `redrawn/species_panels/*`: tree + ANI / AAI / dDDH / ΔG+C. dDDH from TYGS, GGDC formula 2 (d4); n.d. = pair not compared; "<80" = ANI below the skani/FastANI detection limit; ΔG+C = reference minus MAG.
- `redrawn/combined_tree/*`: all MAGs + 5 isolates + GTDB r232 context. GTDB-Tk 2.7.2 bac120/ar53, IQ-TREE 3.1.1 (LG+F+G4). The bacterial and archaeal trees are midpoint-rooted and joined at the root for display only. Grey triangles are collapsed reference clades without Etoliko genomes.
