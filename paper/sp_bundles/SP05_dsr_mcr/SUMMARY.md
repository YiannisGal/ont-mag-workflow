# SP05 summary: identity of the unbinned dsrAB contigs, and the *Methanogenium* mcrA gap

**QC: PASS.**
- **dsr control:** in the DsrA, DsrB and concatenated DsrAB trees, the oxidative references (*Allochromatium*, *Chlorobium*) are split from all reductive references by one edge (UFBoot 100 in each tree; `QC_dsr.txt`).
- **mcrA control:** all 9 step-15 CONFIRMED proteins are reproduced (`QC_mcr.txt`).

**Deviations from the prompt:**
- No step-85 folder exists. The dsrAB references are step 58's 82 sequences, the set steps 58/59 used.
- No GTDB or NR protein set is local. The protein LCA uses UniProt 2026_03 (SwissProt + TrEMBL) with NCBI taxdump lineages.
- Trees: IQ-TREE 3.1.1, LG+C20+F+G, UFBoot 1000.

## H2 — who carries the two unbinned dsrAB operons?

Both operons are **reductive**. In every tree, all 5 nearest references are reductive sulfate reducers (`dsr_placement.tsv`).

| contig | nearest references (DsrA / DsrB / DsrAB) | UFBoot | protein LCA | lineage |
|---|---|---|---|---|
| LIQ_contig_12880 (6.4 kb) | *Desulfosporosinus* sp. DB / *D. meridiei* / *D. meridiei*; all 5 nearest are *Desulfosporosinus* or relatives | 100 (DsrA), 92 (DsrB) | genus *Desulfosporosinus* for all 5 proteins | ***Desulfosporosinus*** (Peptococcaceae / Desulfitobacteriia) |
| LIQ_contig_19299 (20.7 kb) | *Desulfofalx alkaliphila* / *Desulfotomaculum* sp. Lac2 / *D. alkaliphila*; the rest *Desulforamulus*, *Pelotomaculum* | 40–58 | Eubacteriales (NCBI); top hits *Desulforamulus*, *Desulfofundulus*, *Desulfocucumis*, *Pelotomaculum* | **Desulfotomaculia** clade (*Desulfofalx* / *Desulforamulus* group); genus unresolved |

For LIQ_contig_19299, 16 of its 21 proteins resolve to order Eubacteriales and 2 to class Clostridia; one resolves only to root, one only to kingdom, and one has no hit (`contig_protein_LCA.tsv`).

The spore-germination protein K06313 is **not** called above threshold on LIQ_contig_19299 (`work/K06313_check.txt`).

## H2 — could the operons belong to ETO_LIQ_MAG_022 (JBAIYS01, Desulfohalotomaculaceae)?

Tests are in `mag022_compatibility.tsv`.

**LIQ_contig_12880: no.** Its tetranucleotide profile is outside MAG_022's own-contig range, and its nearest MAG is LIQ_MAG_001. Its coverage (6.7×) is below MAG_022's 5–95% range (9.1–17.6×). It is a separate *Desulfosporosinus* organism.

**LIQ_contig_19299: three of four tests are compatible:**

| test | result |
|---|---|
| tetranucleotide profile | compatible; MAG_022 is the nearest of the 29 MAGs |
| coverage | compatible (11.4×) |
| marker redundancy | compatible; its only KOs are DsrA and DsrB, which MAG_022 lacks |
| GC | **incompatible**: 0.441, just above MAG_022's 95th percentile (0.438) |

The prompt requires all four, so the rule says **not attributable** to MAG_022. It stays a candidate. Wording for the thesis: "consistent with MAG_022 in composition, coverage and phylogenetic neighbourhood (Desulfotomaculia), but not attributed". The GC difference is small (0.003) and could be the operon's own composition. Re-binning with LIQ_contig_19299 forced into MAG_022 would settle it. This is weaker than step 58's "most parsimoniously MAG_022's missing operon".

Neither placement says anything about abundance or activity in situ. The contigs come from the enrichment assembly, and the in-situ signal is the 70/68 read-level counts.

## H3 — why is ETO_LIQ_MAG_017 (*Methanogenium*) missing from mcrA_verdicts.tsv?

The two contigs carrying its mcr genes (LIQ_contig_3836 and LIQ_contig_13851) are really in MAG_017's bin, so this is not the Vol XIII wrong-genome error.

**The reason is nomination, not phylogeny.** Each mcrA copy is split into two ORFs. Joined as step 54 does, both still score below the KOfam K00399 threshold of 775.53:
- LIQ_contig_3836_37+38: 759.9
- LIQ_contig_13851_5+6: 674.5

The fragments alone score 138–628 (`work/mag017_nomination.tsv`). So step 15 never nominated them.

**In the tree, both joined proteins are genuine methyl-coenzyme M reductase.** Placed in the step-15 mcrA tree, both are CONFIRMED MCR with the step-15 rule, and their nearest reference is MCR WP_460005979.1 (patristic distance 0.121 and 0.340). No ACR (alkane-oxidising) reference is nearest (`mag017_mcrA_verdict.tsv`).

**Wording for the thesis:** "MAG_017 carries two frameshift-fragmented mcrA loci that place with MCR. Joined, they fall just below the curated KOfam threshold, which is why automated nomination missed them. mcrA is present but fragmented, not absent." The 10th methanogen mcrA is supported by the tree, not by KOfam.

## Files

`dsr_placement.tsv`, `trees/` (dsrA, dsrB, dsrAB and mcrA treefiles and logs), `contig_protein_LCA.tsv`, `mag022_compatibility.tsv`, `mag017_mcrA_verdict.tsv`, `QC_dsr.txt`, `QC_mcr.txt`, `inputs_md5.tsv`, `versions.tsv`, `COMMANDS.log`.
