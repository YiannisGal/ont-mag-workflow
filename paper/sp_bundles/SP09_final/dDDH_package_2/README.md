# dDDH package 2: SeqCode-eligible medium-quality MAGs (GGDC 3.0, DSMZ)

Prepared 2026-10-22. Companion to `dDDH_package.tar.gz` (the 6 high-quality candidates); same procedure. These three MAGs are MIMAG medium quality but meet the SeqCode type-genome limits (>90% complete, <5% contamination).

| folder | MAG | GTDB genus | completeness/contamination | references |
|---|---|---|---|---|
| ETO_LIQ_MAG_005 | *Methanofollis* sp. | g__Methanofollis | 98.5 / 0.7 | 12 *Methanofollis* (6 type strains) + 2 Etoliko MAGs |
| ETO_LIQ_MAG_011 | *Methanofollis* sp. | g__Methanofollis | 90.5 / 2.2 | 12 *Methanofollis* (6 type strains) + 2 Etoliko MAGs |
| ETO_LIQ_MAG_008 | Lachnospiraceae bacterium | g__DIGJ01 | 94.9 / 2.4 | all 5 DIGJ01 species (no type material exists) |

The *Methanofollis* folders also include the other Etoliko *Methanofollis* MAGs as references (`ETO_LIQ_MAG_00x_Etoliko_MAG.fasta`): ETO_LIQ_MAG_007 from package 1, plus 005 or 011. This tests whether the three are one species or several. Those files are local only, so **upload them as FASTA** (they have no accession). `reference_accessions.txt` lists only the public references.

## How to run (https://ggdc.dsmz.de/ggdc.php)
1. Enter your e-mail address.
2. Query: `<MAG>_query.fasta`.
3. References: upload everything in `references/`. Alternatively, paste `reference_accessions.txt` **and** upload the `*_Etoliko_MAG.fasta` files.
4. One job per MAG (3 jobs).
5. Save each result as `<MAG>_ggdc.csv` and put it in `SP09_final/dDDH_results/` (with the 6 from package 1).

Use **formula 2** dDDH (robust for MAGs and draft genomes). 70% = species boundary.
