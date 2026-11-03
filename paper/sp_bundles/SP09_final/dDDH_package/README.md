# dDDH package for the SP09 candidate MAGs (GGDC 3.0, DSMZ)

Prepared 2026-10-22. One folder per candidate MAG. These are the same reference panels that will appear in the species-level tree figures (`A11_species_panels`, same style as YGCC01 Figure 5).

| folder | MAG | references | note |
|---|---|---|---|
| ETO_SED_MAG_001 | *Sedimentibacter* sp. | 28 | 3 type strains + 25 nearest by AAI |
| ETO_LIQ_MAG_001 | Cytobacillaceae (candidate new genus) | 28 | one per Cytobacillaceae genus (type species preferred) + nearest by AAI |
| ETO_LIQ_MAG_003 | MTP4 sp. | 6 | all GTDB MTP4 species; no type material exists |
| ETO_LIQ_MAG_007 | *Methanofollis* sp. | 12 | all; 6 type strains |
| ETO_LIQ_MAG_009 | *Chakrabartyella* sp. | 16 | all; 4 type strains |
| ETO_LIQ_MAG_010 | *Methanosarcina* sp. | 28 | 13 type strains + nearest by AAI |

Each folder contains:
- `<MAG>_query.fasta`: the MAG (query genome)
- `references/<accession>_<species>.fasta`: reference genomes (from NCBI, as used in SP09)
- `reference_accessions.txt`: the same references as NCBI accessions, one per line
- `reference_list.tsv`: species, genus, type status and AAI to the MAG

## How to run (https://ggdc.dsmz.de/ggdc.php)
1. Enter your e-mail address. Results arrive by e-mail, usually within hours.
2. **Query genome:** upload `<MAG>_query.fasta`.
3. **Reference genomes:** either upload the FASTA files from `references/` or paste the accessions from `reference_accessions.txt` (GGDC accepts NCBI accession numbers). Up to 75 references per job, so one job per MAG is enough.
4. Submit one job per MAG (6 jobs).
5. Download the **CSV/XLSX result** from the e-mail link and save it as `<MAG>_ggdc.csv` in `SP09_final/dDDH_results/`.

## What is used
**Formula 2** dDDH (recommended by DSMZ for draft genomes and MAGs), with its confidence interval, plus the probability of DDH ≥ 70% and the G+C difference GGDC reports. The 70% cut-off defines a species. Each value goes into the dDDH column of the species-level figures and into the A8 protologue table.

Note: these are MAGs (all ≥ 90% complete, ≤ 5% contamination). Formula 2 is the one that is robust to incomplete genomes. Use it, not formulas 1 or 3.
