# SP03 summary: 16S links from genomes to amplicon OTUs, with depth profiles

**Status:** QC passes (`QC.txt`), with two user decisions made on 2026-09-29.

- **Arc_Otu22 control accepted at 97.1%.** The primer-free Arc519F/915R region is 415 positions, the same as in ms v2. Otu22 is the top hit (407/417 aligned positions) and Otu674 the second (95.9%). The shortfall from ms v2's 411/415 is one cluster of about 8 mismatches at gene positions 841–850 in the supplied medaka-polished Nanopore gene.
- **Archaea in the bacterial V3–V4 dataset are treated as expected,** because 341F/805R also amplifies archaea. Methanogen genes match archaeal or unassigned OTUs there (Otu892, 2304, 2978, 3645, 2269, 806). These matches are listed in `QC.txt` and `links.tsv`. No methanogen matches an OTU that SILVA classifies as bacterial.
- **The 0–5 cm sanity check passes.** Its top five OTUs are four Bacillaceae and *Synechococcus*. None links to a methanogen, and none links at ≥97% to any genome.

## Other decisions and deviations

- The Illumina tables are the **full normalised** ones (`otutab_norm.txt`: 3,674 bacterial and 2,544 archaeal OTUs).
  - SILVA taxonomy exists only for the 1,025 + 395 trimmed OTUs, so linked OTUs outside that set carry no taxonomy.
  - Plastid exclusion therefore covers only assigned OTUs.
- **Archaeal sample depths** were taken from the bacterial sample with the same YG number, because no archaeal metadata file exists.
- **ONT** uses `Etoliko_138_2/ONT/final_output`: 26 samples rarefied to 1,000 reads. 9 of 14 depths are represented, with 1–5 samples each.
- **Replicates per depth are not always 5.** The metadata gives 6 extractions for 20–25 cm and 4 for 25–30 cm.
- **Primers were removed from both sides.** The OTU FASTAs still contain their primers, so they were trimmed before matching. The genome regions require both primers to be found.
- **Isolate 16S genes come from barrnap** on each isolate's assembly. The supplied `all_strains_16S.fasta` holds YGCC01 and Sporosarcina reference genomes, not the other four isolates. The supplied `YGCC01_16S.fasta` is 2 bp shorter than the barrnap copy.
- **vsearch 2.32.0** was used. The system version is 2.27, below the 2.28 the spec requires.

## Which genomes cannot be linked

- **No 16S gene** (reported as "unlinkable"): MAG_011 (*Methanofollis*), 019, 020 (*Lacrimispora*), 021 and 023.
- **Only a short fragment** (184 and 216 bp): MAG_005 and MAG_024, so table03's list of 7 holds in practice. MAG_024's fragment still reaches the "related" class against ONT.
- **A 16S gene but no link at ≥94.5%** in any dataset: MAG_002, 007, 014, 022, 025 and SED_MAG_001.

## Exact links (≥99.0% over ≥95% of the OTU; these are `figure_data_Fig6_2.tsv`)

| genome | dataset | OTU | OTU rank | depth where it is most abundant |
|---|---|---|---|---|
| YGCA07, YGCA10 (*Paraclostridium*) | V3–V4 | Otu1641 | 2544 | 5–10 cm |
| YGCA07, YGCA10 (*Paraclostridium*) | ONT | C99 (*P. bifermentans*) | 116 | 60–65 cm |
| LIQ_MAG_012 / SED_MAG_002 (*Chakrabartyella*) | ONT | C284 (*Anaerostignum*) | 242 | 60–65 cm |
| LIQ_MAG_018 (*Methanosarcina*) | V4–V5 | Otu909 | 328 | 0–5 cm |
| LIQ_MAG_003 (MTP4) | V4–V5 | Otu448 | 549 | 0–5 cm |
| LIQ_MAG_008 | V3–V4 | Otu2078 | 2713 | — |
| LIQ_MAG_010 / 018, LIQ_MAG_017 | V3–V4, archaea in the bacterial set | Otu892, Otu2304 / Otu806 | 773–2298 | — |

YGCC01 has no exact link. Its best match is ONT C141 (*Paenisporosarcina*, rank 14, 5.87% at 5–10 cm), but only at 95.0% ("related").

## What this means

**The linked amplicon populations are rare.**
- Only 12 genome–OTU pairs are exact, involving 10 genomes and 9 distinct OTUs.
- Almost all of those OTUs rank below 100 in their dataset. Most make up less than 0.1% of amplicons at their most abundant depth; the ONT clusters C99 and C284 reach about 0.7% and 1.3%.
- The dominant amplicon OTUs, such as the 0–5 cm Bacillaceae, have no genome counterpart.

**Depth placement is therefore thin.** Only these genomes can be placed on the depth axis: *Paraclostridium* YGCA07/10, the *Chakrabartyella* MAGs through ONT, and *Methanosarcina* MAG_018 and MTP4 MAG_003 through the archaeal primers.

**The Illumina and ONT profiles of one genome can disagree.** YGCA07 is most abundant at 5–10 cm in V3–V4 but at 60–65 cm in ONT. The ONT depth coverage is sparse, and each platform has its own primer and copy-number biases.

**Limits on interpretation.**
- A ≥99% match over about 415 bp supports the same genus, possibly the same species, within one core. It does not establish that the genome's population is the amplicon population.
- The values are amplicon shares, biased by rRNA copy number and primer choice.
- They say nothing about activity.

## Files and provenance

The chain for each value is: genome FASTA (md5 in `rrna/genomes.tsv`) → barrnap (`rrna/*.gff`) → region (`regions/`, cutadapt logs) → vsearch command (`COMMANDS.log`) → `hits/*.tsv` → `links.tsv` → OTU table (md5 in `inputs_md5.tsv`) → `depth_profiles.tsv` → `figure_data_Fig6_2.tsv`.
