# SP04 summary: order-level classification of archaeal V4–V5 OTUs

**Method:** 2,544 archaeal OTUs (the full normalised table, primers removed) were classified two ways against GTDB SSU R232:
- vsearch 2.32 `--sintax`, bootstrap cut-off 0.8;
- BLAST top hit, accepted at ≥90% identity over ≥90% of the OTU.

GTDB file md5 8393492c…, matching the official MD5SUM. Depth aggregation follows ms v2: extractions summed per depth, then relative abundance. Archaeal depths use the matching YG number (user decision, SP03).

**QC: PASS.** Arc_Otu22 falls in Asgardarchaeota / Thorarchaeia / Thorarchaeales (sintax confidence 1.00; BLAST 98.6%), not Thermoplasmatota. Reads unclassified at order level make up 1.7–9.7% per depth (`QC.txt`).

## Are Methanomassiliicoccales in the amplicons, and where?

Yes, but at trace level. They occur at every depth, at 0.01–0.02% of archaeal reads at 0–15 cm and 0.05–0.17% from 15 cm down, with the maximum at 25–30 cm (0.17%). Below 15 cm they are the most abundant methanogen order in most layers.

## Which methanogen orders dominate the deep layers?

None. All methanogen orders together make up less than 0.5% of archaeal amplicons at every depth (`order_depth_profiles.tsv`):

| order | share of archaeal reads | where |
|---|---|---|
| Methanomassiliicoccales | 0.05–0.17% | 15–65 cm |
| Methanofastidiosales | up to 0.15% | 25–30 cm |
| Methanosarcinales | 0.17% | 0–5 cm |
| Methanomicrobiales | 0.10% | 0–5 cm |
| Methanobacteriales | ≤ 0.08% | — |

The archaeal library is dominated by non-methanogens:
- **Anstonellales** (DPANN): 16–63%
- **DHVEG-1** (Thermoplasmata, Marine Benthic Group D-like): 5–44%
- **Nitrososphaerales**: 3–31%
- **Sigyniarchaeales** and **Woesearchaeales** (Asgard / DPANN)
- **Bathyarchaeales**: up to 9% at 55–60 cm

"Deep-layer methanogens are affiliated with methylotrophic lineages" can only be said in a narrow form: among a very small methanogen amplicon fraction, the methyl-reducing Methanomassiliicoccales are the most frequent order below 15 cm. Taxonomy does not establish methylotrophy for any other order. This is consistent with the metagenome, where Methanomassiliicoccales appear at read level only. It does not contradict the hydrogenotrophic *Methanofollis* MAG, which has no 16S gene and cannot be placed on the amplicon depth axis.

## The 25–30 cm "unassigned Thermoplasmata" spike (C47)

Thermoplasmata rise from 6–12% of archaeal reads at 0–20 cm to **46.9% at 25–30 cm**, then stay at 19–27% down to 60 cm (`thermoplasmata_25_30cm_breakdown.tsv`).

At 25–30 cm, **92.6% of the Thermoplasmata signal is order DHVEG-1** (60 OTUs, 43.5% of all archaeal reads). This is a non-methanogenic Marine Benthic Group D-type lineage. Methanomassiliicoccales contribute 0.4%.

A further 1.5% comes from 11 OTUs that SILVA labels Thermoplasmata but GTDB places in the bacterial order Enterobacterales. These are classifier-conflicting OTUs and are listed in `archaeal_OTU_orders.tsv`.

## Links to methanogen MAGs (from SP03, ≥97%)

| MAG | OTU | identity | OTU order | OTU total reads |
|---|---|---|---|---|
| MAG_003 (MTP4) | Otu448 | 100% (exact) | Methanosarcinales | 21 |
| MAG_010 (*Methanosarcina*) | Otu448 | 97.2% | Methanosarcinales | 21 |
| MAG_010 (*Methanosarcina*) | Otu909 | 98.7% | Methanosarcinales | 46 |
| MAG_016 (*Methanobacterium_C*) | Otu1839 | 97.2% | Methanobacteriales | 2 |

MAG_018 (*Methanosarcina*) links exactly (100%) to Otu909 (46 reads); see `mag_links_ge97.tsv`.

All the linked OTUs are rare. The *Methanofollis* MAGs (005, 007, 011) and *Methanogenium* MAG_017 have no link at ≥97%.

## Limits

- Order membership from a ~415-bp region is L1.
- Only Methanomassiliicoccales can be called "affiliated with methyl-reducing methanogens".
- Shares are amplicon shares, with primer and copy-number bias.

## Files

`archaeal_OTU_orders.tsv` (OTUs with ≥100 reads), `order_depth_profiles.tsv`, `thermoplasmata_25_30cm_breakdown.tsv`, `mag_links_ge97.tsv`, `QC.txt`, `versions.tsv`, `COMMANDS.log`.
