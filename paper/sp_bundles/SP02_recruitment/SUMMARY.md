# SP02 summary: identity-filtered competitive read recruitment

**Status:** QC PASS (`QC.txt`).

The run uses two competitive indices:
- **Index A:** the 29 MAGs, plus YGCA01_Clean, YGCA01_Contaminant, YGCA10 and YGCC01, plus the 11 background genomes from `all_12_genomes.txt`.
- **Index B:** the same, with YGCA07 in place of YGCA10.

Reads mapped with minimap2 map-ont; coverage computed with CoverM at read identity ≥90/95/97/99% (≥50% of each read aligned) and with no filter. 4.47% of the 6,827,767 sediment reads map (expected ~2–6%).

## Calibration (user decision, 2026-09-29, option 2)

The spec rule (breadth ≥ 0.50 at ≥95% identity and histogram mode ≥ 97%) failed its own positive control. Five of the in-situ ten MAGs lost 0.04–0.61 breadth at 95%.

Read-identity histograms show why:
- The sediment's own populations (SED_MAG_001–004) peak at 98–99% identity. This gives REF_MODE = 99, the identity of a same-population read in this Nanopore library.
- Those populations still keep only 73–94% of their reads at ≥95%. The 95% filter therefore removes true same-population reads, and low-depth genomes lose breadth.

The calibrated rules are:

| call | rule |
|---|---|
| **Same-species present** | histogram mode ≥ REF_MODE − 2 (97) **and** breadth ≥ 0.50 at ≥90% |
| **Related population only** | breadth ≥ 0.10 at ≥90% with a lower mode |
| **Not detected** | anything else |

The spec rule is also reported for every genome (`spec_call` in `isolate_calls.tsv`).

**Control 1 was re-tested unfiltered, as step 26 was run.** All 10 in-situ MAGs reproduce their workspace_2 breadth within 0.01.

## QC

**Positive control 2 (*Chakrabartyella*):** reads split between SED_MAG_002 (breadth 0.944, 8.4% of mapped reads) and LIQ_MAG_012 (0.954, 9.4%) at ≥95%.

**Negative control:** MAG_016 reaches at most 0.004 breadth at any threshold, under the 0.02 limit.

**Finding:** three culture-derived "in-situ" MAGs peak at only 93–94% identity in the sediment reads, so they classify as related populations only:

| MAG | sediment histogram mode | breadth at ≥90% | breadth at ≥95% |
|---|---|---|---|
| LIQ_MAG_011 (*Methanofollis*) | 94% | 0.66 | 0.28 |
| LIQ_MAG_015 | 94% | 0.82 | 0.49 |
| LIQ_MAG_025 | 93% | 0.83 | 0.51 |

The sediment carries their relatives, not the cultured genomes. Their "in situ" label should be qualified.

## Isolates (sediment library; values in `isolate_calls.tsv`)

**YGCA01_Contaminant — same-species present (calibrated).**
- Histogram mode 98%. Breadth 0.52 at ≥90%, 0.46 at ≥95%. Depth 0.76×. 1.23% of mapped reads.
- The strict spec rule would say related population only, because the 95% breadth is just below 0.50.

**YGCA01_Clean — related population only.**
- Histogram mode 96% (below 97). Breadth 0.51 at ≥90%, 0.42 at ≥95%. Depth 0.68×.
- This is a borderline case: the mode is 1% below the cut-off.

**YGCA10 — related population only.**
- Histogram mode 79%. Breadth 0.18 at ≥90%, 0.15 at ≥95%, 0.06 at ≥99%.
- The reads come from distant relatives of this *Paraclostridium*.

**YGCA07 (index B) — related population only.**
- Histogram mode 79%. Breadth 0.18 at ≥90%.
- These are the same numbers as YGCA10, as expected for a near-clone.

**YGCC01 — not detected.**
- Breadth 0.0008 at ≥90%, from 127 reads with modal identity 66%.
- The sediment library holds no DNA from this *Sporosarcina* above its detection floor.

**All isolates are detected as same-species in the culture library** (breadth ≥ 0.86 at ≥90%, modes 97–98%), except YGCC01 (not detected) and YGCA01_Contaminant (related population only).

**Coverage pile-ups:** a genome is flagged when its observed breadth is below half the Poisson expectation. No isolate is flagged at any threshold. Pile-ups appear only in the unfiltered runs:
- **LIQ_MAG_001 (Cytobacillaceae):** breadth 0.15 at 4.1× depth in the sediment, against 0.98 expected. Its sediment reads sit on local regions (such as rRNA or mobile elements), not across the genome.
- **LIQ_MAG_013** is also flagged.
- **Six Sporosarcina background genomes** are flagged, from the same kind of local stacking.

## Interpretation limits

- Recruitment shows that DNA similar to the reference is present at the stated identity. It says nothing about activity or viability.
- "Not detected" means below this library's detection floor (1 read in ~6.8 M), not absent from Etoliko.
- The library pools all depths, so it gives no depth information.
- Relative abundance is a share of mapped reads, not of cells.
- In the main run, background genomes cover only the relatives of YGCC01 (the `all_12_genomes.txt` set). The `bg_extended/` re-run adds 39 *Clostridium_J*, *Terrisporobacter* and *Paraclostridium* GTDB representatives. Calls are unchanged, and it shows that YGCA10/07's remaining reads are same-species DNA at very low coverage (see the last section).

## Provenance

- reads (md5 in `versions.tsv`) + index (md5 in `versions.tsv`; per-genome md5 in `refs/genomes.tsv`)
- → minimap2 command (`COMMANDS.log`) → `bam/`
- → CoverM filter/genome at each threshold (`cov/`)
- → `recruitment_by_identity.tsv` and `identity_histograms.tsv`
- → `isolate_calls.tsv` → this summary.

## Re-run with Clostridium_J / Terrisporobacter / Paraclostridium backgrounds (`bg_extended/`)

**Background set:** GTDB R232 species representatives of the three genera, downloaded from NCBI (43 of 45 obtained).
- **Excluded as the isolates' own species (fastANI ≥ 95%):**
  - *C. sulfidigenes* (97.1% to YGCA01_Clean)
  - *T. petrolearius* (95.8% to YGCA01_Contaminant)
  - *P. benzoelyticum* (98.8% to YGCA10)
  - *P. dentum* (95.2% to YGCA10)
- **Kept:** 39, added to the 11 Sporosarcina backgrounds (`bg_extended/bg/kept.tsv`).

**QC passes again:** the in-situ ten are reproduced unfiltered, MAG_016 stays at 0.004, and 4.75% of sediment reads map.

**Calls are unchanged; breadth barely moves:**

| isolate | sediment breadth at ≥90% | read-identity mode | reads left |
|---|---|---|---|
| YGCA01_Clean | 0.505 → 0.497 | 96% | — |
| YGCA01_Contaminant | 0.522 → 0.514 | 98% | — |
| YGCA10 | 0.184 → 0.158 | 79% → **98%** | 4,071 → 1,043 |
| YGCA07 | 0.177 → 0.150 | 79% → 98% | — |
| YGCC01 | not detected | — | — |

**What the backgrounds absorbed:** about 3,000 low-identity reads from YGCA10. They were *Paraclostridium* and *Clostridium_J* relatives, not YGCA10.

**What remains for YGCA10/07 looks like true same-species DNA at very low coverage.**
- The remaining ~1,000 reads peak at 98% identity.
- Their breadth equals the Poisson expectation for their depth at every threshold: 0.143 observed vs 0.155 expected at 0.17× depth (≥95%), and 0.054 vs 0.054 at ≥99%.
- So the reads are spread across the genome, not piled up.
- The calibrated rule still says "related population only", because 0.17× depth cannot reach 0.50 breadth.

**The two YGCA01 genomes** also follow the Poisson expectation (observed about 83–87% of expected). Their signals are spread across the genome too.

**Corrected reading for the thesis:**

| isolate | reading |
|---|---|
| YGCA01_Contaminant | same species present |
| YGCA01_Clean | borderline (mode 96%) |
| YGCA10/YGCA07 | same-species DNA present at very low coverage (about 0.2×, ~1,000 reads); below the breadth needed for a formal call |
| YGCC01 | not detected |

The background genomes themselves recruit little (breadth ≤ 0.05 at ≥90%). So the sediment's related populations are not the sequenced GTDB species either.
