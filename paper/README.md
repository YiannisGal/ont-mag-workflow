# paper/ — analysis scripts for the Etoliko Lagoon metagenome study

This folder archives the scripts that produced the results reported in the
manuscript submitted to *Microorganisms*. It is a record, not a pipeline.

## Archived as run

Each script here is the version that produced the output used in the paper.
None has been refactored, and they will not run as-is on another machine:
they expect the original workspace layout, local databases and conda
environments. For a portable, re-runnable workflow, use the Snakemake
workflow at the top of this repository (see below).

`MANIFEST.md5` gives the md5 of every original file **before** sanitisation,
so each archived file can be matched to the file that ran. Scripts with an md5
stamp in their name (e.g. `47_master_2c00bba.sh`) were run through
`scripts/run.sh`, which refused to run a script whose contents did not match
that stamp.

## Layout

| path | contents |
|---|---|
| `scripts/` | final version of each analysis step, original filenames (`01_…` to `97_…`, figure and table scripts) |
| `scripts/superseded/` | earlier versions that the lab notebook discusses (e.g. step 16, replaced by 16b) |
| `sp_bundles/SP01…SP10_manuscript/` | the later, self-contained analysis bundles: scripts and their notes only, no outputs |
| `config/00_config.template.sh` | the workspace configuration, with every path replaced by a placeholder |
| `config/00_config_patch.sh`, `config/env_resolved.sh` | the environment patch and the resolved tool environments sourced by every step |
| `versions.tsv` | tool versions recorded at run time |
| `SCRIPT_INDEX.tsv` | one row per script (see below) |
| `SANITISATION_LOG.tsv` | every path substitution: file, line, original, replacement |
| `MANIFEST.md5` | md5 of every original file |

## How scripts map to the workspace

Every script sources `00_config.sh` and reads and writes under `${WS2}`, in a
folder named after its step number: `03_assemble.sh` writes `${WS2}/03_assembly/`,
`45_kofam_cf52d21.sh` writes `${WS2}/45_kofam_full/`, and so on. Steps with a
letter suffix (16b, 34b, 43k) add to or replace the folder of their parent
step. `${WS1}` is the first-pass workspace that holds the original assemblies
and liquid-culture reads.

## What the Snakemake workflow re-implements

The top-level `workflow/` covers the generic steps: assembly and polishing
(03–04, Flye + Medaka), read mapping (05), binning and consensus (06–07,
MetaBAT2, VAMB, SemiBin2, DAS Tool), bin quality (08, CheckM2 + GUNC),
taxonomy and phylogenomics (09, 12, GTDB-Tk + IQ-TREE), annotation (10, Bakta),
KOfam (11b, 45), read-level profiling (13, 22, 35: Kraken2, SingleM, Kaiju),
cross-mapping (26, CoverM), viruses (18, geNomad), CAZymes (30, dbCAN),
BGCs (33, antiSMASH), ANI (40, skani), and metabolic models and SMETANA
(42–43, gapseq + SMETANA). The study-specific checks, controls and targeted
gene analyses (steps 14–97 beyond those, and the SP bundles) exist only here.

## Reading SCRIPT_INDEX.tsv

| column | meaning |
|---|---|
| `step` | step number, or the SP bundle name |
| `script` | path under `scripts/` or `sp_bundles/` |
| `status` | `final` (produced the output on disk) or `superseded` |
| `what_it_does` | from the script's own header; caveats in `[brackets]` |
| `main_inputs` / `main_outputs` | workspace folders read and written |
| `lab_notebook_volume` | notebook volume(s) that discuss the step |
| `manuscript_table_figure` | supplementary figure (`S_…`) or table (`T_…`) built from this step's output; `—` when no figure or table uses it; `unknown` when it could not be determined |

The manuscript column was filled by matching each script's output folder and
files against the sources listed for each figure and table in the manuscript
package. It covers figures and tables only, not numbers quoted in the text.

Known caveats:

- Step 16 was never re-run. 16b replaces it.
- `15_marker_trees.sh` uses `wget`, which failed on the analysis server. Use `15_marker_trees_curl.sh` instead.
- 41f's `archaea_kofam_markers.tsv` has no score threshold.

## Sanitisation rule

The only edits to the archived scripts are mechanical replacements of absolute
machine paths with the variables that `00_config.sh` defines (`${WS2}`,
`${WS1}`, `${DB_DIR}`, `${TMPDIR}`, …), plus `echo`/`print` lines and comments
that only told the operator what to check or report next. No command,
parameter, threshold or file name was changed. Every path substitution is
listed in `SANITISATION_LOG.tsv`. Every script still passes `bash -n` or a
Python parse after substitution.
