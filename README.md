# ONT long-read metagenome → MAG workflow

A Snakemake workflow that reproduces the genome-resolved analysis of the
Etoliko Lagoon shotgun metagenomes (Galiatsatos et al., *Microorganisms*, in prep.)
and runs it on **anyone's Oxford Nanopore metagenomes**.

```
reads ─► merge + dedup ─► Flye --meta ─► Medaka ─► self + cross mapping (minimap2)
                                                        │
              MetaBAT2 ◄──── differential coverage ─────┼──► VAMB      SemiBin2
                 └──────────────────► DAS_Tool ◄───────────────┘           │
                                         │◄────────────────────────────────┘
          CheckM2 · GUNC · barrnap · tRNAscan-SE ─► MIMAG + filter ─► <PREFIX>_<SAMPLE>_MAG_<NNN>
                                         │
   ┌──────────────┬──────────────────────┼─────────────────────┬──────────────────────┐
 GTDB-Tk + RED   skani ANI     detection (CoverM identity    annotation (Prodigal,     phylogeny
                 (+ isolates)  sweep + read-identity mode)   KofamScan, dbCAN3,        (GTDB-Tk markers
                                                             geNomad, opt. antiSMASH,   → IQ-TREE)
 read profiling: SingleM · Kraken2/Bracken · Kaiju           Bakta, gapseq+SMETANA)
 optional: competitive DIAMOND gene-panel read search (e.g. mcrA vs mtmB/mtaB)
```

## 1. Install

You need conda/mamba and Snakemake ≥ 8:

```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake=8 pandas
mamba activate snakemake
```

Every step installs its own pinned tool environment on first run (`workflow/envs/`).

## 2. Databases (once)

```bash
bash workflow/scripts/download_databases.sh /path/to/databases
```

Then put the paths into `config/config.yaml → databases`. Required for the core
workflow: CheckM2, GUNC (proGenomes 2.1), GTDB-Tk (R232). Required for annotation:
KOfam, dbCAN, geNomad. Kraken2, Kaiju, Bakta and antiSMASH are optional — leave
the path empty and that step is skipped.

## 3. Describe your samples

**`config/samples.tsv`** — one row per sample

| column | meaning |
|---|---|
| `sample` | short name, letters/digits only (`SED`, `LIQ`, `Site3`) — used in contig and MAG names |
| `role` | `assemble` = assembled, binned and used as a coverage source for binning; `detect_only` = only recruited to the MAGs (e.g. public metagenomes, other sites) |

**`config/units.tsv`** — one row per sequencing run / FASTQ

| column | meaning |
|---|---|
| `sample` | which sample the run belongs to |
| `unit` | run name |
| `fastq` | path to `.fastq` or `.fastq.gz` |
| `profile` | `yes` = include in read-level profiling (SingleM/Kraken2/Kaiju/gene panel); `no` = use only for assembly/mapping (e.g. an amplicon run) |

All runs of a sample are concatenated and de-duplicated by read ID. Reads
≥ `min_len_assembly` (default 1 kb) go into Flye; **all** reads are used for
polishing, mapping and detection.

## 4. Run

```bash
snakemake --workflow-profile profiles/default -n          # dry run: shows what will be done
snakemake --workflow-profile profiles/default --cores 32   # real run
```

On a SLURM cluster add `--executor slurm` (plugin `snakemake-executor-plugin-slurm`).
The workflow is resumable: re-running only redoes missing or outdated steps.

## 5. Key settings (`config/config.yaml`)

| setting | default | note |
|---|---|---|
| `assembly.medaka_model` | `r1041_e82_400bps_sup_v5.0.0` | **must** match your flow cell and basecalling model |
| `assembly.flye_read_type` | `--nano-hq` | R10.4.1 / Dorado sup; use `--nano-raw` for older data |
| `binning.min_contig_len` / `min_bin_size` | 2000 / 200 000 | as used for Etoliko |
| `mag_filter` | ≥ 50 % complete, ≤ 10 % contaminated, GUNC pass | MIMAG medium quality and up |
| `detection.*` | identities 90/95/97/99, ≥ 50 % of read aligned | calibrate the thresholds on your own data (see below) |
| `modules.*` | — | switch optional parts on/off |

## 6. Outputs (`results/`)

| path | content |
|---|---|
| `qc/read_stats.tsv`, `qc/assembly_stats.tsv` | read and assembly statistics |
| `mapping/<asm>/flagstat/` | % of each read set mapping to each assembly |
| `mags/mag_summary.tsv` | every DAS_Tool bin: CheckM2, GUNC, rRNA, tRNA types, MIMAG, selected yes/no |
| `mags/MAG_CROSSWALK.tsv` | final MAG ID ↔ working bin name ↔ winning binner |
| `mags/fasta/` | the MAG catalogue |
| `taxonomy/gtdbtk_summary.tsv` | GTDB-Tk classification (bacteria + archaea, incl. RED) |
| `taxonomy/skani_*.tsv` | ANI among MAGs and against your reference/isolate genomes |
| `detection/coverm_all.tsv` | per MAG × read set: relative abundance, depth, breadth at each identity, read-identity mode, detection call |
| `read_profiling/` | SingleM profile + phylum abundance (filled coverage), Kraken2/Bracken, Kaiju |
| `annotation/summary/ko_matrix.tsv`, `cazy_matrix.tsv` | KO and CAZy family counts per MAG |
| `annotation/{kofam,dbcan,genomad,antismash,bakta}/` | full per-MAG results |
| `phylogeny/trees/` | IQ-TREE trees (`bac120`, `ar53`) |
| `gene_panel/` | competitive read counts per gene label |

## 7. Notes on interpretation (lessons from the Etoliko analysis)

* **Detection is not abundance.** A genome is called a *same-species population*
  only when the mode of its read-identity distribution is ≥ 97 % **and** breadth at
  ≥ 90 % identity is ≥ 0.5; *related population* when breadth at ≥ 90 % is ≥ 0.1 at
  lower identity. Check where reads from known-present genomes peak in your
  library (98–99 % for R10.4.1 sup) and adjust `detection.*`.
* **PCR-amplified libraries** distort abundance. Use them for assembly and
  detection, not for quantitative claims.
* **KofamScan curated thresholds only.** A flat E-value cut-off inflated hits by
  ~5× overall and ~29× for transporters.
* **Competitive read searches need decoys.** In the gene panel include all sibling
  families (e.g. MttB-superfamily methyltransferases alongside *mtmB/mtbB/mttB/mtaB*);
  without them most hits are mis-assigned.
* **skani reports only pairs above its screen** (~80 % ANI); a missing pair means
  "below the lowest value reported", not zero.
* Gene absence is only meaningful for genomes ≥ 90 % complete.

## 8. Not automated (case-specific steps of the original study)

Linking MAG 16S genes to amplicon OTUs (VSEARCH), species/genus delineation for
novel taxa (AAI with EzAAI, POCP, dDDH via the TYGS web server), HydDB hydrogenase
classification, marker-gene trees (e.g. McrA with MAFFT/trimAl/IQ-TREE), and
the choice of medium for SMETANA. The outputs above are the inputs for these.

## 9. Software versions used for the Etoliko analysis

Flye 2.9.3 · Medaka 2.1.1 · minimap2 2.28 · SAMtools 1.19.2 · MetaBAT2 2.18 ·
VAMB 5.0.4 · SemiBin2 2.2.1 · DAS_Tool 1.1.7 · CheckM2 1.1.0 · GUNC 1.0.3 ·
barrnap 0.9 · tRNAscan-SE 2.0.13 · GTDB-Tk 2.7.2 (GTDB R232) · CoverM 0.7.0 ·
IQ-TREE 3.1.1 · Kraken2 2.1.3 · Bracken 3.0.1 · Kaiju 1.10.1 · DIAMOND 2.1.9 ·
Bakta 1.9.4 · antiSMASH 8.0.4. Where bioconda lacks the exact build, the env file
says so.

## Citation

If you use this workflow, please cite the Etoliko metagenome paper and the tools above.
