# How Combined_analysis called its marker matrix (and what "31" / "34" / "36" are)

Source tree (read-only): `${COMBINED_WS}`

## What the numbers count

- `00_inputs/metagenomes/mags/` holds **29 MAG FASTAs**. Classification is in `provenance_31.tsv`.
- `00_inputs/metagenomes/mags/MAG_QC.tsv`: 47 data rows x 19 data columns
- `00_inputs/metagenomes/mags/MAG_TAXONOMY.tsv`: 29 data rows x 38 data columns
- `05_all_31_mags_analysis/all_36_genomes_metabolic_matrix.tsv`: 23 data rows x 34 data columns
- `05_all_31_mags_analysis/integrated_taxonomy_metabolism.tsv`: 48 data rows x 39 data columns
- `05_all_31_mags_analysis/mag_summary_curated.tsv`: 0 data rows x 59 data columns
- `05_all_31_mags_analysis/master_mag_catalog_annotated.tsv`: 48 data rows x 16 data columns

"31" appears only in the step name (`06_expand_to_all_31_mags.sh`, `05_all_31_mags_analysis/`), in RUN_05 and in Chapter 6 prose. No file, table or loop in Combined_analysis enumerates 31 genomes: the matrix is built from `glob(mags/*.fa)` (29) plus 5 hard-coded isolates = 34 columns, written to a file named `all_36_genomes_...`. The MAG_QC.tsv/master catalog rows are the workspace_2 union-bin QC table, not a MAG set.

## Source of the species-style MAG names

`scripts/02_run_read_recruitment.sh` lines 60-63 (verbatim):

```bash
append_genome "$DIR_INPUT/metagenomes/mags/LIQ_das_vamb_002.fa" "MAG_Lysinibacillus_sphaericus" "MAG"
append_genome "$DIR_INPUT/metagenomes/mags/LIQ_das_semibin_003.fa" "MAG_Syntrophomonas_sp" "MAG"
append_genome "$DIR_INPUT/metagenomes/mags/LIQ_das_vamb_017.fa" "MAG_Methanofollis_formosanus" "MAG"
append_genome "$DIR_INPUT/metagenomes/mags/LIQ_das_semibin_016.fa" "MAG_Methanofollis_aquaemaris" "MAG"
```

The same mapping is re-typed in `04_build_metabolic_matrix.sh` (`mag_file_map`). These are hand labels on four workspace_2 bins; they are not GTDB-Tk R232 names from 08d.

## Marker calling in 06_expand_to_all_31_mags.sh (the 34-genome matrix)

Verbatim, lines 35-94:

```python
MARKERS = {
    # Methanogenesis
    "mcrA": ["K00399"], "mcrB": ["K00401"], "mcrG": ["K00402"],
    "mtrA": ["K00577"], "mer":  ["K00320"], "mtd":  ["K00319"],
    # Sulfur
    "dsrA": ["K11180"], "dsrB": ["K11181"], "aprA": ["K00394"], "aprB": ["K00395"], "sat":  ["K00958"],
    # Hydrogen & Formate (IHT)
    "fdhA": ["K00123"], "hydA": ["K00436"], "echA": ["K14083"],
    # Nitrogen & Ureolysis
    "ureA": ["K01428"], "ureB": ["K01429"], "ureC": ["K01430"],
    "nifH": ["K02588"], "nifD": ["K02586"], "nifK": ["K02591"],
    # Beta-oxidation
    "fadA": ["K00632"], "fadB": ["K01782"], "fadD": ["K01897"]
}

ko_to_marker = {ko: m for m, kos in MARKERS.items() for ko in kos}

mag_files = sorted(glob.glob(os.path.join(mag_input_dir, "*.fa")))
mag_ids = [os.path.splitext(os.path.basename(f))[0] for f in mag_files]

isolate_ids = [
    "YGCC01_Sporosarcina", "YGCA07_P_benzoelyticum", "YGCA10_P_benzoelyticum",
    "YGCA01_C_sulfidigenes", "YGCA01_T_petrolearius"
]

all_targets = isolate_ids + mag_ids
matrix = {g: {m: 0 for m in MARKERS} for g in all_targets}

# Parse Kofam for MAGs
kofam_hits_dir = os.path.join(src_meta, "11_kofam", "hits")
for mag in mag_ids:
    tbl_path = os.path.join(kofam_hits_dir, f"{mag}.tbl")
    if os.path.exists(tbl_path):
        with open(tbl_path) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split()
                if len(parts) >= 3 and parts[2] in ko_to_marker:
                    matrix[mag][ko_to_marker[parts[2]]] += 1

# Parse Culturomics GFFs
gff_map = {
    "YGCC01_Sporosarcina":    "YGCC01.gff",
    "YGCA07_P_benzoelyticum": "YGCA07.gff",
    "YGCA10_P_benzoelyticum": "YGCA10.gff",
    "YGCA01_C_sulfidigenes":  "YGCA01_C.gff",
    "YGCA01_T_petrolearius":  "YGCA01_Contam.gff"
}

for iso_id, gff_name in gff_map.items():
    gff_path = os.path.join(cult_annot_dir, gff_name)
    if os.path.exists(gff_path):
        with open(gff_path) as f:
            for line in f:
                if line.startswith("#"): continue
                for m in MARKERS:
                    if f"gene={m}" in line or f"product={m}" in line or f" {m} " in line:
                        matrix[iso_id][m] += 1

df = pd.DataFrame(matrix)
```

04_build_metabolic_matrix.sh (the 9-genome matrix) uses the same MARKERS dict, lines 34-129.

## The hit tables it read

MAG hits came from `workspace_2/11_kofam/hits/<bin>.tbl`, produced by `scripts/11b_kofam_archaea.sh` line 274:

```bash
  run_tool hmmsearch --cpu "$THREADS" --tblout "$out" -E 1e-5 \
      "$K/kofam_targets.hmm" "$faa" > /dev/null 2>&1 \
```

## Method, stated plainly

1. **MAGs:** every row of an `hmmsearch -E 1e-5 --tblout` table whose KO is in MARKERS adds 1 to the cell. **No KOfam curated threshold** and no score-type or length filter. Cells are raw row counts, so any homolog above E=1e-5 counts.
2. **Isolates:** substring match of `gene=<marker>`, `product=<marker>` or ` <marker> ` against each Prokka GFF line. This is name matching, not KO assignment. It is a different method from the MAGs, so the columns are not comparable.
3. **KO labels are wrong in several rows** (checked against `11_kofam/ko_list`):

| label used | KO used | ko_list definition of that KO |
|---|---|---|
| mcrA | K00399 | methyl-coenzyme M reductase alpha subunit [EC:2.8.4.1] |
| mcrB | K00401 | methyl-coenzyme M reductase beta subunit [EC:2.8.4.1] |
| mcrG | K00402 | methyl-coenzyme M reductase gamma subunit [EC:2.8.4.1] |
| mtrA | K00577 | tetrahydromethanopterin S-methyltransferase subunit A [EC:7.2.1.4] |
| mer | K00320 | 5,10-methylenetetrahydromethanopterin reductase [EC:1.5.98.2] |
| mtd | K00319 | methylenetetrahydromethanopterin dehydrogenase [EC:1.5.98.1] |
| dsrA | K11180 | dissimilatory sulfite reductase alpha subunit [EC:1.8.1.22] |
| dsrB | K11181 | dissimilatory sulfite reductase beta subunit [EC:1.8.1.22] |
| aprA | K00394 | adenylylsulfate reductase, subunit A [EC:1.8.99.2] |
| aprB | K00395 | adenylylsulfate reductase, subunit B [EC:1.8.99.2] |
| sat | K00958 | sulfate adenylyltransferase [EC:2.7.7.4] |
| fdhA | K00123 | formate dehydrogenase major subunit [EC:1.17.1.9] |
| hydA | K00436 | NAD-reducing hydrogenase large subunit [EC:1.12.1.2] |
| echA | K14083 | trimethylamine---corrinoid protein Co-methyltransferase [EC:2.1.1.250] |
| ureA | K01428 | urease subunit alpha [EC:3.5.1.5] |
| ureB | K01429 | urease subunit beta [EC:3.5.1.5] |
| ureC | K01430 | urease subunit gamma [EC:3.5.1.5] |
| nifH | K02588 | nitrogenase iron protein NifH |
| nifD | K02586 | nitrogenase molybdenum-iron protein alpha chain [EC:1.18.6.1] |
| nifK | K02591 | nitrogenase molybdenum-iron protein beta chain [EC:1.18.6.1] |
| fadA | K00632 | acetyl-CoA acyltransferase [EC:2.3.1.16] |
| fadB | K01782 | 3-hydroxyacyl-CoA dehydrogenase / enoyl-CoA hydratase / 3-hydroxybutyryl-CoA epimerase [EC:1.1.1.35 4.2.1.17 5.1.2.3] |
| fadD | K01897 | long-chain-fatty-acid---CoA ligase [EC:6.2.1.3] |

Label errors visible above: ureA/ureC are swapped (K01428 is ureC, K01430 is ureA); `hydA`=K00436 is the NiFe hoxH, not [FeFe] HydA; `fadB`=K01782 is the multifunctional fadJ/fadB of E. coli-type bacteria; `echA`=K14083 is the trimethylamine methyltransferase MttB, not the Ech hydrogenase (K14086); `fdhA`=K00123 is the Fdh alpha-subunit, and without a threshold it also collects NADH:quinone oxidoreductase and other Complex-I/molybdopterin paralogs.

Conclusion: the Combined_analysis matrix is presence of any E<1e-5 homolog for MAGs and a name search for isolates. It is not a curated-threshold KO call and must not be cited as one. SP01 replaces it (`markers_long.tsv`, `markers_wide.tsv`).
