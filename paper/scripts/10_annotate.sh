#!/usr/bin/env bash
# =============================================================================
# 10_annotate.sh  --  Bakta annotation of all 29 MAGs (bacteria AND archaea)
#
# PURPOSE
#   Produce the gene calls and functional assignments that every downstream
#   analysis reads. This runs on ALL 29 MAGs — the 21 bacteria are not a
#   background to the archaea, they are the other half of the syntrophy
#   argument: bacterial fermentation supplies the H2, formate and acetate that
#   the eight methanogens consume. Both halves need equal annotation depth.
#
# WHY BAKTA WITH THE FULL DATABASE
#   The full DB (36.7 GB, schema v5.1, 2024-01-19) carries the UniRef/PSC
#   protein-sequence-cluster layer that the light DB (1.3 GB) lacks. Without it
#   most CDS come back as "hypothetical protein" and the metabolic
#   reconstruction at step 11 has nothing to work with. This is why the full DB
#   was downloaded at step 02.
#
# TAXONOMY-AWARE
#   Genus and species from step 09 are passed to Bakta per genome. This affects
#   gene naming and the reference protein set used, so annotations are
#   comparable to published genomes of the same lineage rather than generic.
#
# INDEPENDENT CHECK ON THE GENE-LENGTH QUESTION
#   Bakta calls genes with Prodigal/Pyrodigal independently of CheckM2. The
#   summary below reports mean CDS length per MAG, which can be compared with
#   CheckM2's Average_Gene_Length. Two independent gene callers agreeing makes
#   the 08c conclusion (281 aa is biology, not error) much harder to argue with.
#
# INPUT    07_mags/bins/*.fa   +   07_mags/MAG_TAXONOMY.tsv
# OUTPUT   09_annotation/<working_name>/   (gff3, gbff, faa, ffn, tsv, json)
#          09_annotation/annotation_summary.tsv
#
# RUNTIME  2-4 h for 29 genomes. RUN IN tmux.
# USAGE    ws2 && bash scripts/10_annotate.sh 2>&1 | tee logs/10_annotate_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="10_annotate"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
F="$WS2/07_mags"; A="$WS2/09_annotation"; mkdir -p "$A"

echo "######## 10 — BAKTA ANNOTATION ########"
echo "  started: $(date -Is)  threads: $THREADS"
runlog "10_annotate started (Bakta full DB)"

[ -d "${BAKTA_DB:-/nonexistent}" ] || { echo "  !! BAKTA_DB not a directory: ${BAKTA_DB:-UNSET}"; exit 1; }
if [ -s "$BAKTA_DB/version.json" ]; then
  echo "  [db] $BAKTA_DB"
  grep -oE '"(type|major|minor)": *"?[A-Za-z0-9]+"?' "$BAKTA_DB/version.json" | tr '\n' ' '; echo
fi
[ -s "$BAKTA_DB/psc.dmnd" ] && echo "  [db] psc.dmnd present — full DB confirmed" \
  || echo "  !! psc.dmnd MISSING — this looks like the LIGHT db; annotations will be poor"

record_version "bakta" "${ENV_bakta:-NB_Genome}" bakta --version

TAX="$F/MAG_TAXONOMY.tsv"
[ -s "$TAX" ] || { echo "  !! missing $TAX — run 09_taxonomy.sh first"; exit 1; }
NMAG=$(ls -1 "$F/bins"/*.fa 2>/dev/null | wc -l)
echo "  genomes to annotate: $NMAG"

# --------------------------------------------------- per-genome taxonomy lookup
taxfield() {  # taxfield <working_name> <column_name>
  awk -F'\t' -v n="$1" -v c="$2" '
    NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
    $h["working_name"]==n{print $h[c]}' "$TAX"
}

# ---------------------------------------------------------------- ANNOTATE
i=0
for f in "$F/bins"/*.fa; do
  name=$(basename "$f" .fa); i=$((i+1))
  out="$A/$name"
  if [ -s "$out/${name}.tsv" ]; then
    printf '  [%2d/%2d] %-28s already done\n' "$i" "$NMAG" "$name"; continue
  fi
  genus=$(taxfield "$name" genus)
  sp=$(taxfield "$name" species)
  magid=$(taxfield "$name" mag_id)
  # GTDB species strings look like "Genus species" — keep only the epithet
  epithet=$(echo "$sp" | awk '{print $2}')

  # Bakta enforces the INSDC locus-tag format: ALPHANUMERIC ONLY, no
  # underscores. ETO_LIQ_MAG_001 is rejected; ETOLIQMAG001 is accepted.
  ltag=$(echo "${magid:-$name}" | tr -cd '[:alnum:]' | tr '[:lower:]' '[:upper:]' | cut -c1-16)
  [ -z "$ltag" ] && ltag="MAG$(printf '%03d' "$i")"

  # Pass taxonomy only when GTDB actually assigned it. Placeholder GTDB genus
  # names (JBAIYS01, UBA6094, JAAYQI01 ...) are alphanumeric codes, not Latin
  # binomials, and are safe to pass but meaningless — keep them for traceability.
  taxargs=()
  [ -n "$genus" ]   && taxargs+=(--genus "$genus")
  [ -n "$epithet" ] && taxargs+=(--species "$epithet")
  [ -n "$magid" ]   && taxargs+=(--strain "$magid")

  printf '  [%2d/%2d] %-28s %-18s %-12s tag=%s\n' \
    "$i" "$NMAG" "$name" "${genus:-<no genus>}" "${epithet:-sp.}" "$ltag"
  mkdir -p "$out"
  run_tool bakta \
      --db "$BAKTA_DB" \
      --output "$out" \
      --prefix "$name" \
      "${taxargs[@]}" \
      --locus-tag "$ltag" \
      --threads "$THREADS" \
      --skip-plot \
      --force \
      "$f" > "$out/bakta.log" 2>&1 \
    || { echo "        !! bakta failed — last lines of the log:"
         tail -5 "$out/bakta.log" | sed 's/^/           /'; }
done

# ------------------------------------------------------------------ SUMMARY
echo
echo "######## ANNOTATION SUMMARY ########"
python3 - "$A" "$F" <<'PY'
import sys, os, csv, glob, statistics as st

A, F = sys.argv[1], sys.argv[2]
tax = {r['working_name']: r for r in csv.DictReader(open(f"{F}/MAG_TAXONOMY.tsv"), delimiter='\t')}

rows = []
for d in sorted(glob.glob(os.path.join(A, '*'))):
    if not os.path.isdir(d): continue
    name = os.path.basename(d)
    tsv = os.path.join(d, f'{name}.tsv')
    faa = os.path.join(d, f'{name}.faa')
    if not os.path.exists(tsv): continue

    cds = hyp = trna = rrna = 0
    with open(tsv) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 8: continue
            t = p[1]
            if t == 'cds':
                cds += 1
                # Bakta TSV: seqid,type,start,stop,strand,locus_tag,gene,PRODUCT,dbxrefs
                # The product is column 8 (index 7). p[-1] is dbxrefs and never
                # contains 'hypothetical' — using it reported 0.0% everywhere.
                product = p[7].lower() if len(p) > 7 else ''
                if 'hypothetical' in product or 'uncharacterized' in product:
                    hyp += 1
            elif t == 'tRNA':  trna += 1
            elif t == 'rRNA':  rrna += 1

    lens = []
    if os.path.exists(faa):
        cur = 0
        for line in open(faa):
            if line.startswith('>'):
                if cur: lens.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur: lens.append(cur)

    t = tax.get(name, {})
    rows.append(dict(
        mag_id=t.get('mag_id', ''), working_name=name, sample=t.get('sample', ''),
        domain=t.get('domain', ''), phylum=t.get('phylum', ''), genus=t.get('genus', ''),
        n_cds=cds, n_hypothetical=hyp,
        pct_hypothetical=f"{100*hyp/cds:.1f}" if cds else '',
        n_trna=trna, n_rrna=rrna,
        mean_cds_len_aa=f"{st.mean(lens):.1f}" if lens else '',
        median_cds_len_aa=f"{st.median(lens):.1f}" if lens else '',
        checkm2_avg_gene_length=t.get('avg_gene_length', '')))

if not rows:
    print("  !! no Bakta output parsed"); raise SystemExit(1)

with open(os.path.join(A, 'annotation_summary.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(rows)

print(f"  {'MAG ID':<17} {'dom':<4} {'genus':<18} {'CDS':>6} {'hyp%':>6} "
      f"{'meanAA':>7} {'CheckM2':>8} {'tRNA':>5} {'rRNA':>5}")
for r in sorted(rows, key=lambda r: (r['sample'] != 'SED', r['mag_id'])):
    print(f"  {r['mag_id']:<17} {(r['domain'] or '?')[:3]:<4} {(r['genus'] or '-')[:18]:<18} "
          f"{r['n_cds']:>6} {r['pct_hypothetical']:>6} {r['mean_cds_len_aa']:>7} "
          f"{(r['checkm2_avg_gene_length'] or '')[:8]:>8} {r['n_trna']:>5} {r['n_rrna']:>5}")

# ------------------------------- independent gene-length cross-check
pairs = [(float(r['mean_cds_len_aa']), float(r['checkm2_avg_gene_length']))
         for r in rows if r['mean_cds_len_aa'] and r['checkm2_avg_gene_length']]
print()
print("  GENE LENGTH: BAKTA vs CheckM2 — independent gene callers")
if pairs:
    b = [x for x, _ in pairs]; c = [y for _, y in pairs]
    d = [x - y for x, y in pairs]
    print(f"    Bakta  median mean-CDS : {st.median(b):.1f} aa")
    print(f"    CheckM2 median AGL     : {st.median(c):.1f} aa")
    print(f"    median difference      : {st.median(d):+.1f} aa")
    print(f"    workspace_1 baseline   : 244 aa")
    print()
    if abs(st.median(d)) < 25:
        print("    >> The two callers agree. The gene-length values are a property of")
        print("       the genomes, not an artefact of one tool's calling. This")
        print("       independently supports the audit correction in §0b.")
    else:
        print("    >> The callers DISAGREE by more than 25 aa. Investigate before")
        print("       using either number: check translation tables and whether")
        print("       Bakta and CheckM2 filter short ORFs differently.")

# hypothetical fraction is the test of whether the full DB is really in use
hyps = [float(r['pct_hypothetical']) for r in rows if r['pct_hypothetical']]
if hyps:
    print()
    print(f"  HYPOTHETICAL PROTEINS: median {st.median(hyps):.1f}% "
          f"(range {min(hyps):.1f}-{max(hyps):.1f}%)")
    if st.median(hyps) > 45:
        print("    !! High. Expected ~25-40% for environmental genomes with the FULL")
        print("       Bakta DB. If this is near 60%+, check that psc.dmnd was used —")
        print("       it would mean the light DB was picked up instead.")
    else:
        print("    Within the expected range for environmental genomes on the full DB.")
PY

# ------------------------------------------------------------- VALIDATION
echo
echo "######## VALIDATION GATE ########"
done_n=$(ls -1d "$A"/*/ 2>/dev/null | while read -r d; do
           n=$(basename "$d"); [ -s "$d/$n.tsv" ] && echo x; done | wc -l)
check "every MAG annotated ($done_n/$NMAG)" $([ "$done_n" -eq "$NMAG" ] && echo 0 || echo 1)
check "annotation summary written" $([ -s "$A/annotation_summary.tsv" ] && echo 0 || echo 1)
# Bakta writes BOTH <name>.faa and <name>.hypotheticals.faa — count only the
# main protein file, or this reports 2x NMAG and fails spuriously.
faa_n=$(ls -1d "$A"/*/ 2>/dev/null | while read -r d; do
          n=$(basename "$d"); [ -s "$d/$n.faa" ] && echo x; done | wc -l)
check "protein FASTA per MAG ($faa_n)" $([ "$faa_n" -eq "$NMAG" ] && echo 0 || echo 1)
gff_n=$(ls -1 "$A"/*/*.gff3 2>/dev/null | wc -l)
check "GFF3 per MAG ($gff_n)" $([ "$gff_n" -eq "$NMAG" ] && echo 0 || echo 1)
empty=$(find "$A" -name '*.faa' -size -1k 2>/dev/null | wc -l)
check "no empty protein files ($empty)" $([ "$empty" -eq 0 ] && echo 0 || echo 1)

runlog "10_annotate complete ($done_n/$NMAG)"
echo
echo "######## 10 DONE — next: 11_metabolism.sh ########"
