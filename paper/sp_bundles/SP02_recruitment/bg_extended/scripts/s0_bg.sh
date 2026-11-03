#!/usr/bin/env bash
# bg_extended step 0 -- GTDB R232 species representatives of g__Clostridium_J, g__Terrisporobacter, g__Paraclostridium
# (taxonomy: ${SHARED_DB_DIR}/GTDB/release232/taxonomy/gtdb_taxonomy.tsv), downloaded from NCBI with `datasets`.
# Any representative at >=95% ANI (fastANI 1.34) to an isolate is EXCLUDED: it is the isolate's own species and would
# take the isolate's reads rather than absorbing relatives'. Output: bg/kept.tsv, bg/ani.tsv, bg/genomes/*.fna
set -euo pipefail
cd "$OUT/bg"
TAX=${SHARED_DB_DIR}/GTDB/release232/taxonomy/gtdb_taxonomy.tsv
FANI=$HOME/.conda/envs/coverm_env/bin/fastANI
MA=$B/UniqueELT_Solid/Manuscript_Analysis/Final_Genomes
grep -E "g__(Clostridium_J|Terrisporobacter|Paraclostridium);" "$TAX" | awk -F'\t' '{a=substr($1,4); s=$2; sub(/.*s__/,"",s); print a"\t"s}' > reps.tsv
echo "representatives: $(wc -l < reps.tsv)"
if [ ! -s genomes.done ]; then
  cut -f1 reps.tsv > acc.txt
  if ! unzip -tq reps.zip >/dev/null 2>&1; then
    for i in 1 2 3 4 5; do datasets download genome accession --inputfile acc.txt --include genome --filename reps.zip && break; sleep 60; done
  fi
  rm -rf ncbi_dataset genomes; unzip -q -o reps.zip; mkdir -p genomes
  while IFS=$'\t' read -r a s; do
    f=$(ls ncbi_dataset/data/GC?_${a#GC?_}/*.fna 2>/dev/null | head -1 || true)   # NCBI may file a GCF request under its GCA twin
    [ -n "$f" ] && cp "$f" genomes/$a.fna || echo "MISSING $a"
  done < reps.tsv
  ls genomes/*.fna | wc -l > genomes.done
fi
ls genomes/*.fna > bg_list.txt
: > ani.tsv
for g in YGCA01_Clean YGCA01_Contaminant YGCA10 YGCA07 YGCC01; do
  $FANI -q $MA/$g.fasta --rl bg_list.txt -o ani_$g.tmp -t 8 2>/dev/null
  awk -v g=$g '{print g"\t"$0}' ani_$g.tmp >> ani.tsv; rm -f ani_$g.tmp
done
$PY - <<'PY'
import os, csv
sp = dict(l.rstrip('\n').split('\t') for l in open('reps.tsv'))
best = {}
for l in open('ani.tsv'):
    g, q, r, ani = l.split('\t')[:4]; a = os.path.basename(r)[:-4]; ani = float(ani)
    if a not in best or ani > best[a][1]: best[a] = (g, ani)
with open('kept.tsv', 'w') as o:
    o.write('accession\tspecies\tfasta\tmax_ANI_to_isolate\tclosest_isolate\tdecision\n')
    for a, s in sp.items():
        p = f'{os.getcwd()}/genomes/{a}.fna'
        if not os.path.exists(p): o.write(f'{a}\t{s}\t-\t-\t-\tnot downloaded\n'); continue
        g, ani = best.get(a, ('-', 0.0))
        d = 'EXCLUDED (same species as isolate, ANI>=95)' if ani >= 95 else 'kept'
        o.write(f'{a}\t{s}\t{p}\t{ani if ani else "<80 (no fastANI output)"}\t{g}\t{d}\n')
n = sum(1 for l in open('kept.tsv') if l.rstrip().endswith('\tkept'))
print('kept', n, 'of', len(sp)); print(open('kept.tsv').read())
PY
