#!/usr/bin/env bash
# SP01 step 4a -- NifH (K02588) and NifD (K02586) hits above threshold placed in a reference tree
# (mafft --auto, FastTree -lg), each classified by its nearest labelled reference (patristic distance).
# Reference labels:
#   NIF_I   reviewed NifH/NifD from Pseudomonadota(alpha/beta/gamma), Cyanobacteriota, Actinomycetota
#   NIF_III reviewed NifH/NifD from all other lineages (anaerobes, methanogen diazotrophs)
#   ANF / VNF  alternative nitrogenases (gene:anfH/anfD, vnfH/vnfD)
#   CFB     CfbC/CfbD, coenzyme F430 biosynthesis (gene:cfbC/cfbD)
#   NFL     group IV Nfl: M. jannaschii MJ0879/MJ0880 (non-diazotroph; Staples et al. 2007 NflH/NflD)
#           and R. capsulatus IsrH/IsrD (nitrogenase-like isethionate reductase)
#   BCH     group V outgroup: BchL/ChlL/BchX (H-like) and BchN/ChlN/BchY (D-like), reviewed
set -euo pipefail
N=$OUT/nif; mkdir -p "$N"; cd "$N"
F="accession,gene_names,organism_name,lineage,sequence"
up(){ curl -sf -m 180 "https://rest.uniprot.org/uniprotkb/search?format=tsv&fields=$F&size=500&query=$2" | tail -n +2 | sed "s/^/$1\t/"; }
if [ ! -s refs.tsv ]; then
  { up H_NIF '(gene:nifH)+AND+(reviewed:true)'; up D_NIF '(gene:nifD)+AND+(reviewed:true)'
    up H_ANF '(gene:anfH)'; up D_ANF '(gene:anfD)'; up H_VNF '(gene:vnfH)'; up D_VNF '(gene:vnfD)'
    up H_CFB '(gene:cfbC)'; up D_CFB '(gene:cfbD)'
    up H_NFL '(accession:Q58289)+OR+(accession:D5AKX6)'
    up D_NFL '(accession:Q58290)+OR+(accession:D5AKX5)'
    up H_BCH '((gene:bchL)+OR+(gene:chlL)+OR+(gene:bchX))+AND+(reviewed:true)'
    up D_BCH '((gene:bchN)+OR+(gene:chlN)+OR+(gene:bchY))+AND+(reviewed:true)'
  } > refs.tsv.part
  mv refs.tsv.part refs.tsv
fi
$PY - <<'PY'
import os, csv
from collections import Counter
OUT = os.environ['OUT']; N = f'{OUT}/nif'
G1 = ('Alphaproteobacteria', 'Betaproteobacteria', 'Gammaproteobacteria', 'Cyanobacteriota', 'Cyanobacteria', 'Actinomycetota', 'Actinobacteria')
refs = {'H': [], 'D': []}; seen = set()
for l in open(f'{N}/refs.tsv'):
    f = l.rstrip('\n').split('\t')
    if len(f) < 6: continue
    tag, acc, genes, org, lin, seq = f[:6]
    sub, lab = tag.split('_')
    if lab == 'NIF':
        if acc in ('Q58289', 'Q58290') or 'Methanocaldococcus' in org:   # non-diazotroph group IV annotated "nifH"
            lab = 'NFL'
        else:
            lab = 'NIF_I' if any(g in lin for g in G1) else 'NIF_III'
    if (sub, acc) in seen: continue
    seen.add((sub, acc)); refs[sub].append((f'REF_{lab}__{acc}', seq))
want = {}
for r in csv.DictReader(open(f'{OUT}/markers/calls_detail.tsv'), delimiter='\t'):
    if r['KO'] in ('K02588', 'K02586') and r['call'] in ('call', 'fragment'):
        want[r['prot']] = 'H' if r['KO'] == 'K02588' else 'D'
seqs = {}; name = None
for l in open(f'{OUT}/markers/all34.faa'):
    if l.startswith('>'): name = l[1:].strip(); seqs[name] = [] if name in want else None
    elif seqs.get(name) is not None: seqs[name].append(l.strip())
for sub in ('H', 'D'):
    q = [p for p, s in want.items() if s == sub]
    with open(f'{N}/nif{sub}_in.faa', 'w') as o:
        for acc, s in refs[sub]: o.write(f'>{acc}\n{s}\n')
        for p in q: o.write(f">Q__{p.replace('|', '__')}\n{''.join(seqs[p])}\n")
    print(sub, 'refs', dict(Counter(a.split('__')[0] for a, _ in refs[sub])), 'queries', len(q))
PY
for s in H D; do
  $PIXI/mafft --auto --thread "$CPU" "nif${s}_in.faa" > "nif${s}.aln" 2> "mafft_${s}.err"
  $PIXI/FastTree -lg "nif${s}.aln" > "nif${s}.tree" 2> "fasttree_${s}.err"
done
$PY "$OUT/scripts/s4_nif_classify.py"
