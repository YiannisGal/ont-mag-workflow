"""CAT-like contig classification: per protein LCA of top-10 DIAMOND hits (UniProt TrEMBL, NCBI taxdump), contig phylum = majority of
classified proteins (>50%). Targets: NCBI Chloroflexota and Thermodesulfobacteriota (= GTDB Desulfobacterota)."""
import sys, re, collections
TD = sys.argv[1]; PROT = sys.argv[2] if len(sys.argv) > 2 else 'sed.faa'; HITS = sys.argv[3] if len(sys.argv) > 3 else 'sed_vs_trembl.tsv'; par, rank, name = {}, {}, {}
for l in open(f'{TD}/nodes.dmp'):
    f = l.split('\t|\t'); par[f[0]] = f[1]; rank[f[0]] = f[2]
for l in open(f'{TD}/names.dmp'):
    f = l.split('\t|\t')
    if f[3].startswith('scientific name'): name[f[0]] = f[1]
def lin(t):
    out = []
    while t in par and t != '1': out.append(t); t = par[t]
    return out[::-1]
hits = collections.defaultdict(list)
for l in open(HITS):
    f = l.rstrip('\n').split('\t'); m = re.search(r'OX=(\d+)', f[6])
    if m and len(hits[f[0]]) < 10: hits[f[0]].append(m.group(1))
phy = collections.defaultdict(collections.Counter); nprot = collections.Counter()
for l in open(PROT):
    if l.startswith('>'): nprot[l[1:].split()[0].rsplit('_', 1)[0]] += 1
for p, ts in hits.items():
    L = [lin(t) for t in ts]; common = [x for x in L[0] if all(x in l for l in L[1:])]
    ph = next((name[x] for x in common if rank.get(x) == 'phylum'), None)
    if ph: phy[p.rsplit('_', 1)[0]][ph] += 1
T = {'Chloroflexota', 'Thermodesulfobacteriota', 'Desulfobacterota'}
with open('contig_phylum.tsv', 'w') as o, open('target_contigs.tsv', 'w') as t:
    o.write('contig\tn_proteins\tn_phylum_classified\tphylum\tsupport\n'); t.write('contig\tphylum\tsupport\n')
    for c, n in nprot.items():
        cnt = phy.get(c, collections.Counter()); tot = sum(cnt.values())
        top, k = (cnt.most_common(1)[0] if cnt else ('none', 0))
        call = top if tot and k / tot > 0.5 else 'unclassified'
        o.write(f'{c}\t{n}\t{tot}\t{call}\t{(k / tot if tot else 0):.2f}\n')
        if call in T: t.write(f'{c}\t{call}\t{k / tot:.2f}\n')
print(open('target_contigs.tsv').read().count('\n') - 1, 'target contigs')
