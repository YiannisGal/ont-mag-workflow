"""LCA of the top-10 DIAMOND blastx hits (bitscore, SwissProt + TrEMBL pooled) per read; NCBI taxdump lineage from OX= (as SP05 s4_lca.py)."""
import re, sys, os
from collections import defaultdict
U = os.path.dirname(os.path.abspath(__file__)); TD = sys.argv[1]
par, rank, name = {}, {}, {}
for l in open(f'{TD}/nodes.dmp'):
    f = l.split('\t|\t'); par[f[0]] = f[1]; rank[f[0]] = f[2]
for l in open(f'{TD}/names.dmp'):
    f = l.split('\t|\t')
    if f[3].startswith('scientific name'): name[f[0]] = f[1]
def lin(t):
    out = []
    while t in par and t != '1': out.append(t); t = par[t]
    return out[::-1]
hits = defaultdict(list)
for db in ('sprot', 'trembl'):
    for l in open(f'{U}/hits_{db}.tsv'):
        f = l.rstrip('\n').split('\t'); m = re.search(r'OX=(\d+)', f[6])
        if m: hits[f[0]].append((float(f[5]), float(f[2]), m.group(1), f[6].split(' OS=')[0]))
RK = ['domain', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
ids = [l.split('\t')[1] for l in open(f'{U}/../reads/read_ids.tsv') if not l.startswith('gene\t')]
with open(f'{U}/read_LCA.tsv', 'w') as o:
    o.write('read_id\tn_hits\tuniprot_top_hit\tuniprot_top_pident\tuniprot_LCA_rank\tuniprot_LCA\tuniprot_LCA_lineage\n')
    for q in ids:
        top = sorted(hits.get(q, []), key=lambda x: -x[0])[:10]
        if not top: o.write(f'{q}\t0\tnone\tnone\tnone\tno hit\tnone\n'); continue
        L = [lin(t) for *_, t, _ in [(a, b, c, d) for a, b, c, d in top]]
        common = [x for x in L[0] if all(x in l for l in L[1:])]; lca = common[-1] if common else '1'
        lg = ';'.join(f"{rank[x]}:{name.get(x,'?')}" for x in common if rank.get(x) in RK)
        o.write(f"{q}\t{len(top)}\t{top[0][3]} [{name.get(top[0][2], top[0][2])}]\t{top[0][1]:.1f}\t{rank.get(lca,'-')}\t{name.get(lca,'root')}\t{lg}\n")
print('LCA written')
