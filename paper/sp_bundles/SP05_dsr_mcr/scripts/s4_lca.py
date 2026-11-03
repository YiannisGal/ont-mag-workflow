"""SP05 step 2 -- LCA of the top-10 DIAMOND hits (by bitscore, SwissProt + TrEMBL pooled) for every contig protein; NCBI taxdump
lineages from the OX= taxid in each UniProt title. Also the KOfam call per protein (K06313 spore-germination check on 19299)."""
import os, re, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__)); import kofam
OUT = os.environ['OUT']
par, rank, name = {}, {}, {}
for l in open(f'{OUT}/lca/taxdump/nodes.dmp'):
    f = l.split('\t|\t'); par[f[0]] = f[1]; rank[f[0]] = f[2]
for l in open(f'{OUT}/lca/taxdump/names.dmp'):
    f = l.split('\t|\t')
    if f[3].startswith('scientific name'): name[f[0]] = f[1]
def lin(t):
    out = []
    while t in par and t != '1': out.append(t); t = par[t]
    return out[::-1]
hits = defaultdict(list)
for db in ('uniprot_sprot', 'uniprot_trembl'):
    p = f'{OUT}/lca/hits_{db}.tsv'
    if not os.path.exists(p): continue
    for l in open(p):
        f = l.rstrip('\n').split('\t'); m = re.search(r'OX=(\d+)', f[6])
        if m: hits[f[0]].append((float(f[5]), float(f[2]), m.group(1), f[6].split(' OS=')[0].split(' ', 1)[-1]))
H = kofam.load()
prots = [l[1:].strip() for l in open(f'{OUT}/work/contig_proteins.faa') if l.startswith('>')]
RK = ['superkingdom', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
with open(f'{OUT}/contig_protein_LCA.tsv', 'w') as o:
    o.write('protein\tKO\tKO_definition\tn_hits\ttop_hit_pid\ttop_hit_product\tLCA_rank\tLCA_name\tLCA_lineage\ttop_hit_organism\n')
    for p in prots:
        top = sorted(hits.get(p, []), key=lambda x: -x[0])[:10]
        bp = kofam.best_pass(H, p)
        if not top:
            o.write(f"{p}\t{bp[0] if bp else '-'}\t{bp[3] if bp else '-'}\t0\t-\t-\t-\tno hit\t-\t-\n"); continue
        L = [lin(t) for _, _, t, _ in top]
        common = [x for x in L[0] if all(x in l for l in L[1:])]
        lca = common[-1] if common else '1'
        lg = ';'.join(f"{rank[x]}:{name.get(x,'?')}" for x in common if rank.get(x) in RK)
        top_org = name.get(top[0][2], top[0][2])
        o.write(f"{p}\t{bp[0] if bp else '-'}\t{bp[3] if bp else '-'}\t{len(top)}\t{top[0][1]:.1f}\t{top[0][3]}\t{rank.get(lca,'-')}\t{name.get(lca,'root')}\t{lg}\t{top_org}\n")
print(open(f'{OUT}/contig_protein_LCA.tsv').read())
k6 = [p for p in prots if p.startswith('LIQ_contig_19299') and (kofam.best_pass(H, p) or ('',))[0] == 'K06313']
print('K06313 (spore germination) on LIQ_contig_19299:', k6 or 'none above threshold')
open(f'{OUT}/work/K06313_check.txt', 'w').write(f'K06313 on LIQ_contig_19299: {k6 or "none above threshold"}\n')
