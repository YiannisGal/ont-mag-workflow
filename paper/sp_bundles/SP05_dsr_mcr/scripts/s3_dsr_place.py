"""SP05 step 3 -- DsrA / DsrB / DsrAB placement (IQ-TREE LG+C20+F+G, UFBoot 1000).
Per query: 5 nearest references (patristic), nearest genus, oxidative-vs-reductive call, and the UFBoot support of the smallest
clade joining the query to its nearest reference. CONTROL: some internal edge must separate exactly the oxidative references
(Allochromatium, Chlorobium) from every reductive reference, in each tree."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__)); from tree import T
OUT = os.environ['OUT']; rows = []; qc = []; ok_all = True
for m in ('dsrA', 'dsrB', 'dsrAB'):
    p = f'{OUT}/trees/{m}.treefile'
    if not os.path.exists(p): qc.append(f'[FAIL] {m}: no tree'); ok_all = False; continue
    t = T(open(p).read()); L = t.leaves()
    refs = [n for n in L if n.startswith('REF_')]; oxi = {n for n in refs if n.startswith('REF_OXI')}; red = set(refs) - oxi
    good = any((s & set(refs)) == oxi or (set(refs) - s) == oxi for s, sup in t.splits())
    sup = [sup for s, sup in t.splits() if (s & set(refs)) == oxi or (set(refs) - s) == oxi]
    qc.append(f"[{'PASS' if good else 'FAIL'}] {m}: oxidative refs ({len(oxi)}: Allochromatium/Chlorobium) separated from {len(red)} reductive refs by one edge"
              f"{' (UFBoot ' + str(sup[0]) + ')' if sup else ''}")
    ok_all &= good
    for q in [n for n in L if n.startswith('QUERY')]:
        ds = sorted((t.dist(L[q], L[r])[0], r, t.dist(L[q], L[r])[1]) for r in refs)
        near = ds[:5]
        genus = [r.split('|')[1].split('_')[0] if not r.split('|')[1].startswith(('uncultured', 'Candidatus')) else '_'.join(r.split('|')[1].split('_')[:2]) for _, r, _ in near]
        mrca = near[0][2]; s = t.sup[mrca]
        call = 'OXIDATIVE (reverse-type)' if near[0][1].startswith('REF_OXI') else 'REDUCTIVE'
        rows.append((q.split('|', 1)[1], m, near[0][1].split('|')[1], (near[0][1].split('|') + ['-', '-'])[2], f'{near[0][0]:.3f}', ';'.join(genus),
                     f"{sum(r.startswith('REF_RED') for _, r, _ in near)}/5 reductive", str(s or '-'), call))
with open(f'{OUT}/dsr_placement.tsv', 'w') as o:
    o.write('query\ttree\tnearest_ref_organism\tnearest_ref_acc\tpatristic_dist\tnearest5_genera\tnearest5_reductive\tUFBoot_mrca_query_nearest\ttype_call\n')
    for r in rows: o.write('\t'.join(r) + '\n')
qc.append('QC_STATUS: ' + ('PASS' if ok_all else 'FAIL (placements reported, but the oxidative/reductive split is not recovered)'))
open(f'{OUT}/QC_dsr.txt', 'w').write('\n'.join(qc) + '\n')
print('\n'.join(qc)); print(open(f'{OUT}/dsr_placement.tsv').read())
sys.exit(0 if ok_all else 1)
