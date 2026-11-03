"""SP05 step 4b -- MAG_017 McrA verdict in the mcrA tree, with the step-15 rule (Newick-adjacency MCR vs ACR) and the nearest
reference by patristic distance. CONTROL: the 9 proteins CONFIRMED in 15_marker_trees/mcrA/mcrA_verdicts.tsv stay CONFIRMED."""
import os, re, sys, csv
sys.path.insert(0, os.path.dirname(__file__)); from tree import T
OUT, WS2 = os.environ['OUT'], os.environ['WS2']
tr = open(f'{OUT}/trees/mcrA.treefile').read().strip(); t = T(tr); L = t.leaves()
def adj(leaf, window=400):
    i = tr.find(leaf); seg = tr[max(0, i - window): i + window]
    nb = re.findall(r'(MCR|ACR)\|', seg); c = {g: nb.count(g) for g in set(nb)}
    top = sorted(c, key=lambda g: -c[g]); return ('CONFIRMED' if top and top[0] == 'MCR' else 'REJECTED' if top and top[0] == 'ACR' else 'unclear'), c
refs = [n for n in L if n.startswith(('MCR|', 'ACR|'))]
nom = {r['protein'].replace('+', '_'): r for r in csv.DictReader(open(f'{OUT}/work/mag017_nomination.tsv'), delimiter='\t')}
rows = []
for q in [n for n in L if n.startswith('QUERY')]:
    v, c = adj(q)
    d = sorted((t.dist(L[q], L[r])[0], r) for r in refs)[:5]
    near = d[0][1].split('|')[0]
    pn = q.split('|')[-1]
    rows.append((q, v, ','.join(f'{k}:{n}' for k, n in sorted(c.items())), d[0][1], f'{d[0][0]:.3f}', ','.join(r.split('|')[0] for _, r in d),
                 nom.get(pn, {}).get('K00399_score', '-'), nom.get(pn, {}).get('nomination', 'step-15 original query')))
with open(f'{OUT}/mag017_mcrA_verdict.tsv', 'w') as o:
    o.write('query\tstep15_rule_verdict\tadjacent_groups\tnearest_ref\tpatristic_dist\tnearest5_groups\tK00399_score\tnomination\n')
    for r in rows: o.write('\t'.join(r) + '\n')
orig = [r['query'] for r in csv.DictReader(open(f'{WS2}/15_marker_trees/mcrA/mcrA_verdicts.tsv'), delimiter='\t') if r['verdict'] == 'CONFIRMED']
got = {r[0].split('|', 1)[1]: r[1] for r in rows}
miss = [q for q in orig if got.get(q) != 'CONFIRMED']
ok = len(orig) == 9 and not miss
qc = [f"[{'PASS' if ok else 'FAIL'}] mcrA control: {len(orig) - len(miss)}/{len(orig)} step-15 CONFIRMED proteins reproduced; not reproduced: {miss or 'none'}"]
for r in rows:
    if 'MAG_017' in r[0]: qc.append(f'[RESULT] {r[0]}: {r[1]}; nearest {r[3]} ({r[4]}); K00399 {r[6]} -> {r[7]}')
qc.append('QC_STATUS: ' + ('PASS' if ok else 'FAIL'))
open(f'{OUT}/QC_mcr.txt', 'w').write('\n'.join(qc) + '\n'); print('\n'.join(qc))
sys.exit(0 if ok else 1)
