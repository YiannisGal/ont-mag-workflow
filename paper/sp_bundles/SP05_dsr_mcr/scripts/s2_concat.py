"""SP05 step 2b -- concatenated DsrAB alignment (trimmed A + trimmed B). A taxon is concatenated only when it has exactly one DsrA
and one DsrB: references by organism name, queries by contig. Others stay in the single-gene trees only."""
import os
from collections import defaultdict
OUT = os.environ['OUT']
def rf(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith('>'): n = l[1:].strip(); d[n] = []
        else: d[n].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}
A, Bb = rf(f'{OUT}/trees/dsrA.trim'), rf(f'{OUT}/trees/dsrB.trim')
def key(n): return n.split('|')[1].rsplit('_', 1)[0] if n.startswith('QUERY') else n.split('|')[0] + '|' + n.split('|')[1]
ga, gb = defaultdict(list), defaultdict(list)
for n in A: ga[key(n)].append(n)
for n in Bb: gb[key(n)].append(n)
n = 0
with open(f'{OUT}/trees/dsrAB.trim', 'w') as o, open(f'{OUT}/trees/dsrAB_pairs.tsv', 'w') as t:
    for k in sorted(set(ga) & set(gb)):
        if len(ga[k]) == 1 and len(gb[k]) == 1:
            a, b = ga[k][0], gb[k][0]
            nm = ('QUERY|' + k) if a.startswith('QUERY') else a.split('|')[0] + '|' + a.split('|')[1]
            o.write(f'>{nm}\n{A[a]}{Bb[b]}\n'); t.write(f'{nm}\t{a}\t{b}\n'); n += 1
print('concatenated taxa', n)
