"""SP05 step 4a -- step-15 nomination for MAG_017 using the step-54 fragment-joined McrA: KOfam K00399 score of each joined
protein vs the curated threshold, then the step-15 tree input (15_marker_trees/mcrA/mcrA_all.faa: 60 MCR + 36 ACR refs +
the 9 original queries) plus MAG_017's joined proteins."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__)); import kofam
OUT, WS2 = os.environ['OUT'], os.environ['WS2']
H = kofam.load()
rows = []
for l in open(f'{OUT}/work/mag017_mcrA.faa'):
    if l.startswith('>'):
        n = l[1:].strip()
        s = {k: (sc, t) for k, sc, t, ok, d in H.get(n, [])}
        sc, t = s.get('K00399', (0.0, None))
        rows.append((n, f'{sc:.1f}', str(t), 'passes' if t and sc >= t else 'below threshold'))
open(f'{OUT}/work/mag017_nomination.tsv', 'w').write('protein\tK00399_score\tthreshold\tnomination\n' + ''.join('\t'.join(r) + '\n' for r in rows))
print('\n'.join('\t'.join(r) for r in rows))
seqs = {}; n = None
for l in open(f'{OUT}/work/mag017_mcrA.faa'):
    if l.startswith('>'): n = l[1:].strip(); seqs[n] = ''
    else: seqs[n] += l.strip()
with open(f'{OUT}/trees/mcrA_in.faa', 'w') as o:
    o.write(open(f'{WS2}/15_marker_trees/mcrA/mcrA_all.faa').read())
    for n, s in seqs.items():
        if 'joined' in n: o.write(f'>QUERY|ETO_LIQ_MAG_017|{n}\n{s}\n')
