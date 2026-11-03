"""SP01 step 2b -- KofamScan-style detail table for the isolates (all hits kept, '*' = passes curated threshold).
Score used is the one ko_list names: 'full' -> full-sequence score, 'domain' -> best-domain score
(what exec_annotation does). Also applies the >=60% of profile length (ko_list mlen) rule."""
import os, glob, statistics
from collections import defaultdict
OUT, KOLIST = os.environ['OUT'], os.environ['KOLIST']

ko = {}
with open(KOLIST) as fh:
    h = fh.readline().rstrip('\n').split('\t')
    for l in fh:
        f = l.rstrip('\n').split('\t')
        t = float(f[1]) if f[1] not in ('-', '') else None
        ko[f[0]] = dict(thr=t, st=f[2], mlen=int(f[h.index('mlen')]) if f[h.index('mlen')].isdigit() else 0, defn=f[-1])

plen = {}; name = None
for l in open(f'{OUT}/isolates/isolates.faa'):
    if l.startswith('>'): name = l[1:].strip(); plen[name] = 0
    else: plen[name] += len(l.strip())

rows = []
tbls = sorted(glob.glob(f'{OUT}/isolates/hits/c*.tbl'))
assert len(tbls) == 24, f'expected 24 chunk tables, got {len(tbls)}'
for tb in tbls:
    for l in open(tb):
        if l.startswith('#'): continue
        f = l.split()
        g, k = f[0], f[2]
        info = ko.get(k)
        if info is None: continue
        full, dom = float(f[5]), float(f[8])
        sc = dom if info['st'] == 'domain' else full
        thr = info['thr']
        pas = thr is not None and sc >= thr
        lenok = plen[g] >= 0.6 * info['mlen']
        rows.append((g.split('|')[0], g, k, thr, sc, f[4], info['st'], plen[g], info['mlen'], pas, lenok, info['defn']))

rows.sort(key=lambda r: (r[1], -r[4]))
with open(f'{OUT}/kofam_isolates_detail.tsv', 'w') as o:
    o.write('#pass\tgenome\tgene_name\tKO\tthrshld\tscore\tE-value\tscore_type\tprot_len\tprofile_len\tlen_ge_60pct\tKO_definition\n')
    for r in rows:
        o.write('\t'.join(['*' if r[9] else '', r[0], r[1], r[2], '' if r[3] is None else f'{r[3]:.2f}',
                           f'{r[4]:.1f}', r[5], r[6], str(r[7]), str(r[8]), 'yes' if r[10] else 'no', r[11]]) + '\n')

per = defaultdict(set)
best = {}
for r in rows:
    if r[9] and (r[1] not in best or r[4] > best[r[1]][1]): best[r[1]] = (r[2], r[4])
for g, (k, _) in best.items(): per[g.split('|')[0]].add(k)
print(f'rows {len(rows)}; genes passing {len(best)}')
for g in sorted(per): print(f'  {g:10s} distinct KOs above threshold: {len(per[g])}')
med = statistics.median(len(v) for v in per.values())
assert len(per) == 5 and med >= 200, f'isolate KOfam not genome-wide (median {med})'
print('GENOME-WIDE OK median', med)
