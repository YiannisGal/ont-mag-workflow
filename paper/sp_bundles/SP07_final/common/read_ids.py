"""Rebuild the read-ID lists behind dsrA 70 / dsrB 68 (step 34, best hit per read + TrEMBL dominance filter,
DOMFRAC 0.25, DOMMIN 100, exactly as scripts/34_unbinned_function.sh) and mcrA 15 (step 31, best hit per read)."""
import sys, os
from collections import Counter, defaultdict
WS2, OUT = sys.argv[1], sys.argv[2]
def best_of(p, minlen=0):
    best = {}
    for line in open(p, errors='replace'):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        try: bs = float(f[5]); ln = int(f[3])
        except ValueError: continue
        if ln < minlen: continue   # step 31: MINLEN 30 aa
        if f[0] not in best or bs > best[f[0]][1]: best[f[0]] = (f[1].split('|', 1)[0], bs, f[1])
    return best
best = best_of(f'{WS2}/34_unbinned/hits/SED_reads.tsv')
per = defaultdict(Counter)
for q, (lab, bs, sid) in best.items(): per[lab][sid] += 1
bad = {sid for lab, c in per.items() if sum(c.values()) >= 100 for sid, n in c.items() if '|tr|' in sid and n / sum(c.values()) > 0.25}
ids = defaultdict(list)
for q, (lab, bs, sid) in best.items():
    if sid not in bad and lab in ('dsrA', 'dsrB'): ids[lab].append((q, sid))
b31 = best_of(f'{WS2}/31_readlevel/SED_read_hits.tsv', minlen=30)
ids['mcrA'] = [(q, sid) for q, (lab, bs, sid) in b31.items() if lab == 'mcrA']
exp = {'dsrA': 70, 'dsrB': 68, 'mcrA': 15}
with open(f'{OUT}/read_ids.tsv', 'w') as o:
    o.write('gene\tread_id\tsource_best_hit\tsource_file\n')
    for g in ('dsrA', 'dsrB', 'mcrA'):
        src = '31_readlevel/SED_read_hits.tsv' if g == 'mcrA' else '34_unbinned/hits/SED_reads.tsv'
        for q, sid in sorted(ids[g]): o.write(f'{g}\t{q}\t{sid}\t{src}\n')
ok = all(len(ids[g]) == exp[g] for g in exp)
msg = ' '.join(f'{g}={len(ids[g])}(expected {exp[g]})' for g in exp) + f'  excluded TrEMBL magnets: {sorted(bad)}'
open(f'{OUT}/read_ids_QC.txt', 'w').write(('PASS ' if ok else 'FAIL ') + msg + '\n'); print(('PASS ' if ok else 'FAIL ') + msg)
sys.exit(0 if ok else 1)
