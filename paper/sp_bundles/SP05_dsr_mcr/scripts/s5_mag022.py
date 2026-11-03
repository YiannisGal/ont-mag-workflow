"""SP05 step 3 -- could the dsr contigs belong to ETO_LIQ_MAG_022 (LIQ_das_metabat_022, JBAIYS01, Desulfohalotomaculaceae)?
Tests (each: compatible / incompatible):
  GC        contig GC inside the 5-95% range of MAG_022 contigs (>=2.5 kb, length-weighted)
  TNF       distance of the contig's tetranucleotide profile to the MAG_022 centroid <= 95th percentile of MAG_022's own contigs,
            and MAG_022 is the nearest of the 29 MAG centroids
  coverage  LIQ mean depth (step 58 coverage.tsv, same BAM for all contigs) inside the 5-95% range of MAG_022 contigs
  markers   no KOfam KO above threshold on the contig that MAG_022 already carries (redundancy = incompatible)
Verdict "probably belongs to MAG_022" only if all four are compatible."""
import os, csv, itertools, math, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__)); import kofam
OUT, WS2 = os.environ['OUT'], os.environ['WS2']
def rf(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        else: d[n].append(l.strip().upper())
    return {k: ''.join(v) for k, v in d.items()}
comp = str.maketrans('ACGT', 'TGCA')
K = sorted({min(''.join(k), ''.join(k).translate(comp)[::-1]) for k in itertools.product('ACGT', repeat=4)})
KI = {k: i for i, k in enumerate(K)}
def tnf(s):
    v = [1.0] * len(K)
    for i in range(len(s) - 3):
        w = s[i:i + 4]
        if 'N' in w or len(set(w) - set('ACGT')): continue
        v[KI[min(w, w.translate(comp)[::-1])]] += 1
    t = sum(v); return [x / t for x in v]
def d(a, b): return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
def gc(s): return (s.count('G') + s.count('C')) / max(1, len(s) - s.count('N'))
def q(v, p): v = sorted(v); return v[min(len(v) - 1, int(p * (len(v) - 1)))]
xw = {r['mag_id']: r['working_name'] for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_CROSSWALK.tsv'), delimiter='\t')}
cent = {}; m22 = None
for mid, wn in xw.items():
    c = {k: v for k, v in rf(f'{WS2}/07_mags/bins/{wn}.fa').items() if len(v) >= 2500}
    T = {k: tnf(v) for k, v in c.items()}
    cent[mid] = [sum(x) / len(T) for x in zip(*T.values())]
    if mid == 'ETO_LIQ_MAG_022': m22 = (c, T)
c22, T22 = m22
own = [d(v, cent['ETO_LIQ_MAG_022']) for v in T22.values()]
gcs = [gc(v) for v in c22.values()]
cov = {}
for r in csv.DictReader(open(f'{WS2}/58_dsr_finish/coverage.tsv'), delimiter='\t'):
    cov[r['#rname']] = float(r['meandepth'])
cov22 = [v for k, v in cov.items() if k.startswith(xw['ETO_LIQ_MAG_022'] + '~')]
H = kofam.load()
kos22 = set()
import glob
kl = {}
for l in open(os.environ['KOLIST']):
    f = l.rstrip('\n').split('\t')
    if f[0].startswith('K'):
        try: kl[f[0]] = (float(f[1]), f[2])
        except ValueError: pass
for l in open(f'{WS2}/45_kofam_full/all_hits.tbl'):
    f = l.split()
    if not f[0].startswith(xw['ETO_LIQ_MAG_022'] + '|'): continue
    t = kl.get(f[2])
    if t and (float(f[8]) if t[1] == 'domain' else float(f[5])) >= t[0]: kos22.add(f[2])
ctg = rf(f'{OUT}/work/contigs.fa')
rows = []
for c, s in ctg.items():
    v = tnf(s); dd = d(v, cent['ETO_LIQ_MAG_022'])
    nearest = min(cent, key=lambda m: d(v, cent[m]))
    ckos = {kofam.best_pass(H, p)[0] for p in H if p.startswith(c + '_') and kofam.best_pass(H, p)}
    shared = sorted(ckos & kos22)
    tests = [
      ('GC', f'{gc(s):.3f}', f'MAG_022 5-95%: {q(gcs,.05):.3f}-{q(gcs,.95):.3f}', q(gcs, .05) <= gc(s) <= q(gcs, .95)),
      ('TNF', f'dist {dd:.4f}; nearest MAG centroid {nearest}', f'MAG_022 own-contig 95th pct {q(own,.95):.4f}', dd <= q(own, .95) and nearest == 'ETO_LIQ_MAG_022'),
      ('coverage_LIQ', f"{cov.get('UNBINNED~' + c, float('nan')):.2f}x", f'MAG_022 5-95%: {q(cov22,.05):.2f}-{q(cov22,.95):.2f}x (n={len(cov22)})',
       q(cov22, .05) <= cov.get('UNBINNED~' + c, -1) <= q(cov22, .95)),
      ('marker_redundancy', f"contig KOs {sorted(ckos)}", f"already in MAG_022: {shared or 'none'}", not shared),
    ]
    allok = all(t[3] for t in tests)
    for name, obs, ref, ok in tests: rows.append((c, name, obs, ref, 'compatible' if ok else 'INCOMPATIBLE'))
    rows.append((c, 'VERDICT', '-', '-', 'probably belongs to ETO_LIQ_MAG_022' if allok else 'not attributable to ETO_LIQ_MAG_022'))
with open(f'{OUT}/mag022_compatibility.tsv', 'w') as o:
    o.write('contig\ttest\tobserved\treference\tresult\n')
    for r in rows: o.write('\t'.join(r) + '\n')
print(open(f'{OUT}/mag022_compatibility.tsv').read())
