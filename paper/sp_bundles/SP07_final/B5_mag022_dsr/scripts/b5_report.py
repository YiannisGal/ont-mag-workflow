import os, itertools, math, csv, glob
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
WS2 = os.environ['WS2']
def fa(p):
    d, k = {}, None
    for l in open(p):
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = []
        else: d[k].append(l.strip())
    return {k: ''.join(v).upper() for k, v in d.items()}
K = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
rc = lambda s: ''.join(comp[c] for c in reversed(s))
canon = sorted({min(k, rc(k)) for k in K})
def tnf(s):
    c = dict.fromkeys(canon, 0)
    for i in range(len(s) - 3):
        k = s[i:i + 4]
        if any(x not in 'ACGT' for x in k): continue
        c[min(k, rc(k))] += 1
    t = sum(c.values()) or 1; return [c[k] / t for k in canon]
gc = lambda s: (s.count('G') + s.count('C')) / max(1, sum(s.count(x) for x in 'ACGT'))
dist = lambda a, b: math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
pct = lambda v, q: sorted(v)[min(len(v) - 1, int(q * len(v)))]
mag = fa('work/genomes/MAG022.fa'); tgt = fa('work/contig_19299.fa')['LIQ_contig_19299']
cov = {r['#rname']: float(r['meandepth']) for r in csv.DictReader(open('work/contig_coverage.tsv'), delimiter='\t')}
cent = [sum(x) / len(mag) for x in zip(*[tnf(s) for s in mag.values()])]
own = [dist(tnf(s), cent) for s in mag.values()]; td = dist(tnf(tgt), cent)
g = [gc(s) for s in mag.values()]; cv = [cov[k] for k in mag]
# nearest MAG by TNF centroid among the 29 MAGs
best = []
for p in glob.glob(f'{WS2}/07_mags/bins/LIQ_das_*.fa'):
    m = fa(p); c = [sum(x) / len(m) for x in zip(*[tnf(s) for s in m.values()])]; best.append((dist(tnf(tgt), c), os.path.basename(p)))
best.sort()
cm = {r['Name']: r for r in csv.DictReader(open('work/checkm2/quality_report.tsv'), delimiter='\t')}
gun = {}
for p in glob.glob('work/gunc/GUNC*maxCSS_level.tsv'):
    for r in csv.DictReader(open(p), delimiter='\t'): gun[r['genome']] = r
mk = {}
for p in glob.glob('work/gtdbtk/identify/*bac120_markers_summary.tsv') + glob.glob('work/gtdbtk/identify/gtdbtk.bac120.markers_summary.tsv'):
    for r in csv.DictReader(open(p), delimiter='\t'): mk[r['name']] = r
ln, nm, nc, nmc = open('work/link_stats.txt').read().split()
def dup(name):
    r = mk.get(name, {}); return r.get('multiple_hits', 'NA'), r.get('multiple_unique_genes', r.get('multiple_hits', 'NA'))
# markers that become duplicated: compare per-marker counts if available
dupl = 'see work/gtdbtk/identify'
try:
    import json
    def hits(n):
        out = {}
        for p in glob.glob(f'work/gtdbtk/identify/intermediate_results/marker_genes/{n}/*tophit.tsv'):
            for l in list(open(p))[1:]:
                gene, hs = l.rstrip('\n').split('\t')[:2]
                for h in hs.split(';'):
                    mid = h.split(',')[0]; out[mid] = out.get(mid, 0) + 1
        return out
    a, b = hits('MAG022'), hits('MAG022_plus_19299')
    nd = sorted(k for k in b if b[k] > 1 and a.get(k, 0) <= 1)
    added = sorted(k for k in b if k not in a)
    dupl = f"newly duplicated bac120 markers: {', '.join(nd) if nd else 'none'}; markers added: {', '.join(added) if added else 'none'}"
except Exception as e: dupl = f'not computed ({e})'
c0, c1 = cm.get('MAG022', {}), cm.get('MAG022_plus_19299', {})
rows = [
 ('linking reads LIQ_contig_19299 <-> MAG_022 contig (MAPQ>=20, >=1 kb each side)', f'{ln} reads (19299 linked to any contig: {nm} reads; see B5_links_summary.tsv)', 'yes' if int(ln) > 0 else 'no'),
 ('positive control: MAG_022 contig-contig links', f'{nc} contig pairs bridged by {nmc} reads', 'PASS' if int(nc) > 0 else 'FAIL'),
 ('CheckM2 completeness without -> with', f"{c0.get('Completeness')} -> {c1.get('Completeness')}", 'info'),
 ('CheckM2 contamination without -> with', f"{c0.get('Contamination')} -> {c1.get('Contamination')}", 'yes' if c0 and float(c1['Contamination']) <= float(c0['Contamination']) else 'no'),
 ('GUNC CSS / pass without -> with', f"{gun.get('MAG022',{}).get('clade_separation_score')}/{gun.get('MAG022',{}).get('pass.GUNC')} -> {gun.get('MAG022_plus_19299',{}).get('clade_separation_score')}/{gun.get('MAG022_plus_19299',{}).get('pass.GUNC')}", 'info'),
 ('duplicated single-copy markers (GTDB-Tk bac120)', dupl, 'yes' if 'newly duplicated bac120 markers: none' in dupl else 'check'),
 ('GC', f'{gc(tgt):.3f} vs MAG_022 5-95% {pct(g,.05):.3f}-{pct(g,.95):.3f}', 'yes' if pct(g, .05) <= gc(tgt) <= pct(g, .95) else 'no'),
 ('tetranucleotide distance to MAG_022 centroid', f'{td:.4f} vs own contigs 95th pct {pct(own,.95):.4f}; nearest MAG by TNF: {best[0][1]}', 'yes' if td <= pct(own, .95) else 'no'),
 ('coverage (LIQ reads, mean depth)', f"{cov['LIQ_contig_19299']:.1f}x vs MAG_022 5-95% {pct(cv,.05):.1f}-{pct(cv,.95):.1f}x", 'yes' if pct(cv, .05) <= cov['LIQ_contig_19299'] <= pct(cv, .95) else 'no')]
with open('B5_tests.tsv', 'w') as o:
    o.write('test\tresult\tcompatible\n')
    for r in rows: o.write('\t'.join(r) + '\n')
qc = 'PASS' if int(nc) > 0 else 'FAIL'
linked = int(ln) > 0; cont_ok = rows[3][2] == 'yes'
verdict = 'attributed' if (linked and cont_ok and qc == 'PASS') else ('not attributed' if qc == 'PASS' else 'undecided (positive control failed)')
open('QC.txt', 'w').write(f'positive control (MAG_022 contig-contig links): {nc} pairs, {nmc} reads -> {qc}\nverdict: {verdict}\n')
print(open('B5_tests.tsv').read()); print(open('QC.txt').read())
