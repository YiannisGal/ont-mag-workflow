"""SP03 steps 3-5 -- link classes, depth profiles, Fig 6.2 data, QC.
Link classes (identity = vsearch id; coverage = target (OTU) coverage):
  exact           id >= 99.0 and OTU coverage >= 95%
  same_genus_level 97.0 <= id < 99.0, or id >= 99.0 with OTU coverage < 95%
  related         94.5 <= id < 97.0
  no_link         below 94.5 or no hit
Depth aggregation (ms v2): counts summed over the extractions of a depth, then relative abundance; plastid OTUs
(SILVA o__Chloroplast) removed from every denominator. min/max = per-extraction relative abundance at that depth."""
import os, csv, re, sys
from collections import defaultdict
B, OUT, WS2 = os.environ['B'], os.environ['OUT'], os.environ['WS2']
E = f'{B}/Etoliko_138_2'
def dlab(x):  # '0-05' / '0_5' / '05-10' -> '0-5'
    a, b = re.split(r'[-_]', x.strip()); return f'{int(a)}-{int(b)}'
DEPTHS = [f'{a}-{a+5}' for a in range(0, 70, 5)]

# ---- sample -> depth
bmeta = {r['SampleID']: dlab(r['Depth']) for r in csv.DictReader(open(f'{E}/Bacteria/out_3/qplot_input_fixed/Bacteria Metadata.txt'), delimiter='\t')}
ameta = {k.replace('BA', 'AR'): v for k, v in bmeta.items()}          # user decision: same YG number
ometa = {r['ClusterID']: dlab(r['Depth']) for r in csv.DictReader(open(f'{E}/ONT/final_output/Metadata_minion.txt'), delimiter='\t')}
def taxmap(p):
    t = {}
    for r in csv.DictReader(open(p), delimiter='\t'): t[r['Feature ID']] = r['Taxon']
    return t
DS = {
  'Illumina_bac_V3V4': dict(tab=f'{E}/Bacteria/out_3/otutab_norm.txt', tax=taxmap(f'{E}/Bacteria/out_3/taxonomy_138_2/taxonomy.tsv'), meta=bmeta, hits='hits/V3V4_vs_bac.tsv'),
  'Illumina_arc_V4V5': dict(tab=f'{E}/Archaea/out/otutab_norm.txt', tax=taxmap(f'{E}/Archaea/out/taxonomy_138_2/taxonomy.tsv'), meta=ameta, hits='hits/V4V5_vs_arc.tsv'),
  'ONT_full':          dict(tab=f'{E}/ONT/final_output/cluster_tab_trimmed.txt', tax=taxmap(f'{E}/ONT/final_output/taxonomy.tsv'), meta=ometa, hits='hits/full_vs_ont.tsv'),
}
for name, d in DS.items():
    rows = [l.rstrip('\n').split('\t') for l in open(d['tab']) if not l.startswith('# Constructed')]
    hdr = rows[0][1:]
    miss = [s for s in hdr if s not in d['meta']]
    if miss: sys.exit(f'{name}: samples without a depth in metadata: {miss}')
    counts = {r[0]: [float(x) for x in r[1:]] for r in rows[1:]}
    plast = {o for o in counts if 'o__Chloroplast' in d['tax'].get(o, '')}
    tot = [sum(counts[o][i] for o in counts if o not in plast) for i in range(len(hdr))]
    bydepth = defaultdict(list)
    for i, s in enumerate(hdr): bydepth[d['meta'][s]].append(i)
    d.update(hdr=hdr, counts=counts, plast=plast, tot=tot, bydepth=bydepth)
    ra_all = {o: sum(v) / sum(tot) for o, v in counts.items() if o not in plast}
    d['rank'] = {o: i + 1 for i, (o, _) in enumerate(sorted(ra_all.items(), key=lambda x: -x[1]))}
    d['ra_all'] = ra_all
    print(f'{name}: {len(counts)} features, {len(hdr)} samples, plastid excluded {len(plast)}, '
          f'depths {sorted(bydepth, key=lambda x: int(x.split("-")[0]))}, per-depth n {[len(bydepth[k]) for k in DEPTHS if k in bydepth]}')

# ---- links
def cls(idv, tcov):
    if idv >= 99.0: return 'exact' if tcov >= 95.0 else 'same_genus_level'
    if idv >= 97.0: return 'same_genus_level'
    if idv >= 94.5: return 'related'
    return 'no_link'
gsrc = {l.split('\t')[0]: (l.split('\t')[1], l.split('\t')[3], int(l.split('\t')[6])) for l in list(open(f'{OUT}/rrna/genomes.tsv'))[1:]}
links = []
for name, d in DS.items():
    seen_q = set()
    for l in open(f'{OUT}/{d["hits"]}'):
        f = l.rstrip('\n').split('\t')
        q, t = f[0], f[1]
        g, copy = q.split('|')[0], q.split('|')[1]
        seen_q.add(q)
        if t == '*':
            links.append(dict(genome=g, source=gsrc[g][0], copy=copy, dataset=name, OTU='-', identity='-', aln_len='-', qcov='-', tcov='-', link_class='no_hit_>=90%', taxon=''))
            continue
        idv, tcov = float(f[2]), float(f[7])
        links.append(dict(genome=g, source=gsrc[g][0], copy=copy, dataset=name, OTU=t, identity=f'{idv:.1f}', aln_len=f[3],
                          qcov=f[6], tcov=f[7], link_class=cls(idv, tcov), taxon=d['tax'].get(t, '')))
    # copies with no region (primer not found) are reported, not silently dropped
    qfile = {'Illumina_bac_V3V4': 'regions/V3V4.fa', 'Illumina_arc_V4V5': 'regions/V4V5.fa', 'ONT_full': 'rrna/genome_16S.fa'}[name]
    have = {l[1:].split()[0] for l in open(f'{OUT}/{qfile}') if l.startswith('>')}
    for l in open(f'{OUT}/rrna/genome_16S.fa'):
        if l.startswith('>'):
            q = l[1:].strip()
            if q not in have:
                links.append(dict(genome=q.split('|')[0], source=gsrc[q.split('|')[0]][0], copy=q.split('|')[1], dataset=name, OTU='-',
                                  identity='-', aln_len='-', qcov='-', tcov='-', link_class='no_primer_region', taxon=''))
for g, (src, taxon, n) in gsrc.items():
    if n == 0:
        for name in DS:
            links.append(dict(genome=g, source=src, copy='-', dataset=name, OTU='-', identity='-', aln_len='-', qcov='-', tcov='-',
                              link_class='unlinkable (no 16S in genome)', taxon=''))
cols = ['genome', 'source', 'copy', 'dataset', 'OTU', 'identity', 'aln_len', 'qcov', 'tcov', 'link_class', 'taxon']
with open(f'{OUT}/links.tsv', 'w') as o:
    o.write('genome\tsource\tcopy\tdataset\tOTU\tidentity\taln_len\tquery_cov\tOTU_cov\tlink_class\tOTU_taxonomy_SILVA138.2\n')
    for r in sorted(links, key=lambda r: (r['genome'], r['dataset'], r['copy'], r['OTU'])): o.write('\t'.join(str(r[c]) for c in cols) + '\n')

# ---- depth profiles for every linked OTU (exact / genus / related)
linked = {(r['dataset'], r['OTU']) for r in links if r['link_class'] in ('exact', 'same_genus_level', 'related')}
prof = {}
with open(f'{OUT}/depth_profiles.tsv', 'w') as o:
    o.write('OTU\tdataset\tdepth_cm\tn_extractions\trel_abund_mean\tmin\tmax\toverall_rank\toverall_rel_abund\tplastid\n')
    for ds, otu in sorted(linked):
        d = DS[ds]; v = d['counts'][otu]
        for dep in DEPTHS:
            idx = d['bydepth'].get(dep, [])
            if not idx: continue
            mean = sum(v[i] for i in idx) / sum(d['tot'][i] for i in idx)
            per = [v[i] / d['tot'][i] for i in idx if d['tot'][i] > 0]
            prof[(ds, otu, dep)] = mean
            o.write(f"{otu}\t{ds}\t{dep}\t{len(idx)}\t{mean:.6f}\t{min(per):.6f}\t{max(per):.6f}\t{d['rank'].get(otu, 'NA')}\t"
                    f"{d['ra_all'].get(otu, 0):.6f}\t{'yes' if otu in d['plast'] else 'no'}\n")

# ---- Figure 6.2 data: exact links only; one row per genome x dataset x OTU
with open(f'{OUT}/figure_data_Fig6_2.tsv', 'w') as o:
    o.write('genome\tdataset\tOTU\tcopies_linked\t' + '\t'.join(f'{d}cm' for d in DEPTHS) + '\n')
    ex = defaultdict(set)
    for r in links:
        if r['link_class'] == 'exact': ex[(r['genome'], r['dataset'], r['OTU'])].add(r['copy'])
    for (g, ds, otu), cps in sorted(ex.items()):
        o.write(f'{g}\t{ds}\t{otu}\t{len(cps)}\t' + '\t'.join(f"{prof.get((ds, otu, dep), float('nan')):.6f}" for dep in DEPTHS) + '\n')

# ---- QC
qc = []; ok = True
def ctl(p):
    return [l.rstrip('\n').split('\t') for l in open(f'{OUT}/{p}')]
cf = ctl('hits/control_full_vs_arc.tsv'); cr = ctl('hits/control_V4V5_vs_arc.tsv')
def best(rows, otu):
    x = [float(r[2]) for r in rows if r[1] == otu]; return max(x) if x else None
c22f, c22r, c674 = best(cf, 'Otu22'), best(cr, 'Otu22'), best(cf, 'Otu674') or best(cr, 'Otu674')
top = sorted(cr, key=lambda r: -float(r[2]) if r[2] not in ('*', '') and r[1] != '*' else 0)[:1]
# user decision 2026-09-29: accept 97.1% (supplied gene: 407/417; ms v2 reported 411/415). PASS = Otu22 is the top hit at >=97%
# over the 415-position primer-free region, and Otu674 is also returned.
c1 = (c22r is not None and c22r >= 97.0) and bool(top) and top[0][1] == 'Otu22' and c674 is not None
qc.append(f"[{'PASS' if c1 else 'FAIL'}] Arc_Otu22 control: Otu22 {c22r}% (Arc519F/915R region), {c22f}% (full gene); Otu674 hit {c674}%; "
          f"best region hit {top[0][1] if top else '-'} (ms v2: 411/415 = 99.0%)")
meth = {g for g, (src, taxon, n) in gsrc.items() if taxon.startswith('Methano') or taxon == 'MTP4'}
bad = [r for r in links if r['dataset'] == 'Illumina_bac_V3V4' and r['genome'] in meth and r['link_class'] in ('exact', 'same_genus_level', 'related')
       and 'd__Bacteria' in r['taxon']]
# user decision 2026-09-29: 341F/805R also amplifies archaea, so archaeal/unassigned OTUs in the bacterial dataset are expected;
# only a methanogen matching a SILVA-bacterial OTU stops the run. The expected ones are reported, not hidden.
arcin = sorted({(r['genome'], r['OTU'], r['identity'], r['taxon'] or 'no SILVA assignment (outside trimmed set)') for r in links
                if r['dataset'] == 'Illumina_bac_V3V4' and r['genome'] in meth and r['link_class'] in ('exact', 'same_genus_level', 'related')})
qc.append(f"[INFO] methanogen matches to archaeal/unassigned OTUs in the bacterial V3-V4 dataset (expected, universal primers): {arcin}")
c2 = not bad
qc.append(f"[{'PASS' if c2 else 'FAIL'}] no methanogen links to a SILVA-bacterial V3-V4 OTU: {[(r['genome'], r['OTU'], r['identity']) for r in bad] or 'none'}")
d = DS['Illumina_bac_V3V4']; idx = d['bydepth']['0-5']
top05 = sorted(((sum(d['counts'][o][i] for i in idx), o) for o in d['counts'] if o not in d['plast']), reverse=True)[:5]
lk = defaultdict(set)
for r in links:
    if r['dataset'] == 'Illumina_bac_V3V4' and r['link_class'] in ('exact', 'same_genus_level'): lk[r['OTU']].add(r['genome'])
det = [(o, d['tax'].get(o, '?').split('; ')[-3:], sorted(lk.get(o, []))) for _, o in top05]
c3 = not any(g in meth for _, _, gs in det for g in gs)
qc.append(f"[{'PASS' if c3 else 'FAIL'}] 0-5 cm top-5 bacterial OTUs and genomes linked (>=97%): {det}")
ok = c1 and c2 and c3
qc.append('QC_STATUS: ' + ('PASS' if ok else 'FAIL'))
open(f'{OUT}/QC.txt', 'w').write('\n'.join(qc) + '\n'); print('\n'.join(qc))
sys.exit(0 if ok else 1)
