"""SP04 steps 2-4 -- per-OTU GTDB order (sintax >=0.8; BLAST top hit >=90% id over >=90% of the OTU), depth profiles per order,
25-30 cm Thermoplasmata breakdown, MAG 16S links (from SP03), QC.
Depth aggregation (ms v2): counts summed over the extractions of a depth, then relative abundance of all archaeal-library OTUs
(plastid OTUs excluded; none are expected). min/max = per-extraction share. Reporting floor: OTU total >= 100 reads."""
import os, csv, re, sys, hashlib
from collections import defaultdict
B, WS2, OUT = os.environ['B'], os.environ['WS2'], os.environ['OUT']
E = f'{B}/Etoliko_138_2'
DEPTHS = [f'{a}-{a+5}' for a in range(0, 70, 5)]
def dlab(x): a, b = re.split(r'[-_]', x.strip()); return f'{int(a)}-{int(b)}'
bmeta = {r['SampleID']: dlab(r['Depth']) for r in csv.DictReader(open(f'{E}/Bacteria/out_3/qplot_input_fixed/Bacteria Metadata.txt'), delimiter='\t')}
ameta = {k.replace('BA', 'AR'): v for k, v in bmeta.items()}
rows = [l.rstrip('\n').split('\t') for l in open(f'{E}/Archaea/out/otutab_norm.txt')]
hdr = rows[0][1:]; counts = {r[0]: [float(x) for x in r[1:]] for r in rows[1:]}
silva = {r['Feature ID']: r['Taxon'] for r in csv.DictReader(open(f'{E}/Archaea/out/taxonomy_138_2/taxonomy.tsv'), delimiter='\t')}
plast = {o for o, t in silva.items() if 'o__Chloroplast' in t}
bydepth = defaultdict(list)
for i, s in enumerate(hdr): bydepth[ameta[s]].append(i)
tot = [sum(counts[o][i] for o in counts if o not in plast) for i in range(len(hdr))]
gt = dict(l.rstrip('\n').split('\t') for l in open(f'{OUT}/db/taxonomy.tsv'))
RANKS = 'dpcofgs'
# ---- sintax
sx = {}
for l in open(f'{OUT}/work/sintax.tsv'):
    f = l.rstrip('\n').split('\t') + ['', '', '']
    conf = dict((m.group(1), (m.group(2), float(m.group(3)))) for m in re.finditer(r'([dpcofgs]):([^(,]+)\(([\d.]+)\)', f[1]))
    sx[f[0]] = conf
def srank(o, r):
    v = sx.get(o, {}).get(r)
    return v[0] if v and v[1] >= 0.8 else ''
def sconf(o, r):
    v = sx.get(o, {}).get(r); return f'{v[1]:.2f}' if v else '-'
# ---- BLAST top hit
bl = {}
for l in open(f'{OUT}/work/blast.tsv'):
    q, s, pid, ln, ql, sl, ev, bs, qc = l.rstrip('\n').split('\t')
    if q not in bl or float(bs) > bl[q][1]: bl[q] = (s, float(bs), float(pid), float(qc))
def brank(o, r):
    h = bl.get(o)
    if not h or h[2] < 90 or h[3] < 90: return ''
    tax = {x[0]: x[3:] for x in gt.get(h[0], '').split(';') if len(x) > 3}
    return tax.get(r, '')
METH = {'Methanomassiliicoccales', 'Methanomicrobiales', 'Methanosarcinales', 'Methanocellales', 'Methanotrichales',
        'Methanobacteriales', 'Methanococcales', 'Methanonatronarchaeales', 'Methanofastidiosales'}
# ---- per-OTU table
lab = {}
with open(f'{OUT}/archaeal_OTU_orders.tsv', 'w') as o:
    o.write('OTU\ttotal_reads\trank_overall\tGTDB_phylum_sintax\tGTDB_class_sintax\tGTDB_order_sintax\torder_conf\tGTDB_family_sintax\t'
            'GTDB_genus_sintax\tBLAST_top_hit\tBLAST_pid\tBLAST_qcov\tGTDB_order_BLAST\tGTDB_family_BLAST\tagree_order\tmethanogen_order\tSILVA138.2_msv2\n')
    ra = {k: sum(v) for k, v in counts.items() if k not in plast}
    rank = {k: i + 1 for i, (k, _) in enumerate(sorted(ra.items(), key=lambda x: -x[1]))}
    for otu, v in sorted(counts.items(), key=lambda x: -sum(x[1])):
        ordr = srank(otu, 'o') or brank(otu, 'o')
        ph = srank(otu, 'p') or brank(otu, 'p')
        lab[otu] = ordr if ordr else f'Unclassified order ({ph or "unclassified phylum"})'
        if sum(v) < 100: continue
        h = bl.get(otu, ('-', 0, 0, 0))
        so, bo = srank(otu, 'o'), brank(otu, 'o')
        agree = 'yes' if so and bo and so == bo else ('sintax_only' if so and not bo else ('blast_only' if bo and not so else ('both_unassigned' if not so and not bo else 'CONFLICT')))
        o.write('\t'.join([otu, f'{sum(v):.0f}', str(rank.get(otu, 'NA')), srank(otu, 'p'), srank(otu, 'c'), so, sconf(otu, 'o'),
                           srank(otu, 'f'), srank(otu, 'g'), h[0], f'{h[2]:.1f}', f'{h[3]:.0f}', bo, brank(otu, 'f'), agree,
                           'yes' if (so or bo) in METH else 'no', silva.get(otu, 'not in ms v2 trimmed set')]) + '\n')
# ---- order depth profiles
orders = sorted(set(lab.values()))
with open(f'{OUT}/order_depth_profiles.tsv', 'w') as o:
    o.write('order\tmethanogen_order\tdepth_cm\tn_extractions\trel_abund\tmin\tmax\n')
    for od in orders:
        mem = [k for k, v in lab.items() if v == od and k not in plast]
        for d in DEPTHS:
            idx = bydepth[d]
            ssum = sum(sum(counts[k][i] for k in mem) for i in idx); tt = sum(tot[i] for i in idx)
            per = [sum(counts[k][i] for k in mem) / tot[i] for i in idx if tot[i]]
            o.write(f"{od}\t{'yes' if od in METH else 'no'}\t{d}\t{len(idx)}\t{ssum/tt:.6f}\t{min(per):.6f}\t{max(per):.6f}\n")
# ---- Thermoplasmata at 25-30 cm (SILVA class Thermoplasmata in ms v2 labels, or GTDB class Thermoplasmata)
idx = bydepth['25-30']; tt = sum(tot[i] for i in idx)
def share(k): return sum(counts[k][i] for i in idx) / tt
th = [k for k in counts if 'c__Thermoplasmata' in silva.get(k, '') or srank(k, 'c') == 'Thermoplasmata']
bd = defaultdict(float); bdn = defaultdict(int)
for k in th: bd[lab[k]] += share(k); bdn[lab[k]] += 1
thtot = sum(bd.values())
prof = {d: sum(sum(counts[k][i] for i in bydepth[d]) for k in th) / sum(tot[i] for i in bydepth[d]) for d in DEPTHS}
with open(f'{OUT}/thermoplasmata_25_30cm_breakdown.tsv', 'w') as o:
    o.write('GTDB_order\tn_OTUs\tshare_of_archaeal_reads_25_30cm\tfraction_of_Thermoplasmata_25_30cm\n')
    for od, s in sorted(bd.items(), key=lambda x: -x[1]): o.write(f'{od}\t{bdn[od]}\t{s:.5f}\t{s/thtot:.4f}\n')
    o.write('#Thermoplasmata share of archaeal reads by depth: ' + ', '.join(f'{d} {100*v:.1f}%' for d, v in prof.items()) + '\n')
# ---- unclassified fraction per depth + MAG links + QC
qc = []
unc = {d: sum(sum(counts[k][i] for i in bydepth[d]) for k in counts if lab[k].startswith('Unclassified') and k not in plast) /
          sum(tot[i] for i in bydepth[d]) for d in DEPTHS}
qc.append('[REPORT] unclassified at order level (share of archaeal reads) by depth: ' + ', '.join(f'{d} {100*v:.1f}%' for d, v in unc.items()))
p22, c22, o22 = srank('Otu22', 'p') or brank('Otu22', 'p'), srank('Otu22', 'c') or brank('Otu22', 'c'), lab['Otu22']
asg = bool(re.search(r'Asgard|Prometh|Thor|Loki|Heimdall|Odin', p22 + c22))
ok = asg and 'Thermoplasmatota' not in p22
qc.append(f"[{'PASS' if ok else 'FAIL'}] Arc_Otu22 (Otu22) phylum '{p22}' class '{c22}' order '{o22}' (sintax conf p {sconf('Otu22','p')}, c {sconf('Otu22','c')}); BLAST top {bl.get('Otu22', ('-',0,0,0))[0]} {bl.get('Otu22', ('-',0,0,0))[2]:.1f}%")
lk = []
p = f'{WS2}/99_thesis_audit/SP03_16S_links/links.tsv'
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter='\t'):
        if r['dataset'] == 'Illumina_arc_V4V5' and r['link_class'] in ('exact', 'same_genus_level'):
            lk.append((r['genome'], r['OTU'], r['identity'], r['link_class'], lab.get(r['OTU'], '?'), f"{sum(counts.get(r['OTU'], [0])):.0f}"))
with open(f'{OUT}/mag_links_ge97.tsv', 'w') as o:
    o.write('genome\tOTU\tidentity\tlink_class\tOTU_GTDB_order\tOTU_total_reads\n')
    for x in sorted(set(lk)): o.write('\t'.join(x) + '\n')
qc.append(f'[REPORT] methanogen MAG 16S links >=97% to archaeal OTUs (SP03): {sorted(set(lk))}')
qc.append('QC_STATUS: ' + ('PASS' if ok else 'FAIL'))
open(f'{OUT}/QC.txt', 'w').write('\n'.join(qc) + '\n'); print('\n'.join(qc))
md5 = lambda q: hashlib.md5(open(q, 'rb').read()).hexdigest()
import subprocess
v = lambda c: subprocess.run(c, shell=True, capture_output=True, text=True, errors='replace').stdout.strip().splitlines()[0]
with open(f'{OUT}/versions.tsv', 'w') as o:
    o.write('item\tversion_or_md5\tnote\n')
    o.write(f"vsearch\t{v(os.environ['VS'] + ' --version 2>&1')}\tSP03 pixi env; --sintax --sintax_cutoff 0.8\n")
    o.write(f"blastn\t{v('blastn -version')}\t-task blastn\n")
    o.write(f"GTDB SSU\tmd5 {md5(OUT + '/db/ssu_all_r232.fna.gz')}\tssu_all_r232.fna.gz, data.gtdb.ecogenomic.org release232/232.0\n")
    o.write(f"OTU FASTA (primer-trimmed)\tmd5 {md5(os.environ['OTU'])}\t{os.environ['OTU']}\n")
    o.write(f"OTU table\tmd5 {md5(E + '/Archaea/out/otutab_norm.txt')}\tArchaea/out/otutab_norm.txt\n")
    o.write(f"SILVA ms v2 labels\tmd5 {md5(E + '/Archaea/out/taxonomy_138_2/taxonomy.tsv')}\tSILVA 138.2, trimmed OTU set only\n")
sys.exit(0 if ok else 1)
