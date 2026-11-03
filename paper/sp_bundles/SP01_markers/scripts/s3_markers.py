"""SP01 step 3b -- curated-threshold marker calls for 34 genomes + QC controls.
call = score (ko_list score_type) >= ko_list threshold AND the gene is complete.
LENGTH RULE (user decision 2026-09-29, option 1): the spec's '>=60% of profile length' rule was dropped after it failed the
YGCC01 urease control on full-length UreA (100 aa vs 179-aa profile) and UreG (206 vs 459): KOfam mlen is inflated by fused
homologs. A hit is a 'fragment' only if the gene is truncated: Prodigal partial!=00 (MAGs) or the CDS touches a contig end
(<=3 bp; isolates, Prokka GFF). pct_of_profile is still reported for every hit.
score passes but length fails -> 'fragment' (flagged, never merged, never counted as a call).
Writes markers/calls_detail.tsv (every target hit), markers/calls.tsv (one row per genome x KO),
markers/flat_evalue_inflation.tsv, markers/QC.txt. Exits 1 if any control fails."""
import os, sys, csv, glob, re
from collections import defaultdict
OUT, WS2, B, KOLIST = os.environ['OUT'], os.environ['WS2'], os.environ['B'], os.environ['KOLIST']
MA = f'{B}/UniqueELT_Solid/Manuscript_Analysis'

ko = {}
with open(KOLIST) as fh:
    h = fh.readline().rstrip('\n').split('\t')
    for l in fh:
        f = l.rstrip('\n').split('\t')
        ko[f[0]] = dict(thr=float(f[1]) if f[1] not in ('-', '') else None, st=f[2], mlen=int(f[h.index('mlen')]) if f[h.index('mlen')].isdigit() else 0, defn=f[-1])
targets = [l.rstrip('\n').split('\t') for l in list(open(f'{OUT}/markers/targets.tsv'))[1:]]
gene_of = {}
for k, g in targets: gene_of.setdefault(k, g)

# genomes, completeness, breadth, domain
G = []  # (genome, source, completeness, breadth, domain)
xw = {r['mag_id']: r for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_CROSSWALK.tsv'), delimiter='\t')}
tax = {r['mag_id']: r for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_TAXONOMY.tsv'), delimiter='\t')}
br = {r['MAG']: r['Sediment breadth'] for r in csv.DictReader(open(f'{WS2}/tables/table05_in_situ_detectability.tsv'), delimiter='\t')}
for m in sorted(xw):
    G.append((m, 'MAG', float(xw[m]['completeness']), br.get(m, 'NA'), tax[m]['domain']))
cm2 = {r['Name']: float(r['Completeness']) for r in csv.DictReader(open(f'{MA}/Split_CheckM2/quality_report.tsv'), delimiter='\t')}
ISO_CM2 = {'YGCA01_C': 'YGCA01_Clean', 'YGCA01_Tp': 'YGCA01_Contaminant', 'YGCA07': 'YGCA07', 'YGCA10': 'YGCA10', 'YGCC01': 'YGCC01'}
for g, n in ISO_CM2.items():
    G.append((g, 'isolate', cm2[n], 'NA (isolate)', 'Bacteria'))
GI = {g[0]: g for g in G}
assert len(G) == 34

# protein lengths, contigs
plen = {}; name = None
for l in open(f'{OUT}/markers/all34.faa'):
    if l.startswith('>'): name = l[1:].strip(); plen[name] = 0
    else: plen[name] += len(l.strip())
# truncated genes: MAGs from Prodigal 'partial=' in the step-45 proteome headers
trunc = set()
for l in open(f'{WS2}/45_kofam_full/mags_only.faa'):
    if l.startswith('>'):
        m = re.search(r'partial=(\d\d)', l)
        if m and m.group(1) != '00':
            b, q = l[1:].split()[0].split('|', 1); trunc.add((b, q))
iso_contig = {}; ctg_len = {}
for g, d in (('YGCA01_C', 'YGCA01_C'), ('YGCA01_Tp', 'YGCA01_Contam'), ('YGCA07', 'YGCA07'), ('YGCA10', 'YGCA10'), ('YGCC01', 'YGCC01')):
    for l in open(f'{MA}/Functional_{d}/{d}.gff'):
        if l.startswith('##FASTA'): break
        if l.startswith('##sequence-region'):
            x = l.split(); ctg_len[(g, x[1])] = int(x[3]); continue
        if l.startswith('#'): continue
        f = l.split('\t')
        if len(f) > 8:
            m = re.search(r'locus_tag=([^;\n]+)', f[8])
            if m:
                iso_contig[f'{g}|{m.group(1)}'] = f[0]
                if f[2] == 'CDS' and (int(f[3]) <= 3 or int(f[4]) >= ctg_len.get((g, f[0]), 10**12) - 3):
                    trunc.add((g, m.group(1)))
def contig(p):
    return iso_contig.get(p, '?') if p.split('|')[0] in ISO_CM2 else p.split('|', 1)[1].rsplit('_', 1)[0]

# ---- target hits
xr = {r['working_name']: r['mag_id'] for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_CROSSWALK.tsv'), delimiter='\t')}
xr.update({r['qc_bin_name']: r['mag_id'] for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_CROSSWALK.tsv'), delimiter='\t')})
trunc_ids = {(xr.get(b, b), q) for b, q in trunc}
print('truncated (partial / contig-edge) genes known:', len(trunc_ids))
hits = []
for l in open(f'{OUT}/markers/targets.tbl'):
    if l.startswith('#'): continue
    f = l.split()
    p, k = f[0], f[2]
    i = ko[k]; sc = float(f[8]) if i['st'] == 'domain' else float(f[5])
    sp = i['thr'] is not None and sc >= i['thr']
    lo = (p.split('|')[0], p.split('|', 1)[1]) not in trunc_ids
    call = 'call' if sp and lo else ('fragment' if sp else 'below_threshold')
    hits.append(dict(genome=p.split('|')[0], prot=p, KO=k, score=sc, evalue=float(f[4]), thr=i['thr'], st=i['st'],
                     plen=plen[p], mlen=i['mlen'], call=call))

# ---- genome-wide top passing KO for each target-hit protein (paralog check)
want = {h['prot'] for h in hits if h['call'] != 'below_threshold'}
top = {}
xrev = {}
for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_CROSSWALK.tsv'), delimiter='\t'):
    xrev[r['working_name']] = r['mag_id']; xrev[r['qc_bin_name']] = r['mag_id']
def scan(path, rename):
    for l in open(path, errors='replace'):
        if l.startswith('#'): continue
        f = l.split()
        if len(f) < 9: continue
        p = rename(f[0])
        if p not in want: continue
        i = ko.get(f[2])
        if not i or i['thr'] is None: continue
        sc = float(f[8]) if i['st'] == 'domain' else float(f[5])
        if sc >= i['thr'] and (p not in top or sc / i['thr'] > top[p][1]): top[p] = (f[2], sc / i['thr'])
def mren(x):
    b, q = x.split('|', 1); return f'{xrev.get(b, b)}|{q}'
scan(f'{WS2}/45_kofam_full/all_hits.tbl', mren)
for t in sorted(glob.glob(f'{OUT}/isolates/hits/c*.tbl')): scan(t, lambda x: x)
for h in hits:
    t = top.get(h['prot'])
    h['top_KO'] = t[0] if t else '-'

with open(f'{OUT}/markers/calls_detail.tsv', 'w') as o:
    cols = ['genome', 'prot', 'KO', 'score', 'thr', 'st', 'evalue', 'plen', 'mlen', 'call', 'top_KO']
    o.write('\t'.join(cols) + '\tcontig\n')
    for h in sorted(hits, key=lambda h: (h['genome'], h['KO'], -h['score'])):
        o.write('\t'.join(str(h[c]) for c in cols) + f"\t{contig(h['prot'])}\n")

# ---- one row per genome x KO
by = defaultdict(list)
for h in hits: by[(h['genome'], h['KO'])].append(h)
rows = []
for g in G:
    for k, gene in targets:
        hs = sorted(by.get((g[0], k), []), key=lambda h: -h['score'])
        calls = [h for h in hs if h['call'] == 'call']
        frags = [h for h in hs if h['call'] == 'fragment']
        best = hs[0] if hs else None
        if calls: status = 'present'
        elif frags: status = 'fragment_only'
        elif g[2] >= 90: status = 'absent'
        else: status = 'NA (<90% complete)'
        rows.append(dict(genome=g[0], source=g[1], KO=k, gene=gene, bitscore=f"{best['score']:.1f}" if best else '-',
            threshold=f"{ko[k]['thr']:.2f}" if ko[k]['thr'] is not None else 'none', score_type=ko[k]['st'],
            profile_len=ko[k]['mlen'], passes=status, copies=len(calls), fragments=len(frags),
            contig=';'.join(contig(h['prot']) for h in calls) or '-',
            locus=';'.join(h['prot'].split('|', 1)[1] for h in calls) or '-',
            top_KO_of_called=';'.join(h['top_KO'] for h in calls) or '-',
            completeness=g[2], in_situ_breadth=g[3], domain=g[4], definition=ko[k]['defn']))
with open(f'{OUT}/markers/calls.tsv', 'w') as o:
    o.write('\t'.join(rows[0]) + '\n')
    for r in rows: o.write('\t'.join(str(v) for v in r.values()) + '\n')

# ---- step-16b style flat E-value inflation (E<1e-10, no curated threshold)
with open(f'{OUT}/markers/flat_evalue_inflation.tsv', 'w') as o:
    o.write('KO\tgene\tcurated_calls\tflat_E1e-10_proteins\textra_at_flat_E\tinflation_pct\n')
    tk = tf = 0
    for k, gene in targets:
        hk = [h for h in hits if h['KO'] == k]
        c = sum(h['call'] == 'call' for h in hk); f = sum(h['evalue'] < 1e-10 for h in hk)
        extra = sum(h['evalue'] < 1e-10 and h['call'] != 'call' for h in hk)
        tk += c; tf += extra
        o.write(f"{k}\t{gene}\t{c}\t{f}\t{extra}\t{(100*extra/c) if c else float('nan'):.0f}\n")
    o.write(f"TOTAL\t-\t{tk}\t-\t{tf}\t{100*tf/tk:.0f}\n")

# ---- QC controls
cell = {(r['genome'], r['KO']): r for r in rows}
qc = []; ok = True
arch = [g[0] for g in G if g[4] == 'Archaea']; bact = [g[0] for g in G if g[1] == 'MAG' and g[4] != 'Archaea']
pa = [g for g in arch if cell[(g, 'K00399')]['passes'] == 'present']
pb = [g for g in bact if cell[(g, 'K00399')]['passes'] == 'present']
c1 = len(arch) == 8 and len(bact) == 21 and len(pa) >= 7 and len(pb) == 0
qc.append(f"[{'PASS' if c1 else 'FAIL'}] mcrA: {len(pa)}/{len(arch)} archaea, {len(pb)}/{len(bact)} bacteria  (need >=7/8, 0/21)  archaea without: {sorted(set(arch)-set(pa))}")
ure = ['K01430', 'K01429', 'K01428', 'K03190', 'K03187', 'K03188', 'K03189']
u = {k: cell[('YGCC01', k)]['locus'] for k in ure}
c2 = all(v != '-' and 'JFCCLLDG_' in v for v in u.values())
qc.append(f"[{'PASS' if c2 else 'FAIL'}] YGCC01 urease 7 genes: {u}")
arc = {k: cell[('YGCA01_C', k)]['locus'] for k in ('K01478', 'K00611', 'K00926')}
c3 = {re.sub(r'.*?(BPCFBLIA_\d+).*', r'\1', v) for v in arc.values()} == {'BPCFBLIA_02024', 'BPCFBLIA_02025', 'BPCFBLIA_02026'} \
     and all('BPCFBLIA_0202' in v for v in arc.values())
qc.append(f"[{'PASS' if c3 else 'FAIL'}] YGCA01_C arcA/arcB/arcC: {arc}  (need BPCFBLIA_02024-02026)")
meth = pa
dsr = [g for g in meth for k in ('K11180', 'K11181') if cell[(g, k)]['passes'] == 'present']
c4 = not dsr
qc.append(f"[{'PASS' if c4 else 'FAIL'}] negative control: dsrA/dsrB above threshold in methanogens: {dsr or 'none'}")
ok = c1 and c2 and c3 and c4
qc.append('QC_STATUS: ' + ('PASS' if ok else 'FAIL -- matrix must not be produced'))
open(f'{OUT}/markers/QC.txt', 'w').write('\n'.join(qc) + '\n')
print('\n'.join(qc))
sys.exit(0 if ok else 1)
