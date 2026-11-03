import sys, os, csv, glob
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); D = os.path.join(P, sys.argv[1]); os.chdir(P)
WS2 = os.environ['WS2']
def tsv(p): return list(csv.DictReader(open(p), delimiter='\t'))
if not glob.glob(f'{D}/gtdbtk/gtdbtk.bac120.summary.tsv'):
    open(f'{D}/QC.txt', 'w').write('INVALID: GTDB-Tk bac120 summary missing (bacterial classification failed)\n'); sys.exit('bac120 summary missing')
cm = {r['Name']: r for r in tsv(f'{D}/checkm2/quality_report.tsv')} if os.path.exists(f'{D}/checkm2/quality_report.tsv') else {}
gu = {}
for p in glob.glob(f'{D}/gunc/GUNC*maxCSS*.tsv'):
    for r in tsv(p): gu[r['genome']] = r
gt = {}
for p in glob.glob(f'{D}/gtdbtk/gtdbtk.*.summary.tsv'):
    for r in tsv(p): gt[r['user_genome']] = r['classification']
rank = lambda c, k: next((x[3:] for x in c.split(';') if x.startswith(k + '__')), '') if c else ''
TGT = ('Chloroflexota', 'Desulfobacterota')
rows = []
for b in sorted(cm):
    c = gt.get(b, 'unclassified'); ph = rank(c, 'p'); g = rank(c, 'g')
    comp, cont = float(cm[b]['Completeness']), float(cm[b]['Contamination'])
    smp = 'SED' if b.startswith('SED') else ('LIQ' if b.startswith('LIQ') else 'reassembly')
    rows.append(dict(bin=b, sample=smp, phylum=ph or 'none', genus=g or 'none', classification=c, completeness=comp, contamination=cont,
                     GUNC=f"{gu.get(b,{}).get('clade_separation_score','NA')}/{gu.get(b,{}).get('pass.GUNC','NA')}",
                     MIMAG_medium='yes' if comp >= 50 and cont < 10 else 'no', target_phylum='yes' if ph.startswith(TGT) else 'no'))
k = ['bin', 'sample', 'phylum', 'genus', 'completeness', 'contamination', 'GUNC', 'MIMAG_medium', 'target_phylum', 'classification']
with open(f'{D}/B4_bins.tsv', 'w') as o:
    o.write('\t'.join(k) + '\n')
    for r in rows: o.write('\t'.join(str(r[x]) for x in k) + '\n')
tg = [r for r in rows if r['target_phylum'] == 'yes']
with open(f'{D}/B4_target_bins.tsv', 'w') as o:
    o.write('\t'.join(k) + '\n')
    for r in tg: o.write('\t'.join(str(r[x]) for x in k) + '\n')
hit = [r for r in tg if r['MIMAG_medium'] == 'yes' and r['GUNC'].endswith('True')]
open(f'{D}/eval_flag.txt', 'w').write(''.join(f"STEP2_HIT_{x} " + ','.join(r['bin'] for r in hit if r['sample'] == x) + '\n' for x in ('SED', 'LIQ', 'reassembly') if any(r['sample'] == x for r in hit)) or 'STEP2_NONE\n')
msg = [f'{sys.argv[1]}: {len(rows)} bins (SED {sum(r["sample"]=="SED" for r in rows)}, LIQ {sum(r["sample"]=="LIQ" for r in rows)}); Chloroflexota/Desulfobacterota bins: {len(tg)}; of these MIMAG medium + GUNC pass: {len(hit)}']
if sys.argv[1] == 'step1':
    pub = [r for r in tsv(f'{WS2}/07_mags/MAG_TAXONOMY.tsv') if r['sample'] == 'SED']
    ok = 0
    for m in pub:
        mt = [r for r in rows if r['genus'] == m['genus'] and abs(r['completeness'] - float(m['completeness'])) <= 5]
        ok += bool(mt)
        msg.append(f"  {m['mag_id']} ({m['genus'] or m['family']}, {m['completeness']}%): {'reproduced by ' + mt[0]['bin'] if mt else 'NOT reproduced'}")
    msg.insert(1, f'reproduction QC: {ok}/{len(pub)} published SED MAGs reproduced (same GTDB genus, completeness +-5) -> {"PASS" if ok == len(pub) else "FAIL"}')
open(f'{D}/QC.txt', 'w').write('\n'.join(msg) + '\n'); print('\n'.join(msg))
