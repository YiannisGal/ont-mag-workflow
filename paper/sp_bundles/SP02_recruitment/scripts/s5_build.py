"""SP02 step 5 -- recruitment_by_identity.tsv, identity_histograms.tsv, calls, QC, versions.tsv.
expected breadth = 1 - exp(-mean depth) (Poisson). flag 'pileup' when expected >= 0.20 and observed < 0.5 x expected.
CALIBRATED RULE (user decision 2026-09-29, option 2). The spec's rule (breadth >= 0.50 at >=95%) is biased by Nanopore
read error: the four sediment-binned MAGs (SED_MAG_001-004, i.e. the sediment's own populations) peak at 98-99% read
identity but keep only 73-94% of their reads at >=95%. The library's same-population identity is therefore the median
histogram mode of those four MAGs (REF_MODE), and:
same-species present = histogram mode >= REF_MODE - 2 AND breadth >= 0.50 at >=90% identity
related population   = breadth >= 0.10 at >=90% identity, histogram mode below REF_MODE - 2
The spec rule is still computed and reported (spec_call) for comparison.
not detected         = otherwise (below this library's detection floor)."""
import os, sys, csv, glob, math, subprocess
from collections import defaultdict
OUT, WS2 = os.environ['OUT'], os.environ['WS2']
role = {}; glen = {}
for r in csv.DictReader(open(f'{OUT}/refs/genomes.tsv'), delimiter='\t'):
    role[(r['index'], r['genome'])] = r['role']; glen[(r['index'], r['genome'])] = int(r['length_bp'])
# histograms and modes
mode = {}
with open(f'{OUT}/identity_histograms.tsv', 'w') as o:
    o.write('index\tsample\tgenome\tidentity_bin_pct\treads\n')
    per = defaultdict(dict)
    for p in sorted(glob.glob(f'{OUT}/hist/*.tsv')):
        for r in csv.DictReader(open(p), delimiter='\t'):
            o.write('\t'.join(r.values()) + '\n')
            per[(r['index'], r['sample'], r['genome'])][int(r['identity_bin_pct'])] = int(r['reads'])
    for k, d in per.items(): mode[k] = (max(d, key=d.get), sum(d.values()))
rows = []
for p in sorted(glob.glob(f'{OUT}/cov/*.tsv')):
    ix, s, t = os.path.basename(p)[:-4].split('_'); t = int(t[2:])
    rd = list(csv.reader(open(p), delimiter='\t')); h = rd[0]
    ci = {k: next(i for i, c in enumerate(h) if c.endswith(k)) for k in ('Covered Fraction', ' Mean', 'Trimmed Mean', 'Relative Abundance (%)')}
    for r in rd[1:]:
        g = r[0]
        if g == 'unmapped':
            rows.append(dict(index=ix, genome='unmapped', role='-', sample=s, identity_threshold=t, breadth='NA', mean_depth='NA',
                             trimmed_mean='NA', rel_abundance=float(r[ci['Relative Abundance (%)']]), expected_breadth='NA', flag='-', hist_mode='NA', hist_reads='NA'))
            continue
        b, m, tm, ra = (float(r[ci[k]]) for k in ('Covered Fraction', ' Mean', 'Trimmed Mean', 'Relative Abundance (%)'))
        e = 1 - math.exp(-m)
        flag = 'pileup (observed << Poisson expectation)' if e >= 0.20 and b < 0.5 * e else 'ok'
        md, nr = mode.get((ix, s, g), ('NA', 0))
        rows.append(dict(index=ix, genome=g, role=role[(ix, g)], sample=s, identity_threshold=t, breadth=b, mean_depth=m, trimmed_mean=tm,
                         rel_abundance=ra, expected_breadth=round(e, 4), flag=flag, hist_mode=md, hist_reads=nr))
cols = list(rows[0])
with open(f'{OUT}/recruitment_by_identity.tsv', 'w') as o:
    o.write('\t'.join(cols) + '\n')
    for r in sorted(rows, key=lambda r: (r['index'], r['sample'], r['genome'], r['identity_threshold'])):
        o.write('\t'.join(str(round(v, 5)) if isinstance(v, float) else str(v) for v in r.values()) + '\n')
R = {(r['index'], r['sample'], r['genome'], r['identity_threshold']): r for r in rows}
import statistics
REF = ['ETO_SED_MAG_001', 'ETO_SED_MAG_002', 'ETO_SED_MAG_003', 'ETO_SED_MAG_004']
REF_MODE = statistics.median(mode[('A', 'SED', m)][0] for m in REF)
print('REF_MODE (sediment-binned MAGs, A/SED):', REF_MODE, [mode[('A', 'SED', m)][0] for m in REF])
def call(ix, s, g):
    b95 = R[(ix, s, g, 95)]['breadth']; b90 = R[(ix, s, g, 90)]['breadth']; md = mode.get((ix, s, g), ('NA', 0))[0]
    hi = md != 'NA' and md >= REF_MODE - 2
    if hi and b90 >= 0.50: c = 'same-species present'
    elif b90 >= 0.10: c = 'related population only'
    else: c = 'not detected'
    sc = 'same-species present' if (b95 >= 0.50 and md != 'NA' and md >= 97) else ('related population only' if b90 >= 0.10 else 'not detected')
    return c, sc, b90, b95, R[(ix, s, g, 97)]['breadth'], R[(ix, s, g, 99)]['breadth'], md, R[(ix, s, g, 95)]['mean_depth'], R[(ix, s, g, 95)]['rel_abundance'], R[(ix, s, g, 95)]['flag']
with open(f'{OUT}/isolate_calls.tsv', 'w') as o:
    o.write('genome\tindex\tsample\tcall_calibrated\tspec_call\tbreadth90\tbreadth95\tbreadth97\tbreadth99\thist_mode\tdepth95\trelab95\tflag95\n')
    for g, ix in [('YGCA01_Clean', 'A'), ('YGCA01_Contaminant', 'A'), ('YGCA10', 'A'), ('YGCC01', 'A'), ('YGCA07', 'B')] + [(m, 'A') for m in sorted(r['genome'] for r in csv.DictReader(open(f'{OUT}/refs/genomes.tsv'), delimiter='\t') if r['index'] == 'A' and r['role'] == 'MAG')]:
        for s in ('SED', 'LIQ'):
            o.write(f'{g}\t{ix}\t{s}\t' + '\t'.join(str(round(x, 4)) if isinstance(x, float) else str(x) for x in call(ix, s, g)) + '\n')
# ---- QC
qc = []; ok = True
ms = {r['mag_id']: r for r in csv.DictReader(open(f'{WS2}/47_master/MASTER.tsv'), delimiter='\t')}
ten = [m for m, r in ms.items() if r['in_situ_ten'] == 'yes']
d = [(m, float(ms[m]['sed_breadth']), R[('A', 'SED', m, 0)]['breadth']) for m in ten]   # unfiltered = step-26 method (user decision, option 2)
bad = [(m, a, round(b, 4)) for m, a, b in d if abs(a - b) > 0.03]
c1 = len(ten) == 10 and not bad
qc.append(f"[{'PASS' if c1 else 'FAIL'}] in-situ ten reproduce workspace_2 sediment breadth, unfiltered as in step 26 (+-0.03): "
          f"{'all 10' if c1 else bad}; " + ', '.join(f'{m[-6:]} {a:.3f}->{b:.3f}' for m, a, b in d))
ch = {m: R[('A', 'SED', m, 95)] for m in ('ETO_SED_MAG_002', 'ETO_LIQ_MAG_012')}
qc.append('[REPORT] in-situ ten at the spec 95% threshold / histogram mode: ' + ', '.join(f"{m[-6:]} b95 {R[('A','SED',m,95)]['breadth']:.3f} mode {mode[('A','SED',m)][0]}" for m in ten))
qc.append(f'[REPORT] REF_MODE = {REF_MODE}')
qc.append('[REPORT] Chakrabartyella split (A, SED, 95%): ' + '; '.join(f"{m} breadth {r['breadth']:.3f} depth {r['mean_depth']:.2f} relab {r['rel_abundance']:.3f}%" for m, r in ch.items()))
neg = max(R[('A', 'SED', 'ETO_LIQ_MAG_016', t)]['breadth'] for t in (90, 95, 97, 99))
c3 = neg < 0.02
qc.append(f"[{'PASS' if c3 else 'FAIL'}] negative control MAG_016 max breadth over thresholds {neg:.4f} (< 0.02)")
for ix in ('A', 'B'):
    fs = open(f'{OUT}/bam/{ix}_SED.flagstat').read().splitlines()
    tot = int(next(l for l in fs if 'primary' in l and 'mapped' not in l and 'duplicates' not in l).split()[0])
    mp = int(next(l for l in fs if 'primary mapped' in l).split()[0])
    qc.append(f"[REPORT] index {ix}: sediment primary reads mapped {mp}/{tot} = {100*mp/tot:.2f}% (expected ~2-6%)")
ok = c1 and c3
qc.append('QC_STATUS: ' + ('PASS' if ok else 'FAIL'))
open(f'{OUT}/QC.txt', 'w').write('\n'.join(qc) + '\n'); print('\n'.join(qc))
def v(c):
    return subprocess.run(c, shell=True, capture_output=True, text=True, errors='replace').stdout.strip().splitlines()[0]
with open(f'{OUT}/versions.tsv', 'w') as o:
    o.write('item\tversion\tnote\n')
    o.write(f"minimap2\t{v(os.environ['MM2'] + ' --version')}\t-ax map-ont --secondary=no\n")
    o.write(f"samtools\t{v(os.environ['SAM'] + ' --version')}\t\n")
    o.write(f"coverm\t{v(os.environ['COVERM'] + ' --version')}\t~/.conda/envs/coverm_env (as steps 01/08d/26)\n")
    o.write(f"python\t{v(os.environ['PY'] + ' --version 2>&1')}\tworkspace_2 envs/smetana\n")
    for ix in ('A', 'B'): o.write(f"reference index {ix}\tmd5 {v('md5sum ' + OUT + '/refs/' + ix + '.fa').split()[0]}\trefs/genomes.tsv lists every genome md5\n")
    for l in open(f'{OUT}/refs/reads_md5.tsv'): p, m = l.split(); o.write(f'reads\tmd5 {m}\t{p}\n')
sys.exit(0 if ok else 1)
