"""SP01 step 4c -- hydrogenase catalytic subunits classified against the HydDB-labelled reference sets
(Greening lab metabolic marker databases, figshare 13152518 [NiFe] / 13152542 [FeFe], 2020-12-06; labels from HydDB).
HydDB web (services.birc.au.dk/hyddb) did not respond on 2026-09-29 and METABOLIC is not installed.
Rule (Greening-lab practice): DIAMOND blastp best hit; group assigned at >=50% identity and >=50% query coverage;
40-50% identity is reported as 'tentative'; below 40% not reported."""
import os, re, subprocess, csv
from collections import defaultdict
OUT = os.environ['OUT']; H = f'{OUT}/hyd'; D = os.environ['DIAMOND']; CPU = os.environ['CPU']
with open(f'{H}/hyd_refs.faa', 'w') as o:
    for f in ('NiFe_hydrogenase.fasta', 'FeFe_hydrogenase.fasta'):
        for l in open(f'{H}/{f}'):
            if l.startswith('>'):
                m = re.match(r'(\S+)\s+(?:-\s+)?(.*?)\s+-\s+\[(NiFe|FeFe)\] Group (\S+)\s*$', l[1:].rstrip('\n'))
                assert m, l                       # 5 [FeFe] headers lack the first ' - '; the regex takes both forms
                acc, org, cls, grp = m.groups()
                o.write(f">{acc}|[{cls}]_Group_{grp}|{org.replace(' ', '_')}\n")
            else: o.write(l)
def run(cmd):
    print('CMD', ' '.join(cmd), flush=True); subprocess.run(cmd, check=True)
run([D, 'makedb', '--in', f'{H}/hyd_refs.faa', '-d', f'{H}/hyd_refs', '--quiet'])
run([D, 'blastp', '-d', f'{H}/hyd_refs.dmnd', '-q', f'{OUT}/markers/all34.faa', '-o', f'{H}/hits.tsv', '--more-sensitive',
     '-e', '1e-10', '-k', '5', '--threads', CPU, '--quiet', '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
     'evalue', 'bitscore', 'qcovhsp', 'scovhsp'])
best = {}
for l in open(f'{H}/hits.tsv'):
    f = l.rstrip('\n').split('\t')
    q, bs = f[0], float(f[7])
    if q not in best or bs > best[q][0]: best[q] = (bs, f)
rows = []
for q, (bs, f) in best.items():
    pid, qcov = float(f[2]), float(f[8])
    if pid < 40 or qcov < 50: continue
    acc, grp, org = f[1].split('|')
    cls = grp.split('_Group_')[0]; g = grp.split('_Group_')[-1] if '_Group_' in grp else grp
    conf = 'assigned' if pid >= 50 else 'tentative (40-50% id)'
    rows.append((q.split('|')[0], q.split('|', 1)[1], cls, g, acc, org, f'{pid:.1f}', f'{qcov:.0f}', f[4], conf))
rows.sort()
with open(f'{OUT}/hydrogenase_classes.tsv', 'w') as o:
    o.write('genome\tprotein\tclass\tgroup\tbest_ref\tbest_ref_org\tpident\tqcov\tqlen\tconfidence\n')
    for r in rows: o.write('\t'.join(r) + '\n')
per = defaultdict(set)
for r in rows:
    if r[9] == 'assigned': per[r[0]].add(f'{r[2]} {r[3]}')
for g in sorted(per): print(g, sorted(per[g]))
print('hydrogenase subunits reported:', len(rows), 'assigned:', sum(r[9] == 'assigned' for r in rows))
