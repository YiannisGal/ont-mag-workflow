"""SP01 step 5b -- sediment reads per million for ureC and the step-34 controls, step-34 counting rules:
best hit per read (bitscore), label = header prefix; TrEMBL references taking >25% of a label's reads
(label >=100 reads) are excluded as magnets (dataset ASM_ and SwissProt sp| never excluded); DECOY reads are
reported, not counted. Community-level L2 only: a ureC read is not attributed to any genome."""
import os
from collections import Counter, defaultdict
OUT = os.environ['OUT']; U = f'{OUT}/urease'
nreads = int(float(open(f'{U}/SED_nreads.txt').read().split()[0]))
nref = Counter(l[1:].split('|', 1)[0] for l in open(f'{U}/urease_ref.faa') if l.startswith('>'))
best = {}
for l in open(f'{U}/SED_reads.tsv'):
    f = l.rstrip('\n').split('\t')
    q, sid, bs = f[0], f[1], float(f[5])
    if q not in best or bs > best[q][1]: best[q] = (sid.split('|', 1)[0], bs, sid)
per = defaultdict(Counter)
for q, (lab, bs, sid) in best.items(): per[lab][sid] += 1
bad = []
for lab, c in per.items():
    tot = sum(c.values())
    if tot < 100: continue
    for sid, n in c.items():
        if '|tr|' in sid and n / tot > 0.25: bad.append((lab, sid, n, tot))
badset = {b[1] for b in bad}
cnt = Counter(lab for q, (lab, bs, sid) in best.items() if sid not in badset)
asm = Counter(sid for q, (lab, bs, sid) in best.items() if lab == 'ureC' and '|ASM_' in sid)
with open(f'{OUT}/urease_reads.tsv', 'w') as o:
    o.write('sample\tlabel\trole\tn_ref_seqs\treads_best_hit\treads_per_million\tsed_over_rpoB\ttotal_reads\n')
    rpob = cnt.get('CTRLrpoB', 0) or 1
    for lab in ['ureC', 'CTRLrpoB', 'CTRLglnA', 'CTRLtuf', 'DECOYpyrC', 'DECOYhutI', 'DECOYallB']:
        role = 'TARGET' if lab == 'ureC' else ('CONTROL' if lab.startswith('CTRL') else 'DECOY (not counted)')
        n = cnt.get(lab, 0)
        o.write(f'SED\t{lab}\t{role}\t{nref[lab]}\t{n}\t{1e6*n/nreads:.2f}\t{n/rpob:.4f}\t{nreads}\n')
    for b in bad: o.write(f'#excluded_magnet\t{b[0]}\t{b[1]}\t{b[2]} of {b[3]}\n')
    for sid, n in asm.most_common(): o.write(f'#ureC_reads_on_dataset_query\t{sid}\t{n}\n')
print(open(f'{OUT}/urease_reads.tsv').read())
ctrl_ok = all(cnt.get(c, 0) > 0 for c in ('CTRLrpoB', 'CTRLglnA', 'CTRLtuf'))
print('controls detected:', ctrl_ok)
if not ctrl_ok: raise SystemExit('control not detected at read level: the ureC number is uninterpretable')
