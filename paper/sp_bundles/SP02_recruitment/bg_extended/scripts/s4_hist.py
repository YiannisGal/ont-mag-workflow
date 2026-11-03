"""SP02 step 4 -- per-genome read-identity histogram (1% bins) from primary alignments with >=50% of the read aligned.
identity = (aligned columns - NM) / aligned columns, aligned columns = M/=/X + I + D (CoverM's definition of percent identity)."""
import os, sys, re, subprocess
from collections import defaultdict
OUT, SAM = os.environ['OUT'], os.environ['SAM']
ix, s = sys.argv[1], sys.argv[2]
cig = re.compile(r'(\d+)([MIDNSHP=X])')
h = defaultdict(lambda: [0] * 101)
p = subprocess.Popen([SAM, 'view', '-F', '0x904', '-@', '4', f'{OUT}/bam/{ix}_{s}.bam'], stdout=subprocess.PIPE, text=True)
for l in p.stdout:
    f = l.split('\t', 12)
    al = qlen = 0
    for n, op in cig.findall(f[5]):
        n = int(n)
        if op in 'M=XID': al += n
        if op in 'M=XIS': qlen += n
        if op == 'H': qlen += n
    if not al or qlen == 0: continue
    alq = sum(int(n) for n, op in cig.findall(f[5]) if op in 'M=XI')
    if alq / qlen < 0.5: continue
    m = re.search(r'\tNM:i:(\d+)', l)
    if not m: continue
    idt = 100.0 * (al - int(m.group(1))) / al
    h[f[2].split('~', 1)[0]][min(100, max(0, int(idt)))] += 1
if p.wait() != 0: sys.exit('samtools view failed')
with open(f'{OUT}/hist/{ix}_{s}.tsv.part', 'w') as o:
    o.write('index\tsample\tgenome\tidentity_bin_pct\treads\n')
    for g in sorted(h):
        for b, n in enumerate(h[g]):
            if n: o.write(f'{ix}\t{s}\t{g}\t{b}\t{n}\n')
os.rename(f'{OUT}/hist/{ix}_{s}.tsv.part', f'{OUT}/hist/{ix}_{s}.tsv')
print(ix, s, 'genomes with reads:', len(h))
