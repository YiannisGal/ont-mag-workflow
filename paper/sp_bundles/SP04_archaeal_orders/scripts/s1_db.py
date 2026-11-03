"""SP04 step 1 -- GTDB ssu_all_r232 -> vsearch sintax database and plain BLAST FASTA.
GTDB header: >RS_GCF_000001.1~NZ_xxx d__Archaea;p__...;s__Genus species [location=...] [ssu_len=...] ...
sintax header: >id;tax=d:Archaea,p:...,c:...,o:...,f:...,g:...,s:...;   (commas/colons inside names replaced)"""
import gzip, os, re
OUT = os.environ['OUT']
n = 0
with gzip.open(f'{OUT}/db/ssu_all_r232.fna.gz', 'rt') as fh, open(f'{OUT}/db/gtdb_sintax.fa.part', 'w') as s, \
        open(f'{OUT}/db/gtdb_plain.fa', 'w') as p, open(f'{OUT}/db/taxonomy.tsv', 'w') as t:
    seq = []; hid = tax = None
    def flush():
        global n
        if hid is None: return
        sq = ''.join(seq)
        ranks = [x.strip() for x in tax.split(';')]
        st = ','.join(f"{r[0]}:{re.sub(r'[,:; ]', '_', r[3:])}" for r in ranks if len(r) > 3)
        s.write(f'>{hid};tax={st};\n{sq}\n'); p.write(f'>{hid}\n{sq}\n'); t.write(f'{hid}\t{tax}\n'); n += 1
    for l in fh:
        if l.startswith('>'):
            flush(); seq = []
            parts = l[1:].rstrip('\n').split(' ', 1)
            hid = parts[0].replace('~', '__')
            tax = re.split(r' \[', parts[1])[0]
        else: seq.append(l.strip())
    flush()
os.rename(f'{OUT}/db/gtdb_sintax.fa.part', f'{OUT}/db/gtdb_sintax.fa')
print('GTDB SSU sequences:', n)
