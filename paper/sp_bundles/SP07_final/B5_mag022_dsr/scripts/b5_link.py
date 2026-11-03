import re, collections, os
os.chdir(os.path.join(os.path.dirname(__file__), '..', 'work'))
mag = set(l.strip() for l in open('mag022_contigs.txt'))
def span(cig): return sum(int(n) for n, op in re.findall(r'(\d+)([MDN=X])', cig))
parts = collections.defaultdict(dict)   # read -> contig -> max ref span (MAPQ>=20)
for l in open('region_alns.sam'):
    f = l.split('\t'); fl = int(f[1])
    if fl & 0x100: continue
    q, rn, mq, cig = f[0], f[2], int(f[4]), f[5]
    if mq >= 20: parts[q][rn] = max(parts[q].get(rn, 0), span(cig))
    for t in f[11:]:
        if t.startswith('SA:Z:'):
            for sa in t[5:].strip().split(';'):
                if not sa: continue
                c, pos, st, sc, smq, nm = sa.split(',')
                if int(smq) >= 20: parts[q][c] = max(parts[q].get(c, 0), span(sc))
links = collections.Counter(); ctrl = collections.Counter(); out = []
for q, d in parts.items():
    big = {c for c, s in d.items() if s >= 1000}
    if 'LIQ_contig_19299' in big:
        for c in big - {'LIQ_contig_19299'}:
            links[c] += 1; out.append((q, 'LIQ_contig_19299', c, d['LIQ_contig_19299'], d[c], 'MAG_022' if c in mag else 'other'))
    m = sorted(big & mag)
    for i in range(len(m)):
        for j in range(i + 1, len(m)): ctrl[(m[i], m[j])] += 1
with open('../B5_linking_reads.tsv', 'w') as o:
    o.write('read_id\tcontig_a\tcontig_b\tspan_a_bp\tspan_b_bp\tcontig_b_in\n')
    for r in out: o.write('\t'.join(map(str, r)) + '\n')
with open('../B5_links_summary.tsv', 'w') as o:
    o.write('type\tcontig_a\tcontig_b\tn_reads\n')
    for c, n in links.most_common(): o.write(f'19299_link\tLIQ_contig_19299\t{c}\t{n}\n')
    for (a, b), n in ctrl.most_common(): o.write(f'positive_control_MAG022_MAG022\t{a}\t{b}\t{n}\n')
n_mag = sum(n for c, n in links.items() if c in mag)
open('link_stats.txt', 'w').write(f'{n_mag}\t{sum(links.values())}\t{len(ctrl)}\t{sum(ctrl.values())}\n')
print('19299->MAG022 reads', n_mag, '; 19299->any', sum(links.values()), '; control MAG022 pairs', len(ctrl), 'reads', sum(ctrl.values()))
