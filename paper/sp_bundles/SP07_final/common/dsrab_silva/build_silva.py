"""B1 v3: reference trees from dsrAB@SILVA (Jun 2023 release, GTDB r214 taxonomy; 1,319 dsrA / 1,280 dsrB incl. uncultured family-level lineages 1-13).
Outputs <gene>/ for dsrA, dsrB and dsrAB (concatenated DsrA+DsrB from adjacent genes of the same contig).
- translation: best of 6 frames (fewest stops), table 11; sequences with >1 internal stop dropped
- taxonomy: GTDB path; 'Unclassified_*' ranks end the path
- function: Oxidative (reverse-acting rDsrAB) = classes Gammaproteobacteria, Alphaproteobacteria, Chlorobia, Magnetococcia; all others Reductive
  (checked afterwards by monophyly in the tree, reported in function_check.txt)
- MMseqs2 clustering at 90% identity / 80% coverage to keep IQ-TREE tractable"""
import os, gzip, re, subprocess, collections
from Bio.Seq import Seq
H = os.path.dirname(os.path.abspath(__file__)); S = os.path.dirname(os.path.dirname(H))
CE = '${CONDA_ENVS_USER}'
OXI = {'Gammaproteobacteria', 'Alphaproteobacteria', 'Chlorobia', 'Magnetococcia'}
def fa(p):
    d, k = {}, None
    for l in gzip.open(p, 'rt'):
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = []
        else: d[k].append(l.strip())
    return {k: ''.join(v).upper() for k, v in d.items()}
def tr(nt):
    best = None
    for s in (Seq(nt), Seq(nt).reverse_complement()):
        for f in range(3):
            sub = s[f:]; sub = sub[:len(sub) // 3 * 3]; p = str(sub.translate(table=11)).rstrip('*')
            if best is None or p.count('*') < best.count('*'): best = p
    return best if best.count('*') <= 1 else None
def tax(p):
    t = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1].rstrip(';').split(';') for l in gzip.open(p, 'rt') if '\t' in l}
    out = {}
    for k, v in t.items():
        path = []
        for x in v[1:6]:   # phylum..genus (v[0] = domain)
            if x.startswith('Unclassified') or not x: break
            path.append(x)
        fn = 'Oxidative' if len(v) > 2 and v[2] in OXI else 'Reductive'
        out[k] = [fn] + path
    return out
log = collections.Counter(); prot, T = {}, {}
for g in ('dsrA', 'dsrB'):
    nt = fa(f'{H}/{g}.fasta.gz'); T[g] = tax(f'{H}/Reference_{g}_taxonomy_file_GTDB_r214.txt.gz'); prot[g] = {}
    for k, s in nt.items():
        p = tr(s)
        if p is None: log[f'{g}: dropped, >1 internal stop'] += 1; continue
        if k not in T[g]: log[f'{g}: dropped, no taxonomy'] += 1; continue
        prot[g][k] = p.replace('*', 'X')
    log[f'{g}: translated'] += len(prot[g])
# concatenated pairs: same genome + contig, gene indices adjacent
def key(k):
    m = re.match(r'(GC[AF]_\d+)_(.+)_(\d+)$', k); return (m.group(1), m.group(2), int(m.group(3))) if m else None
bidx = collections.defaultdict(dict)
for k in prot['dsrB']:
    x = key(k)
    if x: bidx[(x[0], x[1])][x[2]] = k
prot['dsrAB'] = {}; T['dsrAB'] = {}
for k in prot['dsrA']:
    x = key(k)
    if not x: continue
    cand = [bidx[(x[0], x[1])].get(x[2] + d) for d in (1, -1)]; cand = [c for c in cand if c]
    if len(cand) == 1:
        prot['dsrAB'][k] = prot['dsrA'][k] + prot['dsrB'][cand[0]]; T['dsrAB'][k] = T['dsrA'][k]
log['dsrAB: adjacent dsrA+dsrB pairs'] = len(prot['dsrAB'])
for g in ('dsrA', 'dsrB', 'dsrAB'):
    o = f'{H}/{g}'; os.makedirs(o, exist_ok=True)
    with open(f'{o}/all.faa', 'w') as f:
        for k, p in prot[g].items(): f.write(f'>{k}\n{p}\n')
    subprocess.run([f'{CE}/env_mmseqs/bin/mmseqs', 'easy-cluster', f'{o}/all.faa', f'{o}/clu', f'{o}/mmtmp', '--min-seq-id', '0.9', '-c', '0.8',
                    '--threads', '8'], check=True, capture_output=True)
    reps = sorted({l.split('\t')[0] for l in open(f'{o}/clu_cluster.tsv')})
    with open(f'{o}/{g}_in.faa', 'w') as f, open(f'{o}/taxon.tsv', 'w') as t:
        for k in reps: f.write(f'>{k}\n{prot[g][k]}\n'); t.write(f'{k}\t{";".join(T[g][k])}\n')
    log[f'{g}: references after 90% clustering'] = len(reps)
# DB for concatenated read segments: separate DsrA / DsrB representatives
with open(f'{H}/dsrAB/ab_sep.faa', 'w') as f:
    for pre, g in (('A__', 'dsrA'), ('B__', 'dsrB')):
        for l in open(f'{H}/{g}/{g}_in.faa'): f.write(l.replace('>', '>' + pre) if l.startswith('>') else l)
open(f'{H}/build_log.txt', 'w').write('\n'.join(f'{k}\t{v}' for k, v in sorted(log.items())) + '\n'); print(dict(log))
