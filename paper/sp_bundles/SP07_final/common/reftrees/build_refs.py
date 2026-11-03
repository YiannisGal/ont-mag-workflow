"""Larger reference sets for B1/B2 (v2). usage: build_refs.py <gene>
Pool = SP05 tree references + step-58 dsrAB references (dsr) + UniProt top-10 hits of every sediment read of that gene.
Filters (fixed before any v2 result was seen):
  1. KOfam profile score >= curated threshold (K11180 DsrA, K11181 DsrB, K00399 McrA)
  2. length >= 70% of the median length of the SP05 references
  3. taxonomy: dsr -> NCBI lineage resolved to family (Bacillota SRB genera grouped GTDB-style, see FAM_OVERRIDE);
               mcrA -> archaeal lineage resolved to order (removes bacterial-host 'McrA' entries)
  4. MMseqs2 clustering at 95% identity, 80% coverage (representatives kept; SP05/step-58 references always kept)
ACR references are exempt from the K00399 threshold (alkane-oxidising homologs score below it by design).
Function labels: dsr oxidative = SP05 REF_OXI or a genus with reverse-acting rDsrAB (SOB list); mcrA ACR = SP05 ACR or ACR genus list."""
import sys, os, re, subprocess, collections, datetime
GENE = sys.argv[1]
H = os.path.dirname(os.path.abspath(__file__)); S = os.path.dirname(os.path.dirname(H)); WS2 = os.path.dirname(os.path.dirname(S))
SP05 = f'{WS2}/99_thesis_audit/SP05_dsr_mcr'; EB = f'{S}/env/.pixi/envs/default/bin'; CE = '${CONDA_ENVS_USER}'
O = f'{H}/{GENE}'; os.makedirs(O, exist_ok=True)
KO = {'dsrA': 'K11180', 'dsrB': 'K11181', 'mcrA': 'K00399'}[GENE]
def run(cmd, **kw):
    open(f'{S}/common/COMMANDS.log', 'a').write(f'CMD [{datetime.datetime.now().isoformat(timespec="seconds")}] {" ".join(cmd)}\n')
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    open(f'{S}/common/COMMANDS.log', 'a').write(f'EXIT rc={r.returncode} :: {cmd[0]}\n')
    if r.returncode: raise RuntimeError(r.stderr[-2000:])
    return r
def fa(p):
    d, k = {}, None
    for l in open(p):
        if l.startswith('>'): k = l[1:].rstrip('\n'); d[k] = []
        else: d[k].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}
# ---------- taxonomy
par, rank, name, gix = {}, {}, {}, collections.defaultdict(list)
for l in open(f'{SP05}/lca/taxdump/nodes.dmp'):
    f = l.split('\t|\t'); par[f[0]] = f[1]; rank[f[0]] = f[2]
for l in open(f'{SP05}/lca/taxdump/names.dmp'):
    f = l.split('\t|\t')
    if f[3].startswith('scientific name'): name[f[0]] = f[1]
for t, n in name.items():
    if rank.get(t) == 'genus': gix[n].append(t)
def lin(t):
    out = {}
    while t in par and t != '1': out[rank[t]] = name.get(t, '?'); t = par[t]
    return out
def glin(g):
    for c in (g, f'Candidatus {g}'):
        for t in gix.get(c, []):
            L = lin(t)
            if L.get('domain', L.get('superkingdom')) in ('Bacteria', 'Archaea'): return L
    return {}
DESULFOTOMACULIA = {'Desulfotomaculum', 'Desulforamulus', 'Desulfofalx', 'Pelotomaculum', 'Desulfofundulus', 'Desulfallas', 'Desulfohalotomaculum',
                    'Desulfofarcimen', 'Desulfurispora', 'Desulfocucumis', 'Desulfoscipio', 'Desulfotruncus', 'Desulfovirgula', 'Desulfosoma',
                    'Desulfotomaculum_A', 'Sporotomaculum', 'Thermincola', 'Desulfitibacter', 'Desulfonispora'}
DESULFITOBACTERIACEAE = {'Desulfosporosinus', 'Desulfitobacterium', 'Dehalobacter', 'Syntrophobotulus', 'Desulfosporosinus_A'}
def fam_override(L):
    g = L.get('genus', '')
    if g in DESULFOTOMACULIA: L = dict(L, family='Desulfotomaculia (GTDB class; NCBI Peptococcaceae/Desulfotomaculaceae)')
    if g in DESULFITOBACTERIACEAE: L = dict(L, family='Desulfitobacteriaceae')
    return L
SOB = {'Allochromatium', 'Chlorobium', 'Chlorobaculum', 'Prosthecochloris', 'Pelodictyon', 'Thiocapsa', 'Thiorhodococcus', 'Thiocystis', 'Thiobacillus',
       'Sideroxydans', 'Beggiatoa', 'Thioploca', 'Thiothrix', 'Sulfuricella', 'Thiohalocapsa', 'Halochromatium', 'Marichromatium', 'Thioalkalivibrio',
       'Ruthia', 'Vesicomyosocius', 'Thiomargarita', 'Thermochromatium', 'Lamprocystis', 'Thiodictyon', 'Thiolamprovum', 'Ectothiorhodospira',
       'Thioflavicoccus', 'Halorhodospira', 'Thiohalobacter', 'Thiorhodovibrio', 'Thiolapillus', 'Sedimenticola', 'Thiogranum', 'Thiohalophilus',
       'Magnetospirillum', 'Magnetovibrio', 'Thiodiazotropha', 'Thioglobus', 'Thiomicrorhabdus', 'Thiofilum', 'Achromatium', 'Thiobaca', 'Rhodovulum'}
ACR = {'Syntrophoarchaeum', 'Methanoliparum', 'Argoarchaeum', 'Ethanoperedens', 'Alkanophaga', 'Polytropus', 'Helarchaeum', 'Hadesarchaea'}
RK = ['phylum', 'class', 'order', 'family', 'genus'] if GENE != 'mcrA' else ['order', 'family', 'genus']
# ---------- pool
pool = {}   # leaf -> (seq, function, lineage dict, source)
sp05 = fa(f'{SP05}/trees/{GENE}_in.faa')
magtax = {}
for l in open(f'{WS2}/07_mags/MAG_TAXONOMY.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0].startswith('ETO_'): magtax[f[0]] = f
hdr = open(f'{WS2}/07_mags/MAG_TAXONOMY.tsv').readline().rstrip('\n').split('\t')
gi = hdr.index('genus'); oi = hdr.index('order')
org = {}
if GENE == 'mcrA':
    for fn in ('mcr', 'acr'):
        for l in open(f'{WS2}/15_marker_trees/mcrA/refs/{fn}.faa'):
            if l.startswith('>'):
                m = re.search(r'\[([^\]]+)\]\s*$', l); org[l[1:].split()[0]] = m.group(1) if m else ''
for k, s in sp05.items():
    k = k.split()[0].replace('+', '_'); p = k.split('|')
    if GENE != 'mcrA':
        if p[0] == 'QUERY':
            L = fam_override(glin('Desulfosporosinus')) if 'LIQ_contig_12880' in k else {**{x: y for x, y in glin('Desulforamulus').items() if x in ('phylum', 'class', 'order')}, 'family': 'Desulfotomaculia (GTDB class; NCBI Peptococcaceae/Desulfotomaculaceae)'}
            fn = 'Reductive'
        else:
            toks = [t for t in p[1].split('_') if t and t not in ('Candidatus', 'uncultured')]
            L = fam_override(glin(toks[0])) if toks else {}; fn = 'Oxidative' if p[0] == 'REF_OXI' else 'Reductive'
    else:
        if p[0] == 'QUERY':
            m = magtax.get(p[1]); L = glin(m[gi]) if m and m[gi] else {}
            if not L and m: L = {'order': {'Methanotrichales': 'Methanosarcinales'}.get(m[oi], m[oi])}
            fn = 'MCR'
        else:
            toks = [t for t in org.get(p[1], '').split() if t != 'Candidatus']; L = glin(toks[0]) if toks else {}; fn = p[0]
    pool[k] = (s, fn, L, 'SP05')
if GENE != 'mcrA':
    for k, s in fa(f'{WS2}/58_dsr_finish/refs.faa').items():
        p = k.split()[0].split('|'); toks = [t for t in p[1].split('_') if t and t not in ('Candidatus', 'uncultured')]
        g = toks[0] if toks else ''
        pool[f'S58|{p[1]}|{p[2]}'] = (s, 'Oxidative' if g in SOB else 'Reductive', fam_override(glin(g)) if g else {}, 'step58')
want = set()
reads = {l.split('\t')[1] for l in open(f'{S}/common/reads/read_ids.tsv') if l.startswith(GENE + '\t')}
for db in ('sprot', 'trembl'):
    for l in open(f'{S}/common/uniprot/hits_{db}.tsv'):
        f = l.split('\t')
        if f[0] in reads: want.add(f[1])
for k, s in fa(f'{H}/candidates.faa').items():
    sid = k.split()[0]
    if sid not in want: continue
    m = re.search(r'OX=(\d+)', k); L = lin(m.group(1)) if m else {}
    g = L.get('genus', ''); L = fam_override(L)
    fn = ('Oxidative' if g in SOB else 'Reductive') if GENE != 'mcrA' else ('ACR' if g in ACR else 'MCR')
    pool[f'UP|{sid.split("|")[1]}|{re.sub(r"[^A-Za-z0-9]+", "_", L.get("species", g or "unknown"))}'] = (s, fn, L, 'UniProt')
# ---------- 1. KOfam threshold
with open(f'{O}/pool.faa', 'w') as o:
    for i, (k, v) in enumerate(pool.items()): o.write(f'>s{i}\n{v[0]}\n')
idx = {f's{i}': k for i, k in enumerate(pool)}
thr = float(next(l.split('\t')[1] for l in open(f'{WS2}/11_kofam/ko_list') if l.startswith(KO + '\t')))
run([f'{CE}/SemiBin/bin/hmmsearch', '--cpu', '8', '--tblout', f'{O}/kofam.tbl', '-o', '/dev/null', f'{WS2}/34_unbinned/profiles/{KO}.hmm', f'{O}/pool.faa'])
sc = collections.defaultdict(float)
for l in open(f'{O}/kofam.tbl'):
    if not l.startswith('#'): f = l.split(); sc[f[0]] = max(sc[f[0]], float(f[5]))
med = sorted(len(v[0]) for k, v in pool.items() if v[3] == 'SP05' and not k.startswith('QUERY'))
minlen = 0.7 * med[len(med) // 2]
log = collections.Counter(); keep = {}
for sid, k in idx.items():
    s, fn, L, src = pool[k]
    if k.startswith('QUERY'): keep[k] = pool[k]; log['kept: SP05 query leaf'] += 1; continue
    if sc[sid] < thr and fn != 'ACR': log[f'dropped: KOfam score < {thr} ({src})'] += 1; continue
    if len(s) < minlen: log[f'dropped: length < {minlen:.0f} aa ({src})'] += 1; continue
    if GENE == 'mcrA' and (L.get('domain', L.get('superkingdom')) != 'Archaea' or 'order' not in L): log[f'dropped: not archaeal or no order ({src})'] += 1; continue
    if GENE != 'mcrA' and 'family' not in L: log[f'dropped: no family-level lineage ({src})'] += 1; continue
    keep[k] = pool[k]; log[f'kept ({src})'] += 1
# ---------- 4. cluster UniProt additions at 95%
up = {k: v for k, v in keep.items() if v[3] == 'UniProt'}
with open(f'{O}/up.faa', 'w') as o:
    for k, v in up.items(): o.write(f'>{k}\n{v[0]}\n')
if up:
    run([f'{CE}/env_mmseqs/bin/mmseqs', 'easy-cluster', f'{O}/up.faa', f'{O}/clu', f'{O}/mmtmp', '--min-seq-id', '0.95', '-c', '0.8', '--threads', '8'])
    reps = {l.split('\t')[0] for l in open(f'{O}/clu_cluster.tsv')}
    for k in list(up):
        if k not in reps: del keep[k]; log['dropped: redundant at 95% (UniProt)'] += 1
# ---------- write
with open(f'{O}/{GENE}_in.faa', 'w') as o, open(f'{O}/taxon.tsv', 'w') as t, open(f'{O}/refs_table.tsv', 'w') as r:
    r.write('leaf\tsource\tfunction\ttaxopath\tlength\tkofam_score\n')
    inv = {v: k for k, v in idx.items()}
    for k, (s, fn, L, src) in keep.items():
        o.write(f'>{k}\n{s}\n'); path = [fn] + [L[x] for x in RK if x in L][:len(RK)]
        # keep only a contiguous prefix of ranks
        p2 = [fn]
        for x in RK:
            if x not in L: break
            p2.append(L[x])
        t.write(f'{k}\t{";".join(p2)}\n'); r.write(f'{k}\t{src}\t{fn}\t{";".join(p2)}\t{len(s)}\t{sc[inv[k]]:.1f}\n')
open(f'{O}/build_log.txt', 'w').write('\n'.join(f'{n}\t{k}' for k, n in sorted(log.items())) + f'\nfinal references: {len(keep)}\n')
print(GENE, dict(log), 'final', len(keep))
