"""SP07 B1/B2: place sediment marker reads in the SP05 trees.
usage: place.py <gene> <part_dir>
blastx (frameshift) -> aligned protein segment (>= 60 aa) -> mafft --addfragments --keeplength -> epa-ng -> gappa assign.
QC: leave-one-out re-placement of every SP05 reference (pruned tree + MSA)."""
import sys, os, re, subprocess, json, io, datetime, csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio import Phylo, SeqIO

GENE, PART = sys.argv[1], sys.argv[2]
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WS2 = os.path.dirname(os.path.dirname(S))
SP05 = f'{WS2}/99_thesis_audit/SP05_dsr_mcr'
EB = f'{S}/env/.pixi/envs/default/bin'
W = f'{PART}/work/{GENE}'; os.makedirs(W, exist_ok=True)
LOG = f'{PART}/COMMANDS.log'
CFG = {'dsrAB': dict(model='LG+G4+F', qc='family', ranks=['phylum', 'class', 'order', 'family', 'genus']),
       'dsrA': dict(model='LG+G4+F', qc='family', ranks=['phylum', 'class', 'order', 'family', 'genus']),
       'dsrB': dict(model='LG+G4+F', qc='family', ranks=['phylum', 'class', 'order', 'family', 'genus']),
       'mcrA': dict(model='LG+I+G4', qc='order', ranks=['order', 'family', 'genus'])}[GENE]
REFDIR = os.environ.get('REFDIR')  # v2: larger reference trees built in common/reftrees/<gene>
TREE, TRIM, FAA = ((f'{REFDIR}/{GENE}.treefile', f'{REFDIR}/{GENE}.trim', f'{REFDIR}/{GENE}_in.faa') if REFDIR else
                   (f'{SP05}/trees/{GENE}.treefile', f'{SP05}/trees/{GENE}.trim', f'{SP05}/trees/{GENE}_in.faa'))
MINSEG = 60

def run(cmd, **kw):
    with open(LOG, 'a') as L: L.write(f'CMD [{datetime.datetime.now().isoformat(timespec="seconds")}] {" ".join(cmd)}\n')
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    with open(LOG, 'a') as L: L.write(f'EXIT [{datetime.datetime.now().isoformat(timespec="seconds")}] rc={r.returncode} :: {cmd[0]}\n')
    if r.returncode: raise RuntimeError(f'{cmd[0]} failed: {r.stderr[-2000:]}')
    return r

# ---------------------------------------------------------------- taxonomy (NCBI taxdump from SP05)
par, rank, name, genus_ix = {}, {}, {}, defaultdict(list)
for l in open(f'{SP05}/lca/taxdump/nodes.dmp'):
    f = l.split('\t|\t'); par[f[0]] = f[1]; rank[f[0]] = f[2]
for l in open(f'{SP05}/lca/taxdump/names.dmp'):
    f = l.split('\t|\t')
    if f[3].startswith('scientific name'): name[f[0]] = f[1]
for t, n in name.items():
    if rank.get(t) == 'genus': genus_ix[n].append(t)
def lineage(t):
    out = {}
    while t in par and t != '1':
        out[rank[t]] = name.get(t, '?'); t = par[t]
    return out
def genus_lin(g):
    for cand in (g, f'Candidatus {g}', g.replace('Candidatus ', '')):
        for t in genus_ix.get(cand, []):
            L = lineage(t)
            if L.get('domain', L.get('superkingdom')) in ('Bacteria', 'Archaea'): return L
    return None
def path(func, L):
    p = [func]
    for r in CFG['ranks']:
        if not L or r not in L: break
        p.append(L[r])
    return p

leaves = [c.name for c in Phylo.read(TREE, 'newick').get_terminals()]
magtax = {r['mag_id']: r for r in csv.DictReader(open(f'{WS2}/07_mags/MAG_TAXONOMY.tsv'), delimiter='\t')}
org = {}
if GENE == 'mcrA':
    for fn in ('mcr', 'acr'):
        for l in open(f'{WS2}/15_marker_trees/mcrA/refs/{fn}.faa'):
            if l.startswith('>'):
                m = re.search(r'\[([^\]]+)\]\s*$', l); org[l[1:].split()[0]] = m.group(1) if m else ''
GTDB2NCBI_ORDER = {'Methanotrichales': 'Methanosarcinales'}
tax = {}
for lf in ([] if REFDIR else leaves):   # REFDIR trees carry their own taxon.tsv
    p = lf.split('|')
    if GENE.startswith('dsr'):
        if p[0] == 'QUERY':
            c = p[1].rsplit('_', 1)[0]
            if c == 'LIQ_contig_12880': L = genus_lin('Desulfosporosinus')
            else:  # LIQ_contig_19299: SP05 = Desulfotomaculia clade, genus unresolved -> NCBI lineage to order only
                L = {k: v for k, v in (genus_lin('Desulforamulus') or {}).items() if k in ('phylum', 'class', 'order')}
            tax[lf] = path('Reductive', L)
        else:
            toks = [t for t in p[1].split('_') if t and t not in ('Candidatus', 'uncultured')]
            tax[lf] = path('Reductive' if p[0] == 'REF_RED' else 'Oxidative', genus_lin(toks[0]) if toks else None)
    else:
        if p[0] == 'QUERY':
            m = magtax.get(p[1], {}); L = genus_lin(m.get('genus', '') or '') if m.get('genus') else None
            if not L and m.get('order'): L = {'order': GTDB2NCBI_ORDER.get(m['order'], m['order'])}
            tax[lf] = path('MCR', L)
        else:
            o = org.get(p[1], ''); toks = [t for t in o.split() if t != 'Candidatus']
            tax[lf] = path(p[0], genus_lin(toks[0]) if toks else None)
if REFDIR:
    tax = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1].split(';') for l in open(f'{REFDIR}/taxon.tsv')}
with open(f'{W}/taxon.tsv', 'w') as o:
    for lf in leaves: o.write(f'{lf}\t{";".join(tax[lf])}\n')

# ---------------------------------------------------------------- tree without support labels; reference MSA
def clean_tree(t):
    for c in t.find_clades():
        c.confidence = None
        if not c.is_terminal(): c.name = None
    return t
def unroot(t):  # epa-ng wants an unrooted (trifurcating-root) tree
    r = t.root
    if len(r.clades) == 2:
        a, b = r.clades
        x, y = (a, b) if not a.is_terminal() else (b, a)
        if not x.is_terminal():
            for k in x.clades: k.branch_length = (k.branch_length or 0)
            y.branch_length = (y.branch_length or 0) + (x.branch_length or 0)
            r.clades = list(x.clades) + [y]
    return t
def write_tree(t, p):
    s = io.StringIO(); Phylo.write(t, s, 'newick'); txt = s.getvalue().strip()
    open(p, 'w').write(txt + ('\n' if txt.endswith(';') else ';\n'))
T0 = unroot(clean_tree(Phylo.read(TREE, 'newick'))); write_tree(T0, f'{W}/ref.nwk')
msa = {r.id.replace('+', '_'): str(r.seq) for r in SeqIO.parse(TRIM, 'fasta')}  # IQ-TREE writes '+' as '_' in leaf names
assert set(msa) == set(leaves), 'trim/tree leaf mismatch'
with open(f'{W}/ref.aln', 'w') as o:
    for k in leaves: o.write(f'>{k}\n{msa[k]}\n')

def epa(tree, ref, qry, out, model):
    os.makedirs(out, exist_ok=True)
    try:
        run([f'{EB}/epa-ng', '--tree', tree, '--ref-msa', ref, '--query', qry, '--model', model, '--outdir', out, '--redo', '--threads', '4'])
        return model
    except RuntimeError:
        fb = 'LG+G4+F'
        run([f'{EB}/epa-ng', '--tree', tree, '--ref-msa', ref, '--query', qry, '--model', fb, '--outdir', out, '--redo', '--threads', '4'])
        return fb
def assign(jp, taxf, out):
    os.makedirs(out, exist_ok=True)
    run([f'{EB}/gappa', 'examine', 'assign', '--jplace-path', jp, '--taxon-file', taxf, '--per-query-results', '--out-dir', out, '--allow-file-overwriting'])
    res = defaultdict(list)
    for r in csv.DictReader(open(f'{out}/per_query.tsv'), delimiter='\t'):
        res[r['name']].append((float(r['aLWR']), r['taxopath'].split(';')))
    return res
def top_at(rows, depth):
    acc = defaultdict(float)
    for alwr, tp in rows:
        if len(tp) > depth: acc[tuple(tp[:depth + 1])] = max(acc[tuple(tp[:depth + 1])], alwr)
    if not acc: return None, 0.0
    k = max(acc, key=acc.get); return k[-1], acc[k]

# ---------------------------------------------------------------- QC: leave-one-out re-placement of references
qcd = depth_qc = 1 + CFG['ranks'].index(CFG['qc'])
def loo(lf):
    d = f'{W}/loo/{re.sub(r"[^A-Za-z0-9_.]", "_", lf)}'; os.makedirs(d, exist_ok=True)
    t = clean_tree(Phylo.read(TREE, 'newick')); t.prune(next(c for c in t.get_terminals() if c.name == lf)); write_tree(unroot(t), f'{d}/t.nwk')
    with open(f'{d}/ref.aln', 'w') as o:
        for k in leaves:
            if k != lf: o.write(f'>{k}\n{msa[k]}\n')
    open(f'{d}/q.aln', 'w').write(f'>{lf}\n{msa[lf]}\n')
    with open(f'{d}/taxon.tsv', 'w') as o:
        for k in leaves:
            if k != lf: o.write(f'{k}\t{";".join(tax[k])}\n')
    epa(f'{d}/t.nwk', f'{d}/ref.aln', f'{d}/q.aln', d, CFG['model'])
    rows = assign(f'{d}/epa_result.jplace', f'{d}/taxon.tsv', f'{d}/g').get(lf, [])
    exp_f, exp_q = tax[lf][0], (tax[lf][qcd] if len(tax[lf]) > qcd else None)
    got_f, w_f = top_at(rows, 0); got_q, w_q = top_at(rows, qcd)
    od = 1 + CFG['ranks'].index('order'); exp_o = tax[lf][od] if len(tax[lf]) > od else None; got_o, w_o = top_at(rows, od)
    ev_o = exp_o is not None and any(len(tax[k]) > od and tax[k][od] == exp_o for k in leaves if k != lf)
    evaluable = exp_q is not None and any(len(tax[k]) > qcd and tax[k][qcd] == exp_q for k in leaves if k != lf)
    return dict(leaf=lf, expected_function=exp_f, got_function=got_f, function_aLWR=round(w_f, 3), function_ok=got_f == exp_f,
                expected_qcrank=exp_q or 'none', got_qcrank=got_q or 'none', qcrank_aLWR=round(w_q, 3),
                evaluable=evaluable, qcrank_ok=(got_q == exp_q) if evaluable else None,
                expected_order=exp_o or 'none', got_order=got_o or 'none', order_aLWR=round(w_o, 3), order_evaluable=ev_o, order_ok=(got_o == exp_o) if ev_o else None)
refs = [lf for lf in leaves if not lf.startswith('QUERY')]
with ThreadPoolExecutor(6) as ex: qc = list(ex.map(loo, refs))
with open(f'{PART}/QC_{GENE}_leave_one_out.tsv', 'w') as o:
    k = list(qc[0]); o.write('\t'.join(k) + '\n')
    for r in qc: o.write('\t'.join(str(r[x]) for x in k) + '\n')
ev = [r for r in qc if r['evaluable']]
fq = sum(r['function_ok'] for r in qc) / len(qc); rq = (sum(r['qcrank_ok'] for r in ev) / len(ev)) if ev else 0
qc_line = (f'{GENE}: leave-one-out re-placement of {len(qc)} references -> function clade recovered {fq:.1%}; '
           f'{CFG["qc"]} recovered {rq:.1%} of {len(ev)} evaluable (taxon present elsewhere in tree). '
           f'Threshold 95% -> {"PASS" if fq >= 0.95 and rq >= 0.95 else "FAIL"}')
if CFG['qc'] != 'order':
    evo = [r for r in qc if r['order_evaluable']]; ro = sum(r['order_ok'] for r in evo) / len(evo) if evo else 0
    qc_line += (f'\n{GENE} order-level fallback: order recovered {ro:.1%} of {len(evo)} evaluable -> '
                f'{"PASS" if fq >= 0.95 and ro >= 0.95 else "FAIL"}')
print(qc_line)

# ---------------------------------------------------------------- reads: blastx -> segment
ids = ([l.strip() for l in open(f'{REFDIR}/reads.txt')] if GENE == 'dsrAB' else
       [l.split('\t')[1] for l in open(f'{S}/common/reads/read_ids.tsv') if l.startswith(GENE + '\t')])
recs = {r.id: r for r in SeqIO.parse(f'{S}/common/reads/reads.fa', 'fasta')}
missing = [i for i in ids if i not in recs]
if missing: raise SystemExit(f'QC FAIL: {len(missing)} {GENE} reads not found in extracted reads')
with open(f'{W}/reads.fa', 'w') as o:
    for i in ids: o.write(f'>{i}\n{recs[i].seq}\n')
run([f'{EB}/diamond', 'makedb', '--in', f'{REFDIR}/ab_sep.faa' if GENE == 'dsrAB' else FAA, '-d', f'{W}/refs', '--quiet'])
run([f'{EB}/diamond', 'blastx', '-d', f'{W}/refs.dmnd', '-q', f'{W}/reads.fa', '--frameshift', '15', '--more-sensitive', '--evalue', '1e-10',
     '-k', '0', '--threads', '8', '--quiet', '-o', f'{W}/blastx.tsv',
     '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'qstart', 'qend', 'qseq_translated'])
hits = defaultdict(list)
for l in open(f'{W}/blastx.tsv'):
    f = l.rstrip('\n').split('\t'); hits[f[0]].append(f)
seg = {}
for q, hs in hits.items():
    if GENE == 'dsrAB':   # concatenated: best DsrA segment + best DsrB segment of the same read, each >= MINSEG
        parts = []
        for pre in ('A__', 'B__'):
            h = [f for f in hs if f[1].startswith(pre)]
            if h: b = max(h, key=lambda f: float(f[5])); parts.append((b, re.sub(r'[^A-Za-z]', '', b[8]).upper()))
        if len(parts) == 2 and all(len(x[1]) >= MINSEG for x in parts):
            b = list(parts[0][0]); b[1] = parts[0][0][1] + '+' + parts[1][0][1]; seg[q] = (b, parts[0][1] + parts[1][1])
        continue
    b = max(hs, key=lambda f: float(f[5])); s = re.sub(r'[^A-Za-z]', '', b[8]).upper()
    seg[q] = (b, s)
with open(f'{W}/segments.faa', 'w') as o:
    for q, (b, s) in seg.items():
        if len(s) >= MINSEG: o.write(f'>{q}\n{s}\n')
placed = [q for q, (b, s) in seg.items() if len(s) >= MINSEG]
res, used_model = {}, CFG['model']
if placed:
    out = run([f'{EB}/mafft', '--addfragments', f'{W}/segments.faa', '--keeplength', '--thread', '8', f'{W}/ref.aln'])
    al = {r.id: str(r.seq) for r in SeqIO.parse(io.StringIO(out.stdout), 'fasta')}
    with open(f'{W}/query.aln', 'w') as o:
        for q in placed: o.write(f'>{q}\n{al[q]}\n')
    used_model = epa(f'{W}/ref.nwk', f'{W}/ref.aln', f'{W}/query.aln', f'{W}/epa', CFG['model'])
    res = assign(f'{W}/epa/epa_result.jplace', f'{W}/taxon.tsv', f'{W}/gappa')
    jp = json.load(open(f'{W}/epa/epa_result.jplace'))
    tr = Phylo.read(io.StringIO(re.sub(r'\{(\d+)\}', r'[\1]', jp['tree'])), 'newick')
    edge_leaves = {int(c.comment): {t.name for t in c.get_terminals()} for c in tr.find_clades() if c.comment}
    best_edge = {pl['n'][0] if 'n' in pl else pl['nm'][0][0]: max(pl['p'], key=lambda x: x[2]) for pl in jp['placements']}
else:
    edge_leaves, best_edge = {}, {}

# ---------------------------------------------------------------- UniProt LCA (shared run in common/uniprot)
lca = {}
up = f'{S}/common/uniprot/read_LCA.tsv'
if os.path.exists(up):
    for r in csv.DictReader(open(up), delimiter='\t'): lca[r['read_id']] = r

# ---------------------------------------------------------------- per-read table
mag11 = [lf for lf in leaves if 'ETO_LIQ_MAG_011' in lf]
cols = ['gene', 'read_id', 'read_len', 'best_ref', 'best_ref_pident', 'aln_aa', 'segment_aa', 'placed', 'function', 'function_aLWR']
cols += [f'{r}_top' for r in CFG['ranks']] + [f'{r}_aLWR' for r in CFG['ranks']]
cols += [f'assigned_{CFG["qc"]}', 'assigned_genus', 'best_placement_edge_leaves', 'uniprot_top_hit', 'uniprot_top_pident', 'uniprot_LCA_rank', 'uniprot_LCA', 'uniprot_LCA_lineage']
if GENE == 'mcrA': cols += ['falls_with_MAG_011', 'identity_to_MAG_011']
rows = []
for q in ids:
    b, s = seg.get(q, (None, ''))
    r = dict(gene=GENE, read_id=q, read_len=len(recs[q].seq), best_ref=b[1] if b else 'no hit (1e-10)', best_ref_pident=b[2] if b else 'none',
             aln_aa=b[3] if b else 0, segment_aa=len(s), placed='yes' if q in res else ('no: segment < 60 aa' if b else 'no: no blastx hit'))
    rr = res.get(q, [])
    fn, fw = top_at(rr, 0); r['function'] = fn if fn and fw >= 0.5 else 'unassigned'; r['function_aLWR'] = round(fw, 3)
    for i, rk in enumerate(CFG['ranks'], 1):
        t, w = top_at(rr, i); r[f'{rk}_top'] = t or 'none'; r[f'{rk}_aLWR'] = round(w, 3)
    qc_ok = r[f'{CFG["qc"]}_top'] != 'none' and r[f'{CFG["qc"]}_aLWR'] >= 0.5
    r[f'assigned_{CFG["qc"]}'] = r[f'{CFG["qc"]}_top'] if qc_ok else f'unassigned at {CFG["qc"]} level'
    r['assigned_genus'] = r['genus_top'] if r['genus_top'] != 'none' and r['genus_aLWR'] >= 0.5 else 'unassigned'
    e = best_edge.get(q); el = edge_leaves.get(e[0], set()) if e else set()
    if e and len(el) > len(leaves) / 2: el = set(leaves) - el
    r['best_placement_edge_leaves'] = (';'.join(sorted(el)) if len(el) <= 5 else f'{len(el)} leaves') if e else 'none'
    u = lca.get(q, {})
    for k in ('uniprot_top_hit', 'uniprot_top_pident', 'uniprot_LCA_rank', 'uniprot_LCA', 'uniprot_LCA_lineage'): r[k] = u.get(k, 'not run') or 'none'
    if GENE == 'mcrA':
        r['falls_with_MAG_011'] = 'yes' if (e and any(m in el for m in mag11) and len(el) <= 3) else 'no'
        m11 = [float(f[2]) for f in hits.get(q, []) if 'ETO_LIQ_MAG_011' in f[1]]
        r['identity_to_MAG_011'] = f'{max(m11):.1f}' if m11 else 'no hit'
    rows.append(r)
with open(f'{PART}/{GENE}_read_placements.tsv', 'w') as o:
    o.write('\t'.join(cols) + '\n')
    for r in rows: o.write('\t'.join(str(r[c]) for c in cols) + '\n')
open(f'{PART}/QC_{GENE}.txt', 'w').write(qc_line + f'\nplacement model used: {used_model}; reads {len(ids)}, blastx hits {len(seg)}, placed (segment >= {MINSEG} aa) {len(placed)}\n')
print(open(f'{PART}/QC_{GENE}.txt').read())
