"""SP01 step 4b -- classify each NifH/NifD query by patristic nearest labelled reference in the FastTree tree.
verdict uses the nearest reference; 'ambiguous' unless >=4 of the 5 nearest references agree with the nearest on nitrogenase vs not.
CfbC/CfbD, NflH/NflD (group IV) and Bch (group V) are NOT nitrogenases."""
import os, re
OUT = os.environ['OUT']; N = f'{OUT}/nif'
NITROGENASE = {'NIF_I': 'Nif group I', 'NIF_III': 'Nif group III', 'ANF': 'Anf (group II alt.)', 'VNF': 'Vnf (alt.)'}
NOT = {'CFB': 'CfbC/CfbD (F430 biosynthesis)', 'NFL': 'NflH/NflD (group IV)', 'BCH': 'Bch/Chl (group V)'}

def parse(s):
    """Newick -> (parent, blen, leafname) dicts keyed by node id."""
    s = s.strip().rstrip(';'); parent, blen, name = {}, {}, {}
    stack, cur, nid, i = [], None, 0, 0
    root = nid; parent[root] = None; blen[root] = 0.0; cur = root
    while i < len(s):
        c = s[i]
        if c == '(':
            nid += 1; parent[nid] = cur; blen[nid] = 0.0; stack.append(cur); cur = nid; i += 1
        elif c == ',':
            p = stack[-1]; nid += 1; parent[nid] = p; blen[nid] = 0.0; cur = nid; i += 1
        elif c == ')':
            cur = stack.pop(); i += 1
        else:
            m = re.match(r"([^:,();]*)(?::([-0-9.eE]+))?", s[i:])
            lab, bl = m.group(1), m.group(2)
            if lab and not re.fullmatch(r'[0-9.]+', lab): name[cur] = lab
            if bl: blen[cur] = float(bl)
            i += m.end()
    return parent, blen, name

def dist_all(parent, blen, name):
    kids = {}
    for n, p in parent.items():
        if p is not None: kids.setdefault(p, []).append(n)
    adj = {n: [] for n in parent}
    for n, p in parent.items():
        if p is not None: adj[n].append((p, blen[n])); adj[p].append((n, blen[n]))
    def from_leaf(src):
        d = {src: 0.0}; st = [src]
        while st:
            u = st.pop()
            for v, w in adj[u]:
                if v not in d: d[v] = d[u] + w; st.append(v)
        return d
    return from_leaf

rows = []
for sub, gene in (('H', 'nifH K02588'), ('D', 'nifD K02586')):
    p = f'{N}/nif{sub}.tree'
    parent, blen, name = parse(open(p).read())
    fl = dist_all(parent, blen, name)
    leaves = {v: k for k, v in name.items()}
    refs = [l for l in leaves if l.startswith('REF_')]
    for q in (l for l in leaves if l.startswith('Q__')):
        d = fl(leaves[q])
        near = sorted(refs, key=lambda r: d[leaves[r]])[:5]
        labs = [r[4:].split('__')[0] for r in near]
        top = labs[0]
        is_n = [l in NITROGENASE for l in labs]
        agree = sum(x == (top in NITROGENASE) for x in is_n)
        if agree >= 4:
            verdict = ('NITROGENASE: ' + NITROGENASE[top]) if top in NITROGENASE else ('NOT a nitrogenase: ' + NOT[top])
        else:
            verdict = f'AMBIGUOUS (only {agree}/5 nearest refs agree with nearest on nitrogenase vs not)'
        g, prot = q[3:].split('__', 1)
        rows.append((g, prot.replace('__', '|'), gene, top, near[0].split('__')[1], f'{d[leaves[near[0]]]:.3f}', ','.join(labs), verdict))
rows.sort()
with open(f'{OUT}/nif_tree_verdicts.tsv', 'w') as o:
    o.write('genome\tprotein\tmarker\tnearest_ref_group\tnearest_ref_acc\tpatristic_dist\t5_nearest_groups\tverdict\n')
    for r in rows: o.write('\t'.join(r) + '\n')
for r in rows: print('\t'.join(r))
print('queries classified:', len(rows))
