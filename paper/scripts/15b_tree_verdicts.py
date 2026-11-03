#!/usr/bin/env python3
"""
15b_tree_verdicts.py -- assign each query sequence by PATRISTIC DISTANCE.

TWO EARLIER METHODS FAILED, FOR DIFFERENT REASONS
  1. Newick string adjacency (in 15_marker_trees.sh). Newick order is not
     phylogenetic proximity. It called ETO_LIQ_MAG_006 (Vallitalea, a
     bacterium) as cfbC/cfbD -- an archaeal F430 biosynthesis protein.
  2. Walking up to the nearest reference-containing clade (v1 of this script).
     IQ-TREE writes UNROOTED trees, so "up" heads toward an arbitrary
     trifurcation. Any query attaching near that point returned the whole tree
     as its clade -- which is why every mttB query came back
     "TMA_MttB:40, BETAINE_SF:40" and read as AMBIGUOUS.

THIS METHOD
  Patristic distance -- the summed branch lengths along the path between two
  leaves -- is independent of rooting. For each query, the K nearest reference
  leaves are found and their group composition reported. A query whose nearest
  neighbours are >=PURITY one group is assigned to it.

  Distances are computed by a single traversal of the tree treated as an
  undirected graph, so the arbitrary root plays no part.

  The nearest single reference and its distance are also reported: a query that
  is far from everything (long branch) is flagged, because "closest to MCR" is
  weak evidence if the closest MCR is still very distant.

No dependencies beyond the standard library.
USAGE  python3 scripts/15b_tree_verdicts.py 15_marker_trees [K] [PURITY]
"""
import sys, os, re, heapq
from collections import Counter, defaultdict

GROUPS = {
    'mcrA': {'good': 'MCR', 'bad': 'ACR',
             'note': 'MCR = true methanogen; ACR = anaerobic alkane oxidiser.'},
    'nifH': {'good': 'NIFH', 'bad': 'CFBCD',
             'note': 'NIFH = nitrogenase; CFBCD = coenzyme F430 biosynthesis '
                     '(archaeal; every methanogen needs it for MCR); '
                     'BCHL = protochlorophyllide reductase.'},
    'mttB': {'good': 'TMA_MttB', 'bad': 'BETAINE_SF',
             'note': 'TMA_MttB = methylamine methanogenesis; '
                     'BETAINE_SF = glycine betaine superfamily.'},
}


def parse_newick(text):
    """Return (edges, leaf_names). edges: adjacency dict id -> [(id, length)]."""
    text = text.strip().rstrip(';')
    adj = defaultdict(list)
    names = {}
    nid = [0]

    def new_node():
        nid[0] += 1
        return nid[0]

    root = new_node()
    stack = [root]
    cur = root
    token = ''

    def finish(node, tok):
        tok = tok.strip()
        if not tok:
            return 0.0
        if ':' in tok:
            lab, _, ln = tok.rpartition(':')
            try:
                length = float(ln)
            except ValueError:
                lab, length = tok, 0.0
        else:
            lab, length = tok, 0.0
        lab = lab.strip().strip("'\"")
        if lab and not re.fullmatch(r'\d+(\.\d+)?', lab):
            names[node] = lab
        return length

    for c in text:
        if c == '(':
            child = new_node()
            adj[cur].append((child, None))
            adj[child].append((cur, None))
            stack.append(cur)
            cur = child
            token = ''
        elif c == ',':
            ln = finish(cur, token); token = ''
            _set_len(adj, cur, stack[-1], ln)
            parent = stack[-1]
            sib = new_node()
            adj[parent].append((sib, None))
            adj[sib].append((parent, None))
            cur = sib
        elif c == ')':
            ln = finish(cur, token); token = ''
            _set_len(adj, cur, stack[-1], ln)
            cur = stack.pop()
        else:
            token += c
    finish(cur, token)

    leaves = {n: names[n] for n in names if len(adj[n]) <= 1}
    return adj, leaves, names


def _set_len(adj, a, b, ln):
    for i, (n, l) in enumerate(adj[a]):
        if n == b:
            adj[a][i] = (n, ln)
    for i, (n, l) in enumerate(adj[b]):
        if n == a:
            adj[b][i] = (n, ln)


def distances_from(adj, start):
    """Dijkstra over the tree (all lengths >= 0); returns node -> distance."""
    dist = {start: 0.0}
    pq = [(0.0, start)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, float('inf')):
            continue
        for v, ln in adj[u]:
            nd = d + (ln if ln else 0.0)
            if nd < dist.get(v, float('inf')):
                dist[v] = nd
                heapq.heappush(pq, (nd, v))
    return dist


def group_of(name):
    if name.startswith('QUERY|'):
        return 'QUERY'
    m = re.match(r'^([A-Za-z0-9_]+)\|', name)
    return m.group(1) if m else 'UNLABELLED'


def main(base, K=10, PURITY=0.80):
    for mk, cfg in GROUPS.items():
        tf = os.path.join(base, mk, f'{mk}.treefile')
        print()
        print('=' * 104)
        print(f'  {mk} — PATRISTIC-DISTANCE VERDICTS  (K={K} nearest references, '
              f'purity >={PURITY:.0%})')
        print('=' * 104)
        print(f'  {cfg["note"]}')
        if not os.path.exists(tf):
            print('  no tree found'); continue

        adj, leaves, _ = parse_newick(open(tf).read())
        queries = {n: nm for n, nm in leaves.items() if nm.startswith('QUERY|')}
        refs = {n: nm for n, nm in leaves.items() if not nm.startswith('QUERY|')}
        print(f'  tree: {len(leaves)} leaves — {len(queries)} queries, '
              f'references: {dict(Counter(group_of(v) for v in refs.values()))}')
        if not queries or not refs:
            continue

        print()
        print(f'  {"query":<46}{"K nearest":<30}{"closest ref":<26}{"dist":>7}  verdict')
        rows = []
        for qn, qname in sorted(queries.items(), key=lambda kv: kv[1]):
            dist = distances_from(adj, qn)
            near = sorted(((dist.get(rn, float('inf')), group_of(rname), rname)
                           for rn, rname in refs.items()))[:K]
            comp = Counter(g for _, g, _ in near)
            top, cnt = comp.most_common(1)[0]
            frac = cnt / len(near)
            d0, g0, n0 = near[0]

            if frac >= PURITY and top == cfg['good']:
                verdict = 'CONFIRMED'
            elif frac >= PURITY and top == cfg['bad']:
                verdict = 'REJECTED'
            elif frac >= PURITY:
                verdict = f'other:{top}'
            else:
                verdict = 'AMBIGUOUS'

            # The single NEAREST reference disagreeing with the K-nearest
            # majority is a genuine conflict: it usually means the query is
            # sister to one sequence of group A while sitting inside a cloud of
            # group B. Never let that pass as a clean call.
            if verdict in ('CONFIRMED', 'REJECTED') and g0 != top:
                verdict += f' !conflict(nearest={g0})'

            short = qname.replace('QUERY|', '')
            desc = ', '.join(f'{g}:{c}' for g, c in comp.most_common(3))
            print(f'  {short[:46]:<46}{desc[:30]:<30}{n0.split("|")[-1][:26]:<26}'
                  f'{d0:>7.3f}  {verdict}')
            rows.append((short, verdict, desc, f'{d0:.4f}', n0))

        out = os.path.join(base, mk, f'{mk}_verdicts.tsv')
        with open(out, 'w') as fh:
            fh.write('query\tverdict\tk_nearest_composition\tdist_to_closest_ref\tclosest_ref\n')
            for r in rows:
                fh.write('\t'.join(r) + '\n')

        print()
        print('  ' + '   '.join(f'{k}: {v}' for k, v in Counter(r[1] for r in rows).most_common()))

        # per-MAG roll-up (query names are MAGID|working|locus)
        mags = defaultdict(list)
        for short, verdict, _, _, _ in rows:
            mags[short.split('|')[0]].append(verdict)
        pos = sorted(m for m, v in mags.items() if 'CONFIRMED' in v)
        neg = sorted(m for m, v in mags.items() if 'CONFIRMED' not in v)
        print(f'  MAGs with >=1 CONFIRMED copy ({len(pos)}): {", ".join(pos) or "none"}')
        if neg:
            print(f'  MAGs with none ({len(neg)}): {", ".join(neg)}')
        print(f'  written: {out}')

    print()
    print('  Patristic distance is rooting-independent, so these verdicts do not')
    print('  depend on where IQ-TREE happened to place the arbitrary root. Check the')
    print('  "dist" column: a CONFIRMED call whose closest reference is still far')
    print('  away means a divergent homologue, not a textbook match.')


if __name__ == '__main__':
    base = sys.argv[1] if len(sys.argv) > 1 else '15_marker_trees'
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    P = float(sys.argv[3]) if len(sys.argv) > 3 else 0.80
    main(base, K, P)
