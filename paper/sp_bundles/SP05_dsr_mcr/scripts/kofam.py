"""Shared: KOfam calls from work/kofam_all.tbl with ko_list thresholds (score_type honoured)."""
import os
from collections import defaultdict
def load():
    ko = {}
    with open(os.environ['KOLIST']) as fh:
        fh.readline()
        for l in fh:
            f = l.rstrip('\n').split('\t')
            ko[f[0]] = (float(f[1]) if f[1] not in ('-', '') else None, f[2], f[-1])
    hits = defaultdict(list)
    for l in open(f"{os.environ['OUT']}/work/kofam_all.tbl"):
        if l.startswith('#'): continue
        f = l.split()
        t, st, d = ko.get(f[2], (None, 'full', ''))
        sc = float(f[8]) if st == 'domain' else float(f[5])
        hits[f[0]].append((f[2], sc, t, t is not None and sc >= t, d))
    for k in hits: hits[k].sort(key=lambda x: -(x[1] / x[2] if x[2] else x[1] / 1e6))
    return hits
def best_pass(hits, p):
    for k, sc, t, ok, d in hits.get(p, []):
        if ok: return k, sc, t, d
    return None
