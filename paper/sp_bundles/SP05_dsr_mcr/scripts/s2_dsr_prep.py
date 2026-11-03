"""SP05 step 2a -- label step-58 references and contig proteins as DsrA (K11180) or DsrB (K11181) by KOfam, write tree inputs.
Reference role: Allochromatium / Chlorobium = OXIDATIVE (reverse-type) controls; all others = REDUCTIVE."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__)); import kofam
OUT = os.environ['OUT']
H = kofam.load()
def rf(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith('>'): n = l[1:].strip(); d[n] = []
        else: d[n].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}
refs = rf(f'{OUT}/work/dsr_refs.faa'); q = rf(f'{OUT}/work/contig_proteins.faa')
OXI = ('Allochromatium', 'Chlorobium')
out = {'dsrA': [], 'dsrB': []}; lab = open(f'{OUT}/work/dsr_labels.tsv', 'w')
lab.write('seq\tsource\trole\tassigned\tKO_score\tthreshold\n')
def which(p):
    s = {k: sc for k, sc, t, ok, d in H.get(p, []) if k in ('K11180', 'K11181')}
    if not s: return None, None
    k = max(s, key=s.get); return ('dsrA' if k == 'K11180' else 'dsrB'), s[k]
for k, v in refs.items():
    m, sc = which(k)
    if not m: continue
    org = k.split('|')[1]; role = 'OXIDATIVE' if org.startswith(OXI) else 'REDUCTIVE'
    name = f"REF_{role[:3]}|{org}|{k.split('|')[2]}"
    out[m].append((name, v)); lab.write(f'{k}\tref\t{role}\t{m}\t{sc:.1f}\t-\n')
for k, v in q.items():
    bp = kofam.best_pass(H, k)
    if bp and bp[0] in ('K11180', 'K11181'):
        m = 'dsrA' if bp[0] == 'K11180' else 'dsrB'
        out[m].append((f'QUERY|{k}', v)); lab.write(f'{k}\tquery\t-\t{m}\t{bp[1]:.1f}\t{bp[2]}\n')
lab.close()
for m, L in out.items():
    open(f'{OUT}/trees/{m}_in.faa', 'w').write(''.join(f'>{n}\n{s}\n' for n, s in L))
    print(m, 'refs', sum(n.startswith('REF') for n, _ in L), 'queries', [n for n, _ in L if n.startswith('QUERY')])
if not out['dsrA'] or not any(n.startswith('QUERY') for n, _ in out['dsrA'] + out['dsrB']): sys.exit('no DsrA/DsrB queries found')
