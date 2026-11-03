"""SP01 step 6-7 -- markers_long.tsv, markers_wide.tsv (with M00087, nif tree, hydrogenase and dsr-tree columns), versions.tsv.
Absence is written only for genomes >=90% complete; otherwise 'NA (<90% complete)'.
Requires markers/QC.txt == PASS (the spec forbids producing the matrix if a control failed)."""
import os, csv, subprocess, hashlib, datetime
from collections import defaultdict
OUT, WS2 = os.environ['OUT'], os.environ['WS2']
qc = open(f'{OUT}/markers/QC.txt').read()
assert 'QC_STATUS: PASS' in qc, 'QC controls failed -- matrix not produced'
rows = list(csv.DictReader(open(f'{OUT}/markers/calls.tsv'), delimiter='\t'))
genomes = list(dict.fromkeys(r['genome'] for r in rows))
G = {r['genome']: r for r in rows}
cell = {(r['genome'], r['KO']): r for r in rows}

# ---------------- markers_long.tsv
LONG = ['genome', 'source', 'KO', 'gene', 'bitscore', 'threshold', 'pass', 'copies', 'contig', 'locus', 'completeness',
        'in_situ_breadth', 'score_type', 'profile_len', 'fragments', 'top_KO_of_called', 'domain', 'definition']
with open(f'{OUT}/markers_long.tsv', 'w') as o:
    o.write('\t'.join(LONG) + '\n')
    for r in rows:
        r = dict(r); r['pass'] = r.pop('passes')
        o.write('\t'.join(str(r[c]) for c in LONG) + '\n')

# ---------------- M00087
S1 = ['K00232', 'K00249', 'K00255', 'K06445', 'K09479']
S2_multi = ['K01825', 'K01782', 'K07514', 'K07515', 'K10527']; S2_hyd = ['K01692', 'K07511', 'K13767']; S2_had = ['K00022', 'K07516']
S3 = ['K00632', 'K07508', 'K07509', 'K07513']
def pres(g, k): return cell[(g, k)]['passes'] == 'present'
m87 = {}
with open(f'{OUT}/markers/M00087_completeness.tsv', 'w') as o:
    o.write('genome\tstep1_acylCoA_DH\tstep2_hydratase+HADH\tstep3_thiolase\tsteps_complete\tcompleteness\tcall\n')
    for g in genomes:
        s1 = [k for k in S1 if pres(g, k)]
        s2m = [k for k in S2_multi if pres(g, k)]; s2h = [k for k in S2_hyd if pres(g, k)]; s2d = [k for k in S2_had if pres(g, k)]
        s3 = [k for k in S3 if pres(g, k)]
        st2 = bool(s2m) or (bool(s2h) and bool(s2d))
        n = int(bool(s1)) + int(st2) + int(bool(s3))
        comp = float(G[g]['completeness'])
        call = 'complete (3/3)' if n == 3 else (f'{n}/3' + (' (NA: <90% complete, absence not callable)' if comp < 90 else ''))
        s2txt = ','.join(s2m) if s2m else (f"hyd:{','.join(s2h) or '-'} had:{','.join(s2d) or '-'}")
        m87[g] = (n, call)
        o.write(f"{g}\t{','.join(s1) or '-'}\t{s2txt}\t{','.join(s3) or '-'}\t{n}\t{n/3:.2f}\t{call}\n")

# ---------------- nif, hydrogenases, dsr tree
nif = defaultdict(lambda: {'H': [], 'D': []})
for r in csv.DictReader(open(f'{OUT}/nif_tree_verdicts.tsv'), delimiter='\t'):
    nif[r['genome']]['H' if r['marker'].startswith('nifH') else 'D'].append(r['verdict'])
def nifcall(g):
    h, d = nif[g]['H'], nif[g]['D']
    if not h and not d: return 'no NifH/NifD above threshold'
    hn = [v for v in h if v.startswith('NITROGENASE')]; dn = [v for v in d if v.startswith('NITROGENASE')]
    k = pres(g, 'K02591')
    if hn and dn: return f"nitrogenase H+D in Nif/Anf/Vnf clades ({hn[0].split(': ')[1]}); nifK {'present' if k else 'absent'}"
    parts = [f"{x}:{v.split(': ')[-1]}" for x, vs in (('H', h), ('D', d)) for v in vs]
    return 'NOT nitrogenase / incomplete: ' + '; '.join(parts)
hyd = defaultdict(set)
for r in csv.DictReader(open(f'{OUT}/hydrogenase_classes.tsv'), delimiter='\t'):
    tag = f"{r['class']} {r['group']}" + ('' if r['confidence'] == 'assigned' else '?')
    hyd[r['genome']].add(tag)
DSR = {'ETO_LIQ_MAG_022': 'no dsrAB binned; unbinned LIQ_contig_19299 reductive dsrAB (Desulforamulus/Desulfofalx) attributed ONLY as an open hypothesis (step 59)'}
def dsrcall(g):
    base = 'dsrAB above threshold' if pres(g, 'K11180') or pres(g, 'K11181') else 'no dsrAB in genome'
    return DSR.get(g, base + ('' if G[g]['source'] == 'isolate' else '; steps 55-59: both reductive dsrAB operons are on unbinned contigs'))

# ---------------- markers_wide.tsv
kos = list(dict.fromkeys(r['KO'] for r in rows if r['gene'] != 'M00087_only'))
gene = {r['KO']: r['gene'] for r in rows}
def v(g, k):
    p = cell[(g, k)]['passes']
    return {'present': cell[(g, k)]['copies'], 'absent': '0', 'fragment_only': 'frag'}.get(p, 'NA')
with open(f'{OUT}/markers_wide.tsv', 'w') as o:
    o.write('# cell = copies above curated KOfam threshold (complete genes); 0 = absent in a >=90%-complete genome; NA = not found in a <90%-complete genome (absence not callable); frag = only truncated (partial/contig-edge) hits above threshold. Length rule: user decision 2026-09-29, option 1.\n')

    o.write('\t'.join(['genome', 'source', 'domain', 'completeness'] + [f'{gene[k]}_{k}' for k in kos] +
                      ['M00087_beta_oxidation', 'nif_tree_verdict', 'hydrogenase_groups', 'dsr_tree']) + '\n')
    for g in genomes:
        o.write('\t'.join([g, G[g]['source'], G[g]['domain'], G[g]['completeness']] + [v(g, k) for k in kos] +
                          [m87[g][1], nifcall(g), ';'.join(sorted(hyd[g])) or 'none', dsrcall(g)]) + '\n')

# ---------------- versions.tsv
def sh(c):
    try: return subprocess.run(c, shell=True, capture_output=True, text=True).stdout.strip().splitlines()[0]
    except Exception as e: return f'unavailable ({e})'
def md5(p): return hashlib.md5(open(p, 'rb').read()).hexdigest()
KOL = os.environ['KOLIST']; PIXI = os.environ['PIXI']
V = [
 ('annotation_engine', 'HMMER hmmsearch over KOfam profiles + ko_list thresholds (step-45 method; KofamScan exec_annotation NOT used; approved deviation)'),
 ('hmmsearch', sh('/usr/bin/hmmsearch -h | sed -n 2p')),
 ('ko_list', f'{KOL}  md5={md5(KOL)}  mtime={datetime.datetime.fromtimestamp(os.path.getmtime(KOL)).date()}'),
 ('kofam_profiles', f"{WS2}/11_kofam/profiles.tar.gz extracted at {os.environ['PROFILES']} (27,864 profiles; release date not recorded upstream)"),
 ('diamond', sh('/usr/local/bin/diamond version')),
 ('mafft', sh(f'{PIXI}/mafft --version 2>&1')),
 ('FastTree', sh(f'{PIXI}/FastTree </dev/null 2>&1 | head -1')),
 ('hydrogenase_refs', f"Greening lab HydDB-labelled sets, figshare 13152518 (NiFe, md5 {md5(OUT+'/hyd/NiFe_hydrogenase.fasta')}) and 13152542 (FeFe, md5 {md5(OUT+'/hyd/FeFe_hydrogenase.fasta')}), published 2020-12-06; HydDB web unreachable 2026-09-29"),
 ('nif_refs', 'UniProt REST queries in scripts/s4_nif.sh, fetched 2026-09-29 (nif/refs.tsv)'),
 ('urease_refs', 'UniProt gene:ureC (name-filtered) + reviewed pyrC/hutI/allB decoys, fetched 2026-09-29; controls copied from 34_unbinned/marker_ref.faa'),
 ('M00087', 'KEGG REST get/M00087, fetched 2026-09-29'),
 ('python_env', f"{WS2}/envs/smetana (python {sh(WS2 + '/envs/smetana/bin/python -V 2>&1')}, pandas unused)"),
 ('pixi_env', PIXI.rsplit('/bin', 1)[0] + ' (gapseq_project default; mafft, FastTree)'),
]
with open(f'{OUT}/versions.tsv', 'w') as o:
    o.write('item\tversion_or_source\n')
    for a, b in V: o.write(f'{a}\t{b}\n')
print(open(f'{OUT}/versions.tsv').read())
print('long rows', len(rows), 'genomes', len(genomes))
