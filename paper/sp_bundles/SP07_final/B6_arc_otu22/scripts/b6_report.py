import hashlib, re, os
os.chdir('work')
def fa(p):
    s=''.join(l.strip() for l in open(p) if not l.startswith('>')); return s
def row(p):
    if not os.path.getsize(p): return None
    return open(p).read().split('\n')[0].split('\t')
def span(a):
    t=open(f'{a}_iddef2.aln').read(); m=re.findall(r'^Tgt\s+(\d+)\s+[+-]\s+\S+\s+(\d+)',t,re.M); return m[0][0], m[-1][1]
gene=fa('gene.fa'); out=[]; ctrl=row('ctrl.tsv')
L=['amplicon\tmd5_primerfree\tlength_primerfree\tidentical_positions\tidentical_over_length\tiddef0_pct\tiddef1_pct\tiddef2_pct\taln_cols\tmism\tgaps\tgene_from\tgene_to']
for a in ('ms','server'):
    s=fa(f'{a}_pf.fa'); r={d:row(f'{a}_iddef{d}.tsv') for d in (0,1,2)}
    ids=int(r[2][4]) if r[2] else 0
    L.append('\t'.join(map(str,[f'Arc_Otu22_{a}',hashlib.md5(s.encode()).hexdigest(),len(s),ids,f'{ids}/{len(s)} ({100*ids/len(s):.1f}%)',
        *[(r[d][2] if r[d] else 'no hit') for d in (0,1,2)], r[2][3], r[2][5], r[2][6], *span(a)])))
open('../B6_result.tsv','w').write('\n'.join(L)+'\n')
same=fa('server_pf.fa')==fa('ms_pf.fa'); sv=row('server_vs_ms.tsv')
# window 841-850 of the gene, base by base, from the iddef2 alignment of each amplicon
def window(a,lo=835,hi=856):
    txt=open(f'{a}_iddef2.aln').read(); q=t=''; qs=ts=None
    for m in re.finditer(r'^Qry\s+(\d+)\s+[+-]\s+(\S+)\s+\d+\s*\n.*\nTgt\s+(\d+)\s+[+-]\s+(\S+)',txt,re.M):
        if qs is None: qs,ts=int(m.group(1)),int(m.group(3))
        q+=m.group(2); t+=m.group(4)
    pos=ts; cols=[]
    for qc,tc in zip(q,t):
        if lo<=pos<=hi: cols.append((pos if tc!='-' else f'{pos-1}+',tc,qc))
        if tc!='-': pos+=1
    return cols
W=['# Gene window 835-856 (target = full-length gene coordinates); columns from vsearch iddef2 alignments','gene_pos\tgene\tms_amplicon\tserver_Otu22']
wm,ws=window('ms'),window('server')
for (p,g,m),(_,_,s) in zip(wm,ws): W.append(f'{p}\t{g}\t{m}\t{s}')
open('../B6_window_841_850.tsv','w').write('\n'.join(W)+'\n')
qc=f"positive control gene vs itself: {ctrl[2] if ctrl else 'FAIL'}% ; server Otu22 vs ms amplicon (full, with primers): {sv[2] if sv else 'no hit'}% ({sv[3] if sv else '-'} cols, {sv[4] if sv else '-'} mism); primer-free identical: {same}"
open('../QC.txt','w').write(qc+'\n'); print(qc); print(open('../B6_result.tsv').read()); print('\n'.join(W))
