#!/usr/bin/env bash
# 47_master_table.sh -- one table per genome, and the three questions 46 raised.
#
# Read-only apart from writing into 47_master/. Minutes.
#
# 46 VERIFIED the KOfam run:
#   0 truncated profiles of 27,864 -- the CRC error cost nothing.
#   525 ko_list entries have no profile; those are rRNA/tRNA (K01977-K01986 etc.), which
#   have no protein HMM by design.
#   CROSSCHECK PASSED: mcr in 8/8 archaea, 0/21 bacteria -- the genome-wide run
#   independently reproduces 41f's targeted result. The 4,476-KO table is trustworthy.
#
# THREE THINGS 46 RAISED, and this script answers:
#   A  sulfate reduction: 5 genomes. WHICH ONES. The whole coexistence argument is about
#      sulfate reducers outcompeting methanogens, and we have never named them.
#   B  methylamine methyltransferases: 10 genomes, but only 8 archaea exist. Two bacteria
#      hit that panel -- real, or my regex catching related corrinoid methyltransferases?
#   C  B12-dependent enzymes: 45c said 26 genomes, 46 said 8. SAME DATA, different KO
#      sets: my 46 panel regex dropped ribonucleotide reductase and diol dehydratase.
#      Printed side by side here so the thesis can state which definition it uses.
#
# And it builds MASTER.tsv: one row per genome with taxonomy, both breadths, the call,
# in-situ membership, and every functional panel. That is the table the figures come from.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/47_master"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; command -v "$PY" >/dev/null 2>&1 || PY=python3
TAB="$WS2/45_kofam_full/tables"; KOLIST="$WS2/11_kofam/ko_list"
XW="$WS2/07_mags/MAG_CROSSWALK.tsv"; XAB="$WS2/26_crossmap/cross_abundance.tsv"
TEN="$WS2/41_env_checks/in_situ_10.txt"
echo "47  $(date -Is)"

say "A. WHO are the sulfate reducers, and who does methylamines?"
"$PY" - "$KOLIST" "$TAB/per_genome_kos.json" "$XW" "$XAB" "$TEN" "$OUT" <<'PYEOF'
import sys, os, json, csv, re, collections
kolist, kojson, xw, xab, tenf, out = sys.argv[1:7]
defn={}
for line in open(kolist):
    p=line.rstrip('\n').split('\t')
    if len(p)>=3 and re.fullmatch(r'K\d{5}', p[0]): defn[p[0]]=p[-1]
per={k:set(v) for k,v in json.load(open(kojson)).items()}

bin2mag={};
for r in csv.DictReader(open(xw), delimiter='\t'):
    bin2mag[r.get('working_name','')]=r.get('mag_id','')
meta={}
for r in csv.DictReader(open(xab), delimiter='\t'):
    meta[r['mag_id']]=r
ten={l.strip() for l in open(tenf) if l.strip()} if os.path.exists(tenf) else set()

def panel(pat):
    rx=re.compile(pat, re.I)
    return {k for k,d in defn.items() if rx.search(d or '')}

def show(title, pat, detail=True):
    ks=panel(pat)
    print("\n   --- %s --- (%d KOs)" % (title, len(ks)))
    hits=[(g, sorted(ks & s)) for g,s in per.items() if ks & s]
    hits.sort(key=lambda x:-len(x[1]))
    for g,k in hits:
        m=bin2mag.get(g,'?'); md=meta.get(m,{})
        tag=' [TEN]' if m in ten else ''
        print("     %-24s %-16s %-18s breadth=%-7s %s%s"
              % (g, m, md.get('genus','-'), md.get('sed_breadth','-'),
                 ('%d KOs' % len(k)), tag))
        if detail:
            for kk in k[:6]:
                print("         %s  %s" % (kk, (defn.get(kk) or '')[:74]))
    print("     genomes: %d" % len(hits))
    return {g for g,_ in hits}

show('DISSIMILATORY SULFATE REDUCTION  (dsrAB / aprAB / sat)',
     r'dissimilatory sulfite reductase|adenylylsulfate reductase|sulfate adenylyltransferase')
show('methylamine methyltransferases (mtaB/mtbB/mtmB/mttB)',
     r'methylamine.*methyltransferase|monomethylamine|dimethylamine|trimethylamine.*methyltransferase|\[methyl-Co\(III\)')
show('methanol methyltransferase (mtaB)', r'methanol.*methyltransferase')

print("\n   --- B12-dependent enzymes: the two definitions side by side ---")
WIDE=panel(r'methionine synthase|methylmalonyl-CoA mutase|ribonucleoside-triphosphate reductase|glutamate mutase|ethanolamine ammonia-lyase|diol dehydratase')
NARROW=panel(r'methionine synthase|methylmalonyl-CoA mutase|glutamate mutase|ethanolamine ammonia-lyase')
w=[g for g,s in per.items() if WIDE & s]; n=[g for g,s in per.items() if NARROW & s]
print("     wide  (%d KOs, incl. ribonucleotide reductase + diol dehydratase): %d genomes" % (len(WIDE), len(w)))
print("     narrow(%d KOs): %d genomes" % (len(NARROW), len(n)))
print("     the difference is %d genome(s): %s" % (len(set(w)-set(n)), sorted(set(w)-set(n))[:8]))
print("     >>> 45c used the wide set (26), 46's panel used the narrow one (8).")
print("     >>> Same data. State the definition with the number.")
PYEOF

say "B. MASTER.tsv -- one row per genome"
"$PY" - "$KOLIST" "$TAB/per_genome_kos.json" "$XW" "$XAB" "$TEN" "$WS2" "$OUT" <<'PYEOF'
import sys, os, json, csv, re, glob
kolist, kojson, xw, xab, tenf, ws2, out = sys.argv[1:8]
defn={}
for line in open(kolist):
    p=line.rstrip('\n').split('\t')
    if len(p)>=3 and re.fullmatch(r'K\d{5}', p[0]): defn[p[0]]=p[-1]
per={k:set(v) for k,v in json.load(open(kojson)).items()}
bin2mag={}; qc={}
for r in csv.DictReader(open(xw), delimiter='\t'):
    b=r.get('working_name',''); bin2mag[b]=r.get('mag_id','')
    qc[b]={k:r.get(k,'') for k in ('completeness','contamination','bp','n_contigs','mimag','sample','binner','has_16S')}
meta={r['mag_id']:r for r in csv.DictReader(open(xab), delimiter='\t')}
ten={l.strip() for l in open(tenf) if l.strip()} if os.path.exists(tenf) else set()
tax={}
for f in glob.glob(os.path.join(ws2,'08_taxonomy','**','gtdbtk.*.summary.tsv'), recursive=True):
    rd=csv.DictReader(open(f), delimiter='\t'); c=rd.fieldnames or []
    nn=next((x for x in c if re.search(r'user_genome|genome|bin',x,re.I)), c[0])
    tt=next((x for x in c if re.search(r'classification|taxonom',x,re.I)), None)
    for r in rd:
        if tt: tax[r[nn]]=r[tt]
PANELS={
 'mcr':r'methyl-coenzyme M reductase',
 'methylamine_MT':r'methylamine.*methyltransferase|monomethylamine|dimethylamine|trimethylamine.*methyltransferase|\[methyl-Co\(III\)',
 'methanol_MT':r'methanol.*methyltransferase',
 'fwd_fmd':r'formylmethanofuran dehydrogenase',
 'hydrogenase':r'\[NiFe\].*hydrogenase|hydrogenase large subunit|coenzyme F420 hydrogenase',
 'fdh':r'formate dehydrogenase',
 'dsr_apr_sat':r'dissimilatory sulfite reductase|adenylylsulfate reductase|sulfate adenylyltransferase',
 'nif':r'nitrogenase',
 'ack_pta':r'acetate kinase|phosphate acetyltransferase',
 'B12_biosynth':r'cob\(|cobalamin|cobyrinate|precorrin|cobinamide|cobaltochelatase',
 'B12_dependent_wide':r'methionine synthase|methylmalonyl-CoA mutase|ribonucleoside-triphosphate reductase|glutamate mutase|ethanolamine ammonia-lyase|diol dehydratase',
 'B12_transport':r'\bbtu[BFCD]\b|cobalamin transport',
}
sets={k:{ko for ko,d in defn.items() if re.search(p, d or '', re.I)} for k,p in PANELS.items()}
def rank(t,l):
    m=re.search(r'%s__([^;]*)'%l, t or ''); return (m.group(1) or '').strip() if m else ''
rows=[]
for g in sorted(per):
    m=bin2mag.get(g,''); md=meta.get(m,{}); t=tax.get(g,'')
    r={'working_name':g,'mag_id':m,'domain':'Archaea' if 'd__Archaea' in t else 'Bacteria',
       'phylum':rank(t,'p'),'family':rank(t,'f'),'genus':rank(t,'g'),
       'sample':qc.get(g,{}).get('sample',''),'binner':qc.get(g,{}).get('binner',''),
       'completeness':qc.get(g,{}).get('completeness',''),'contamination':qc.get(g,{}).get('contamination',''),
       'bp':qc.get(g,{}).get('bp',''),'mimag':qc.get(g,{}).get('mimag',''),'has_16S':qc.get(g,{}).get('has_16S',''),
       'sed_breadth':md.get('sed_breadth',''),'liq_breadth':md.get('liq_breadth',''),
       'sed_relabund':md.get('sed_relabund',''),'liq_relabund':md.get('liq_relabund',''),
       'call':md.get('call',''),'in_situ_ten':'yes' if m in ten else 'no',
       'n_KO':len(per[g])}
    for k,s in sets.items(): r[k]=len(s & per[g])
    rows.append(r)
cols=list(rows[0].keys())
p=os.path.join(out,'MASTER.tsv')
with open(p,'w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
    for r in rows: w.writerow(r)
print("   wrote %s  (%d genomes x %d columns)" % (p, len(rows), len(cols)))
print()
hdr=['mag_id','genus','dom','ten','sed_br','nKO','mcr','maMT','dsr','H2ase','fdh','B12b']
print("   %-16s %-20s %-4s %-4s %-7s %5s %4s %5s %4s %6s %4s %5s" % tuple(hdr))
for r in sorted(rows, key=lambda r: -(float(r['sed_breadth']) if r['sed_breadth'] else 0)):
    print("   %-16s %-20s %-4s %-4s %-7s %5d %4d %5d %4d %6d %4d %5d" % (
        r['mag_id'] or '?', (r['genus'] or '-')[:20], r['domain'][:4], r['in_situ_ten'][:3],
        r['sed_breadth'] or '-', r['n_KO'], r['mcr'], r['methylamine_MT'], r['dsr_apr_sat'],
        r['hydrogenase'], r['fdh'], r['B12_biosynth']))
PYEOF
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
