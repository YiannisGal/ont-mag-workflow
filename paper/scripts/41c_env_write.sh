#!/usr/bin/env bash
# 41c -- resolve compounds BY ModelSEED ID (verified against the model's own names),
# work out how SMETANA addresses exchange reactions, then write the environment file.
#
# 41b established: namespace is ModelSEED cpd#####. Its name matching then failed for a
# format reason -- gapseq names carry a compartment suffix ('phosphate-c0') and every
# compound exists in both c0 and e0, so EXACT never fired and SUBSTR returned both plus
# unrelated hits ('h2' -> 'h2o', 'sulfate' -> '3-phosphoadenylylsulfate').
# A medium concerns EXTRACELLULAR compounds only, so e0 is the only compartment that matters.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; mkdir -p "$OUT"
CULT="${CULT:-$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis}"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41c  $(date -Is)"

M="$CULT/gapseq_project/YGCA01_Clean.xml"
[ -f "$M" ] || { no "model not found: $M"; exit 1; }

say "A. compounds by ModelSEED ID, name-verified"
python3 - "$M" "$OUT" <<'PY'
import re,sys,os
m,out=sys.argv[1],sys.argv[2]
t=open(m,encoding='utf-8',errors='replace').read()
sp=dict(re.findall(r'<species[^>]*\bid="([^"]+)"[^>]*\bname="([^"]*)"',t))
def norm(n): return re.sub(r'-[ce]\d+$','',(n or '').strip()).lower()
ext={k:norm(v) for k,v in sp.items() if re.search(r'_e\d+$',k)}
print(f"  ...    {len(sp)} species total, {len(ext)} extracellular (_e#)")
# ID is the assertion; the name is the CHECK on that assertion.
SPEC=[("cpd00001","h2o",        "water",        "base"),
      ("cpd00067","h+",         "proton",       "base"),
      ("cpd00009","phosphate",  "phosphate",    "non-limiting"),
      ("cpd00013","nh3",        "ammonium",     "non-limiting"),
      ("cpd00048","sulfate",    "sulfate",      "non-limiting"),
      ("cpd00239","h2s",        "sulfide",      "present"),
      ("cpd00011","co2",        "CO2",          "present"),
      ("cpd00242","hco3","bicarbonate","present"),
      ("cpd00205","k+",         "potassium",    "base"),
      ("cpd00971","na+",        "sodium",       "base"),
      ("cpd00254","mg",         "magnesium",    "base"),
      ("cpd00063","ca2+",       "calcium",      "base"),
      ("cpd00099","cl-",        "chloride",     "base"),
      ("cpd10515","fe2+",       "iron(II)",     "trace"),
      ("cpd00030","mn2+",       "manganese",    "trace"),
      ("cpd00034","zn2+",       "zinc",         "trace"),
      ("cpd00058","cu2+",       "copper",       "trace"),
      ("cpd00149","co2+",       "cobalt",       "trace"),
      ("cpd00244","ni2+",       "nickel",       "trace"),
      ("cpd11574","molybdate",  "molybdate",    "trace"),
      ("cpd00029","acetate",    "acetate",      "product"),
      ("cpd11640","h2",         "H2",           "product"),
      ("cpd00154","xylose",     "xylose",       "CARBON"),
      ("cpd00158","cellobiose", "cellobiose",   "CARBON"),
      ("cpd00027","d-glucose",  "glucose",      "CARBON"),
      ("cpd00280","galacturonate","galacturonate","CARBON"),
      ("cpd00116","methanol",   "methanol",     "EXCLUDE"),
      ("cpd00441","trimethylamine","TMA",       "EXCLUDE"),
      ("cpd00098","choline",    "choline",      "EXCLUDE"),
      ("cpd00166","cobalamin",  "B12",          "EXCLUDE")]
rows=[];bad=[];absent=[]
print(f"  ...    {'cpd id':<12}{'expected name':<18}{'found (e0)':<22}{'role'}")
for cid,expect,label,role in SPEC:
    hits=[k for k in ext if k.startswith(f"M_{cid}_e") or k.startswith(f"{cid}_e")]
    if not hits:
        absent.append((cid,label)); print(f"           {cid:<12}{expect:<18}{'-- not in this model':<22}{role}"); continue
    nm=ext[hits[0]]
    match = expect in nm or nm in expect
    print(f"           {cid:<12}{expect:<18}{nm:<22}{role}   {'ok' if match else 'MISMATCH'}")
    if not match: bad.append((cid,expect,nm))
    rows.append((cid,nm,role))
if bad:
    print(f"  [FAIL] {len(bad)} ID/name mismatch(es) -- the ID table is wrong, not the model:")
    for c,e,n in bad: print(f"           {c}: expected ~{e!r}, model says {n!r}")
else:
    print("  [OK]   every ID present in this model carries the expected name")
if absent:
    print(f"  ...    {len(absent)} compound(s) absent from THIS isolate model:")
    print("         "+", ".join(f"{c} ({l})" for c,l in absent))
    print("         Expected: one Sporosarcina isolate cannot contain every compound the")
    print("         29-MAG community needs. ModelSEED IDs are global, so absence here is")
    print("         not evidence against the ID -- but it is NOT verification either.")
    print("  [WARN] those IDs remain UNVERIFIED until the real MAG models exist.")
open(os.path.join(out,'compound_id_resolution.tsv'),'w').write(
  "cpd_id\tname_in_model\trole\n"+"".join(f"{c}\t{n}\t{r}\n" for c,n,r in rows))
print(f"  [OK]   wrote {out}/compound_id_resolution.tsv")
PY

say "B. how does SMETANA address an exchange reaction?"
note "gapseq writes EX_cpd00013_e0; SMETANA's media db assumes a convention. Which?"
SM=""
for c in "$WS2/envs/smetana/bin/smetana" "$(command -v smetana 2>/dev/null)"; do
  [ -x "$c" ] && { SM="$c"; break; }; done
if [ -z "$SM" ]; then no "smetana not found"; else
  ok "using $SM"
  "$SM" --help 2>&1 | grep -iE 'flavor|exchange|format|media|mediadb' | sed 's/^/           /'
  note "reframed's own view of this model's exchanges:"
  python3 - "$M" <<'PY' 2>&1 | sed 's/^/           /'
try:
    from reframed import load_cbmodel
    m=load_cbmodel(__import__('sys').argv[1], flavor='fbc2')
    ex=[r for r in m.reactions if r.startswith('R_EX') or r.startswith('EX')]
    print(f"reframed loaded the model: {len(m.reactions)} reactions, {len(ex)} exchanges")
    for r in ex[:5]: print(f"  {r}")
    print("-> the media db 'compound' column must produce exactly these ids")
except Exception as e:
    print(f"reframed could not load it: {type(e).__name__}: {e}")
    print("-> resolve before any SMETANA run; a model it cannot read yields no result")
PY
fi

say "C. methanogenesis markers, with the names Bakta actually uses"
note "41b scored mcrA=0 for all eight archaea including known Methanosarcina."
note "MCR's enzyme name is 'coenzyme-B sulfoethylthiotransferase' -- that was a"
note "naming artifact in my screen, not an absence in the genomes."
python3 - "$WS2" <<'PY'
import sys,os,csv,re
ws2=sys.argv[1]
tax=os.path.join(ws2,'07_mags','MAG_TAXONOMY.tsv')
rows=list(csv.DictReader(open(tax),delimiter='\t'))
nk=next(k for k in rows[0] if k.strip().lower()=='working_name')
dk=next(k for k in rows[0] if k.strip().lower() in ('domain','kingdom'))
gk=next(k for k in rows[0] if k.strip().lower()=='genus')
arch=[r for r in rows if 'archae' in (r[dk] or '').lower()]
MK={'mcrABG':[r'methyl-coenzyme m reductase',r'coenzyme-b sulfoethylthiotransferase',r'\bmcr[abg]\b'],
    'mtaBC  methanol':[r'methanol.*methyltransferase',r'\bmta[bc]\b'],
    'mttBC  TMA':[r'trimethylamine.*methyltransferase',r'\bmtt[bc]\b'],
    'mtmBC  MMA':[r'monomethylamine.*methyltransferase',r'\bmtm[bc]\b'],
    'mtbBC  DMA':[r'dimethylamine.*methyltransferase',r'\bmtb[bc]\b'],
    'fwd/fmd':[r'formylmethanofuran dehydrogenase'],
    'hdr':[r'heterodisulfide reductase'],
    'pyl  amber-tRNA':[r'pyrrolysine',r'\bpyl[bcds]\b']}
print(f"           {'genome':<24}{'genus':<20}"+"".join(f"{k.split()[0]:>10}" for k in MK))
meth={}
for r in arch:
    g=r[nk]; f=os.path.join(ws2,'09_annotation',g,g+'.tsv')
    if not os.path.exists(f): print(f"           {g:<24}(no table)"); continue
    txt=open(f,encoding='utf-8',errors='replace').read().lower()
    line=f"           {g:<24}{r[gk]:<20}"; d={}
    for k,ps in MK.items():
        n=sum(len(re.findall(p,txt)) for p in ps); d[k]=n; line+=f"{n if n else '.':>10}"
    print(line); meth[g]=(r[gk],d)
print()
mcr=[g for g,(_,d) in meth.items() if d['mcrABG']>0]
print(f"  [{'OK' if len(mcr)>=6 else 'FAIL'}]   mcr found in {len(mcr)}/{len(arch)} archaea"
      + ("" if len(mcr)>=6 else "  -- still a naming problem, not biology; needs KOfam"))
mt=[(g,gen,d) for g,(gen,d) in meth.items()
    if d['mttBC  TMA']+d['mtmBC  MMA']+d['mtbBC  DMA']+d['mtaBC  methanol']>0]
print(f"  ...    {len(mt)} archaeon(s) carry a methylotrophic methyltransferase:")
for g,gen,d in mt:
    which=[k.split()[0] for k in ('mtaBC  methanol','mttBC  TMA','mtmBC  MMA','mtbBC  DMA') if d[k]]
    print(f"           {g:<24}{gen:<20}{', '.join(which)}")
print("  ...    Cross these against sed_breadth: the methylotrophic CONSUMER question")
print("         in Part 1 turns on whether any of them is in the in-situ set.")
PY

say "D. aromatic hits: degradation or biosynthesis?"
note "4-hydroxybenzoate and catechol appear in UBIQUINONE BIOSYNTHESIS (ubiC/ubiD/ubiA)."
note "41b's 52 loci almost certainly mix biosynthesis into a degradation claim."
python3 - "$WS2" "$OUT" <<'PY'
import sys,os,csv,re,collections
ws2,out=sys.argv[1],sys.argv[2]
tax=os.path.join(ws2,'07_mags','MAG_TAXONOMY.tsv')
rows=list(csv.DictReader(open(tax),delimiter='\t'))
nk=next(k for k in rows[0] if k.strip().lower()=='working_name')
DEG=re.compile(r'benzoyl-coa reductase|\bbam[bcdefghi]\b|\bbcr[abcd]\b|'
               r'benzoate-*-*coa ligase|\bbad[ai]\b|phenylphosphate (carboxylase|synthase)|'
               r'\bppc[ab]\b|cyclohexa-1,5-diene-1-carbonyl',re.I)
BIO=re.compile(r'ubiquinone|\bubi[a-z]\b|chorismate|menaquinone|\bmen[a-z]\b|'
               r'shikimat|phenylalanine|tyrosine|tryptophan|folate',re.I)
deg=collections.Counter(); bio=collections.Counter(); ex=collections.defaultdict(list)
for r in rows:
    g=r[nk]; f=os.path.join(ws2,'09_annotation',g,g+'.tsv')
    if not os.path.exists(f): continue
    seen=set()
    for line in open(f,encoding='utf-8',errors='replace'):
        p=line.rstrip('\n').split('\t')
        if len(p)<2 or p[0] in seen: continue
        if DEG.search(line) and not BIO.search(line):
            deg[g]+=1; seen.add(p[0]); ex[g].append(line.strip()[:110])
        elif BIO.search(line): bio[g]+=1
print(f"  ...    STRICT anaerobic-aromatic-degradation loci: {sum(deg.values())} across {len(deg)} genomes")
if deg:
    for g,n in deg.most_common(12):
        print(f"           {n:>3}  {g}")
        for e in ex[g][:2]: print(f"                {e}")
else:
    print("           none")
print(f"  ...    biosynthesis-context loci excluded: {sum(bio.values())}")
ins=os.path.join(out,'in_situ_10.txt')
xw=os.path.join(ws2,'MAG_CROSSWALK.tsv')
if os.path.exists(ins) and os.path.exists(xw):
    xr=list(csv.DictReader(open(xw),delimiter='\t'))
    cols=list(xr[0].keys())
    ec=next((k for k in cols if any(str(x[k]).startswith('ETO_') for x in xr)),None)
    wc=next((k for k in cols if any(str(x[k]) in {r[nk] for r in rows} for x in xr)),None)
    m={x[ec]:x[wc] for x in xr} if ec and wc else {}
    tgt={m.get(i.strip(),i.strip()) for i in open(ins)}
    hit={g:n for g,n in deg.items() if g in tgt}
    print(f"  ...    of the in-situ 10: {len(hit)} carry a strict degradation locus")
    for g,n in sorted(hit.items(),key=lambda x:-x[1]): print(f"           {n:>3}  {g}")
    print("  ...    VERDICT: a curated KOfam/HMM pass on bamBCDEFGHI + bcrABCD is"
          if hit else "  ...    VERDICT: no strict degradation evidence in the in-situ set;")
    print("         justified." if hit else "         the phenol pool is not obviously used by the resident community.")
open(os.path.join(out,'aromatic_strict.tsv'),'w').write(
  "genome\tstrict_degradation_loci\n"+"".join(f"{g}\t{n}\n" for g,n in deg.most_common()))
print(f"  [OK]   wrote {out}/aromatic_strict.tsv")
PY
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
