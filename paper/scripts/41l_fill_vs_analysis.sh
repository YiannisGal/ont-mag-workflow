#!/usr/bin/env bash
# 41l -- verify ids against gapseq's OWN ModelSEED database, then build TWO media.
#
# 🔴 DESIGN CORRECTION. 41k said "fill on ETO_SED, then analyse on ETO_SED". That is
# wrong and would have destroyed the result:
#   * gapseq fill FORCES growth on the medium given. A hydrogenotrophic Methanofollis,
#     given no H2/acetate/methanol, would be gap-filled into a sugar fermenter.
#   * If every model is self-sufficient on the analysis medium, no member needs a
#     partner and MIP/MRO collapse toward zero BY CONSTRUCTION. The cross-feeding
#     would be gap-filled away before it could be measured.
#
# So: FILL on sediment + the currencies a member could really obtain in situ
#     (H2, acetate, formate, methanol, TMA, amino acids, B12 -- all genuinely present
#     in the lagoon, produced by other organisms).
#     ANALYSE on sediment alone, where those must come from a partner.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41l  $(date -Is)"

say "A. locate gapseq's ModelSEED compound table"
GS="${GS:-$WS2/envs/gapseq_local}"
CAND=$(find "$GS" -maxdepth 3 \( -name 'seed_metabolites*' -o -name '*metabolites*.tsv' \
       -o -name 'seed_reactions*' \) 2>/dev/null | head -5)
[ -z "$CAND" ] && CAND=$(find "$GS" -maxdepth 4 -name '*.tsv' -size +100k 2>/dev/null | head -5)
if [ -z "$CAND" ]; then no "no ModelSEED table under $GS -- set GS=<gapseq dir>"; else
  echo "$CAND" | sed 's/^/           /'
fi

say "B. verify every id by name, authoritatively"
python3 - "$OUT" $CAND <<'PY'
import sys,os,csv,gzip
out=sys.argv[1]; files=sys.argv[2:]
db={}
for f in files:
    try:
        op=gzip.open if f.endswith('.gz') else open
        with op(f,'rt',encoding='utf-8',errors='replace') as fh:
            head=fh.readline()
            delim='\t' if '\t' in head else ','
            fh.seek(0)
            for r in csv.DictReader(fh,delimiter=delim):
                idc=next((k for k in r if k and k.lower() in ('id','compound','cpd')),None)
                nmc=next((k for k in r if k and k.lower() in ('name','abbreviation')),None)
                if not idc or not nmc: break
                v=str(r[idc]).strip()
                if v.startswith('cpd'): db.setdefault(v,str(r[nmc]).strip())
    except Exception as e:
        print(f"  ...    {os.path.basename(f)}: {type(e).__name__}")
print(f"  ...    ModelSEED compounds loaded: {len(db)}")
WANT={"cpd00001":"H2O","cpd00067":"H+","cpd00009":"Phosphate","cpd00013":"Ammonium",
 "cpd00048":"Sulfate","cpd00239":"H2S","cpd00011":"CO2","cpd00242":"HCO3",
 "cpd00205":"K+","cpd00971":"Na+","cpd00254":"Mg","cpd00063":"Ca2+","cpd00099":"Cl-",
 "cpd10515":"Fe2+","cpd10516":"Fe3+","cpd00030":"Mn2+","cpd00034":"Zn2+",
 "cpd00058":"Cu2+","cpd00149":"Co2+","cpd00244":"Ni2+",
 "cpd00027":"D-Glucose","cpd00154":"D-Xylose","cpd00158":"Cellobiose",
 "cpd00280":"Galacturonate","cpd00076":"Sucrose?","cpd01200":"Cellulose?",
 "cpd00029":"Acetate","cpd11640":"H2","cpd00047":"Formate","cpd00116":"Methanol",
 "cpd00441":"Trimethylamine","cpd00098":"Choline","cpd00166":"B12?","cpd03424":"B12?",
 "cpd00007":"O2","cpd00209":"Nitrate"}
if not db:
    print("  [FAIL] no compound table parsed -- cannot verify. STOP.")
    sys.exit(1)
bad=[]
print(f"  ...    {'id':<12}{'expected':<18}{'ModelSEED says'}")
for cid,exp in WANT.items():
    nm=db.get(cid)
    if nm is None:
        print(f"           {cid:<12}{exp:<18}-- NOT IN DATABASE"); bad.append((cid,exp,None)); continue
    e=exp.rstrip('?').lower(); n=nm.lower()
    okm = e in n or n in e or e.replace('-','') in n.replace('-','')
    print(f"           {cid:<12}{exp:<18}{nm}{'' if okm else '   <-- MISMATCH'}")
    if not okm and not exp.endswith('?'): bad.append((cid,exp,nm))
if bad:
    print(f"  [FAIL] {len(bad)} id(s) do not name what I assumed:")
    for c,e,n in bad: print(f"           {c}: assumed {e!r}, database says {n!r}")
else:
    print("  [OK]   every assumed id names the expected compound")
open(os.path.join(out,'modelseed_id_check.tsv'),'w').write(
    "cpd\tassumed\tmodelseed_name\n"+"".join(
      f"{c}\t{e}\t{db.get(c,'ABSENT')}\n" for c,e in WANT.items()))
print(f"  [OK]   wrote {out}/modelseed_id_check.tsv")
PY

say "C. build the two media"
python3 - "$OUT" <<'PY'
import sys,os
out=sys.argv[1]
chk=os.path.join(out,'modelseed_id_check.tsv')
name={}
if os.path.exists(chk):
    for i,l in enumerate(open(chk)):
        if i==0: continue
        p=l.rstrip('\n').split('\t')
        if len(p)>=3 and p[2]!='ABSENT': name[p[0]]=p[2]
def nm(c,d): return name.get(c,d)
SED=["cpd00001","cpd00067","cpd00009","cpd00013","cpd00048","cpd00239","cpd00011",
     "cpd00242","cpd00205","cpd00971","cpd00254","cpd00063","cpd00099","cpd10515",
     "cpd10516","cpd00030","cpd00034","cpd00058","cpd00149","cpd00244"]
CARB=["cpd00027","cpd00154","cpd00158","cpd00280"]     # glucose, xylose, cellobiose, galacturonate
CURRENCY=["cpd00029","cpd11640","cpd00047"]            # acetate, H2, formate
ROUTE=["cpd00116","cpd00441","cpd00098"]               # methanol, TMA, choline
DEP=["cpd00166"]                                       # B12
AA=["cpd00023","cpd00033","cpd00035","cpd00039","cpd00041","cpd00051","cpd00053",
    "cpd00054","cpd00060","cpd00065","cpd00066","cpd00069","cpd00084","cpd00107",
    "cpd00119","cpd00129","cpd00132","cpd00156","cpd00161","cpd00322"]
carb=[c for c in CARB if c in name] or ["cpd00027","cpd00154"]
analysis = SED+carb
fill     = SED+carb+CURRENCY+ROUTE+DEP+AA
BANNED={"cpd00007","cpd00209","cpd00075","cpd00659","cpd00418"}
for lab,s in (("ETO_FILL",fill),("ETO_SED",analysis)):
    b=[c for c in s if c in BANNED]
    if b: print(f"  [FAIL] {lab} contains oxidant(s) {b}"); sys.exit(1)
print("  [OK]   neither medium contains O2, nitrate or nitrite")
leak=[c for c in analysis if c in set(CURRENCY+ROUTE+DEP+AA)]
if leak: print(f"  [FAIL] analysis medium contains a currency/route/B12/AA: {leak}"); sys.exit(1)
print("  [OK]   analysis medium contains no currency, route compound, B12 or amino acid")
if not set(analysis) <= set(fill):
    print("  [FAIL] analysis medium is not a subset of the fill medium"); sys.exit(1)
print(f"  [OK]   ETO_SED ({len(analysis)}) is a strict subset of ETO_FILL ({len(fill)})")
print(f"  ...    the difference is {len(set(fill)-set(analysis))} compounds -- exactly what"
      " SMETANA must explain by cross-feeding")
def w(n,ids):
    p=os.path.join(out,f"{n}-medium.csv")
    with open(p,'w') as fh:
        fh.write("compounds,name,maxFlux\n")
        for c in sorted(set(ids)): fh.write(f"{c},{nm(c,c)},10\n")
    print(f"  [OK]   {p}  ({len(set(ids))})")
w("ETO_FILL",fill); w("ETO_SED",analysis)
w("ETO_SED_NOSO4",[c for c in analysis if c!="cpd00048"])
p=os.path.join(out,"etoliko_media.tsv")
with open(p,'w') as fh:
    fh.write("medium\tdescription\tcompound\tname\n")
    for n2,ids in (("ETO_SED",analysis),
                   ("ETO_SED_NOSO4",[c for c in analysis if c!="cpd00048"])):
        for c in sorted(set(ids)): fh.write(f"{n2}\tEtoliko sediment, anoxic\t{c}\t{nm(c,c)}\n")
print(f"  [OK]   {p}   (analysis media only -- ETO_FILL is for gapseq, not SMETANA)")
print("\n  ...    carbon sources in the analysis medium: "+", ".join(f"{c} ({nm(c,c)})" for c in carb))
if len(carb)<4:
    print("  [WARN] cellobiose and/or galacturonate could not be id-verified and were")
    print("         dropped. Those are the cellulose and PECTIN arms -- pectin is where")
    print("         the methanol in the headline route comes from. Resolve before the run.")
PY
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
