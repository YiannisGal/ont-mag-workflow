#!/usr/bin/env bash
# 41m -- confirm every compound id from the ModelSEED table's FULL row.
#
# 41l flagged 12 mismatches. All 12 are spurious: the column it read is ModelSEED's
# BiGG-style abbreviation ('pi', 'so4', 'glc-D'), not the full name. But "I recognise
# those abbreviations" is not verification, so this prints every column for every id
# and matches against name AND abbreviation AND aliases.
#
# 41l also printed [FAIL] in section B and then built the media in section C anyway.
# The gate did not gate. Here the build is INSIDE the check.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; mkdir -p "$OUT"; GS="${GS:-$WS2/envs/gapseq_local}"
say(){ printf '\n=== %s ===\n' "$*"; }
echo "41m  $(date -Is)"
SEED="$GS/dat/seed_metabolites.tsv"
[ -f "$SEED" ] || { echo "  [FAIL] $SEED not found"; exit 1; }

say "A. full ModelSEED record for every id in the media"
python3 - "$SEED" "$OUT" <<'PY'
import sys,csv,os
seed,out=sys.argv[1],sys.argv[2]
rows={}
with open(seed,encoding='utf-8',errors='replace') as fh:
    rd=csv.DictReader(fh,delimiter='\t')
    cols=rd.fieldnames
    for r in rd:
        i=r.get('id') or r.get('compound')
        if i: rows[i.strip()]=r
print(f"  ...    columns in seed_metabolites.tsv: {', '.join(cols)}")
print(f"  ...    {len(rows)} compounds")

# id -> the compound it MUST be, given as several acceptable spellings
SPEC={
 "cpd00001":["h2o","water"],              "cpd00067":["h","h+","proton"],
 "cpd00009":["pi","phosphate","orthophosphate"], "cpd00013":["nh4","ammonium","ammonia","nh3"],
 "cpd00048":["so4","sulfate"],            "cpd00239":["h2s","hydrogen sulfide","sulfide"],
 "cpd00011":["co2","carbon dioxide"],     "cpd00242":["hco3","bicarbonate"],
 "cpd00205":["k","k+","potassium"],       "cpd00971":["na1","na+","sodium"],
 "cpd00254":["mg2","mg","magnesium"],     "cpd00063":["ca2","ca2+","calcium"],
 "cpd00099":["cl","cl-","chloride"],      "cpd10515":["fe2","fe2+"],
 "cpd10516":["fe3","fe3+"],               "cpd00030":["mn2","mn2+"],
 "cpd00034":["zn2","zn2+"],               "cpd00058":["cu2","cu2+"],
 "cpd00149":["cobalt2","co2+","cobalt"],  "cpd00244":["ni2","ni2+","nickel"],
 "cpd00027":["glc-d","d-glucose","glucose"], "cpd00154":["xyl-d","d-xylose","xylose"],
 "cpd00158":["cellb","cellobiose"],       "cpd00280":["galur","galacturonate","d-galacturonate"],
 "cpd00029":["ac","acetate"],             "cpd11640":["h2","hydrogen"],
 "cpd00047":["for","formate"],            "cpd00116":["meoh","methanol"],
 "cpd00441":["tma","trimethylamine"],     "cpd00098":["chol","choline"],
 "cpd00166":["adocbl","adenosylcobalamin","cobalamin"],
 "cpd03424":["vitamin b12","cbl","cobalamin"],
 "cpd00007":["o2","oxygen"],              "cpd00209":["no3","nitrate"],
}
def fields(r):
    v=[]
    for k in ('name','abbreviation','formula','aliases'):
        if r.get(k): v.append(str(r[k]))
    return v
bad=[]
for cid,accept in SPEC.items():
    r=rows.get(cid)
    if not r: bad.append((cid,"ABSENT")); print(f"  {cid}  ABSENT FROM DATABASE"); continue
    fs=[f.lower() for f in fields(r)]
    hit=any(a in f or f==a for a in accept for f in fs)
    mark="ok " if hit else "MISMATCH"
    print(f"  [{mark}] {cid}  name={r.get('name','')!r:<18} abbr={r.get('abbreviation','')!r:<18} formula={r.get('formula','')!r}")
    if not hit: bad.append((cid,r.get('name','')))
if bad:
    print(f"\n  [FAIL] {len(bad)} id(s) unconfirmed: {bad}")
    print("         MEDIA NOT WRITTEN.")
    sys.exit(1)
print(f"\n  [OK]   all {len(SPEC)} ids confirmed against the ModelSEED record")

# ---- the build happens HERE, inside the check that gates it ----
def n(c): return rows[c].get('name') or rows[c].get('abbreviation') or c
SED=["cpd00001","cpd00067","cpd00009","cpd00013","cpd00048","cpd00239","cpd00011",
     "cpd00242","cpd00205","cpd00971","cpd00254","cpd00063","cpd00099","cpd10515",
     "cpd10516","cpd00030","cpd00034","cpd00058","cpd00149","cpd00244"]
CARB=["cpd00027","cpd00154","cpd00158","cpd00280"]
EXTRA=["cpd00029","cpd11640","cpd00047","cpd00116","cpd00441","cpd00098","cpd00166"]
AA=["cpd00023","cpd00033","cpd00035","cpd00039","cpd00041","cpd00051","cpd00053",
    "cpd00054","cpd00060","cpd00065","cpd00066","cpd00069","cpd00084","cpd00107",
    "cpd00119","cpd00129","cpd00132","cpd00156","cpd00161","cpd00322"]
AA=[a for a in AA if a in rows]
analysis=SED+CARB; fill=SED+CARB+EXTRA+AA
BANNED={"cpd00007","cpd00209","cpd00075","cpd00659","cpd00418"}
assert not (set(analysis)|set(fill)) & BANNED, "oxidant present"
assert not set(analysis) & set(EXTRA+AA), "currency/route/B12/AA in analysis medium"
assert set(analysis) < set(fill), "analysis medium is not a strict subset of fill medium"
print("  [OK]   assertions: no oxidant; analysis medium clean; strict subset")
def w(name,ids):
    p=os.path.join(out,f"{name}-medium.csv")
    with open(p,'w') as fh:
        fh.write("compounds,name,maxFlux\n")
        for c in sorted(set(ids)): fh.write(f"{c},{n(c)},10\n")
    print(f"  [OK]   {p}  ({len(set(ids))})")
w("ETO_FILL",fill); w("ETO_SED",analysis)
w("ETO_SED_NOSO4",[c for c in analysis if c!="cpd00048"])
p=os.path.join(out,"etoliko_media.tsv")
with open(p,'w') as fh:
    fh.write("medium\tdescription\tcompound\tname\n")
    for lab,ids in (("ETO_SED",analysis),("ETO_SED_NOSO4",[c for c in analysis if c!="cpd00048"])):
        for c in sorted(set(ids)): fh.write(f"{lab}\tEtoliko sediment, anoxic\t{c}\t{n(c)}\n")
print(f"  [OK]   {p}")
print("\n  ...    ETO_SED (analysis) — what the community is given:")
for c in sorted(set(analysis)): print(f"           {c}  {n(c)}")
print(f"\n  ...    present during FILL, absent during ANALYSIS ({len(set(fill)-set(analysis))}):")
for c in sorted(set(fill)-set(analysis)): print(f"           {c}  {n(c)}")
print("           ^ this is the list SMETANA must explain by cross-feeding.")
PY
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
