#!/usr/bin/env bash
# 41j -- build the Etoliko environment BY SUBTRACTION from a medium known to work.
#
# 41i: gapseq's own 36-compound medium -> mro computes. My hand-built 21 -> nothing grows.
# Building a medium from first principles guesses what a model needs. Deriving it from
# one that demonstrably supports growth does not.
#
# Also settles MIP: mip_score returns None when minimal_medium cannot hit min_growth=0.1.
# 41h measured M2 at 0.069 unconstrained. It cannot reach 0.1 on ANY medium.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; S="$OUT/smoke_run"
CULT="${CULT:-$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis}"; G="$CULT/gapseq_project"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41j  $(date -Is)"
PY="$WS2/envs/smetana/bin/python"
MED=$(ls "$G"/*-medium.csv 2>/dev/null | head -1)

say "A. confirm the MIP explanation before moving on"
"$PY" - "$S/M1.xml" "$S/M2.xml" <<'PY' 2>&1 | grep -v UserWarning | grep -v '^  warn'
import sys, warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, FBA
g={}
for lab,p in (('M1',sys.argv[1]),('M2',sys.argv[2])):
    m=load_cbmodel(p,flavor='fbc2'); s=FBA(m); g[lab]=getattr(s,'fobj',None)
    print(f"  ...    {lab} maximum growth, unconstrained: {g[lab]:.4f}")
floor=0.1
bad=[k for k,v in g.items() if v is None or v<floor]
if bad:
    print(f"  [OK]   PREDICTION CONFIRMED: {', '.join(bad)} cannot reach min_growth={floor}")
    print("         mip_score therefore returns None for this pair on every medium.")
    print("         MIP=n/a is a property of THIS control pair, not of the pipeline.")
else:
    print(f"  [FAIL] both members exceed {floor} -- MIP=n/a needs another explanation")
PY

say "B. what is in gapseq's working medium that mine lacks?"
[ -z "$MED" ] && { no "no gapseq medium"; exit 1; }
python3 - "$MED" "$OUT" <<'PY'
import sys,csv,os
med,out=sys.argv[1],sys.argv[2]
rows=list(csv.DictReader(open(med)))
base={str(r['compounds']).split('[')[0].strip():str(r.get('name','')).strip() for r in rows}
mine={"cpd00001","cpd00067","cpd00009","cpd00013","cpd00048","cpd00239","cpd00011",
      "cpd00205","cpd00971","cpd00254","cpd00063","cpd00099","cpd10515","cpd00030",
      "cpd00034","cpd00058","cpd00149","cpd00244","cpd11574","cpd00027","cpd00154"}
missing={k:v for k,v in base.items() if k not in mine}
extra=sorted(mine-set(base))
print(f"  ...    gapseq medium: {len(base)} compounds; mine: {len(mine)}")
print(f"  ...    IN gapseq, ABSENT from mine ({len(missing)}) -- these are what a model expects:")
for k,v in sorted(missing.items()): print(f"           {k}  {v}")
print(f"  ...    in mine, not in gapseq's ({len(extra)}): {', '.join(extra)}")
open(os.path.join(out,'gapseq_base_medium.tsv'),'w').write(
    "cpd\tname\n"+"".join(f"{k}\t{v}\n" for k,v in sorted(base.items())))
print(f"  [OK]   wrote {out}/gapseq_base_medium.tsv")
PY

say "C. Etoliko medium by subtraction, and what gets removed"
python3 - "$MED" "$OUT" <<'PY'
import sys,csv,os
med,out=sys.argv[1],sys.argv[2]
rows=list(csv.DictReader(open(med)))
base={str(r['compounds']).split('[')[0].strip():(str(r.get('name','')).strip(),
      str(r.get('maxFlux','10')).strip()) for r in rows}
# The chapter's exclusions, by ModelSEED id. Anything here that is IN gapseq's base
# gets REMOVED and named, so the removal is visible rather than silent.
EXCLUDE={"cpd00116":"methanol","cpd00441":"trimethylamine","cpd00098":"choline",
         "cpd00166":"cobalamin/B12","cpd03424":"vitamin B12","cpd00029":"acetate",
         "cpd11640":"H2","cpd00047":"formate","cpd00159":"lactate","cpd00211":"butyrate",
         "cpd00036":"succinate","cpd00020":"pyruvate","cpd00023":"L-glutamate",
         "cpd00033":"glycine","cpd00035":"L-alanine","cpd00039":"L-lysine",
         "cpd00041":"L-aspartate","cpd00051":"L-arginine","cpd00054":"L-serine",
         "cpd00060":"L-methionine","cpd00065":"L-tryptophan","cpd00066":"L-phenylalanine",
         "cpd00069":"L-tyrosine","cpd00084":"L-cysteine","cpd00107":"L-leucine",
         "cpd00119":"L-histidine","cpd00129":"L-proline","cpd00156":"L-valine",
         "cpd00161":"L-threonine","cpd00322":"L-isoleucine","cpd00053":"L-glutamine",
         "cpd00132":"L-asparagine"}
ADD={"cpd00048":("sulfate","10"),"cpd00239":("h2s","10"),
     "cpd00011":("co2","10"),"cpd00027":("d-glucose","10"),"cpd00154":("xylose","10")}
removed={k:EXCLUDE[k] for k in base if k in EXCLUDE}
kept={k:v for k,v in base.items() if k not in EXCLUDE}
print(f"  ...    removed from gapseq's base ({len(removed)}):")
for k,v in sorted(removed.items()): print(f"           {k}  {v}")
if not removed: print("           (none -- gapseq's base contains no excluded compound)")
final=dict(kept); final.update(ADD)
print(f"  ...    added ({len(ADD)}): "+", ".join(f"{k} ({v[0]})" for k,v in ADD.items()))
print(f"  ...    ETO_SED final: {len(final)} compounds")
def write(name,items,desc):
    p=os.path.join(out,f"{name}-medium.csv")
    with open(p,"w") as fh:
        fh.write("compounds,name,maxFlux\n")
        for cid,(nm,mf) in sorted(items.items()): fh.write(f"{cid},{nm},{mf}\n")
    print(f"  [OK]   {p}  ({len(items)})")
    return p
write("ETO_SED",final,"primary")
noso4={k:v for k,v in final.items() if k!="cpd00048"}
write("ETO_SED_NOSO4",noso4,"no sulfate")
p=os.path.join(out,"etoliko_media.tsv")
with open(p,"w") as fh:
    fh.write("medium\tdescription\tcompound\tname\n")
    for med_name,items in (("ETO_SED",final),("ETO_SED_NOSO4",noso4)):
        for cid,(nm,_) in sorted(items.items()):
            fh.write(f"{med_name}\tEtoliko in-silico sediment\t{cid}\t{nm}\n")
print(f"  [OK]   {p}")
PY

say "D. does anything actually grow on it?"
note "M1/M2 are Sporosarcina isolates, not our MAGs -- but if THEY cannot grow, the"
note "medium is too poor for anything, and that is worth knowing before the fill run."
"$PY" - "$S/M1.xml" "$S/M2.xml" "$OUT/ETO_SED-medium.csv" <<'PY' 2>&1 | grep -v UserWarning | grep -v '^  warn'
import sys,csv,warnings; warnings.filterwarnings('ignore')
from reframed import load_cbmodel, FBA
from reframed.core.environment import Environment
cpds=[r['compounds'].split('[')[0].strip() for r in csv.DictReader(open(sys.argv[3]))]
for lab,p in (('M1',sys.argv[1]),('M2',sys.argv[2])):
    m=load_cbmodel(p,flavor='fbc2')
    env=Environment.from_compounds(cpds, fmt_func=lambda x: f"R_EX_{x}_e")
    present=[k for k in env if k in m.reactions]
    print(f"  ...    {lab}: {len(present)}/{len(env)} medium compounds exist as exchanges")
    try:
        env.apply(m, inplace=True, warning=False)
        s=FBA(m); v=getattr(s,'fobj',None)
        print(f"  ...    {lab} growth on ETO_SED: {v!r}")
        if v is None or v<1e-9:
            print(f"  [FAIL] {lab} does NOT grow on ETO_SED -- medium still too poor")
        else:
            print(f"  [OK]   {lab} grows on ETO_SED")
    except Exception as e:
        print(f"  [FAIL] {lab}: {type(e).__name__}: {e}")
print("  ...    NOTE: these organisms were never gap-filled on ETO_SED. Growth here is")
print("         encouraging; failure would be decisive. Neither is proof for our MAGs.")
PY
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
