#!/usr/bin/env bash
# 41k -- remove O2 from the medium, and make its absence an ASSERTION rather than an
#        intention. Then emit the fill command, because feasibility can only be tested
#        after the fill, not before.
#
# 41j built ETO_SED by subtracting from gapseq's default base. That base is AEROBIC:
# it contains cpd00007 (O2), plus acetate, lactate and 12 amino acids. My exclusion list
# was written for cross-feeding currencies and never considered oxygen, so O2 survived
# into a medium for a permanently anoxic lagoon. smetana --anaerobic would have hidden
# it: the flag constrains the community run, but gapseq would already have gap-filled
# aerobic respiration into every model.
#
# 41j section D is also void as a control: YGCA01 is an aerobic Sporosarcina filled on an
# amino-acid-rich medium. It cannot grow anoxically without amino acids no matter what
# the medium says. A control whose failure is guaranteed by the organism tests nothing.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"
CULT="${CULT:-$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis}"; G="$CULT/gapseq_project"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41k  $(date -Is)"
MED=$(ls "$G"/*-medium.csv 2>/dev/null | head -1)
[ -z "$MED" ] && { no "no gapseq base medium"; exit 1; }

say "A. rebuild, with oxidants excluded by name"
python3 - "$MED" "$OUT" <<'PY'
import sys,csv,os
med,out=sys.argv[1],sys.argv[2]
rows=list(csv.DictReader(open(med)))
base={str(r['compounds']).split('[')[0].strip():(str(r.get('name','')).strip(),
      str(r.get('maxFlux','10')).strip()) for r in rows}

# 1. OXIDANTS -- the lagoon is permanently anoxic and euxinic. None of these may appear.
OXIDANTS={"cpd00007":"O2","cpd00209":"nitrate","cpd00075":"nitrite","cpd00659":"N2O",
          "cpd00418":"NO","cpd11574":"molybdate?no","cpd00528":"N2"}
OXIDANTS.pop("cpd11574")            # molybdate is a trace nutrient, not an oxidant
# 2. cross-feeding currencies -- supplying these erases the finding
CURRENCY={"cpd00029":"acetate","cpd11640":"H2","cpd00047":"formate","cpd00159":"lactate",
          "cpd00211":"butyrate","cpd00036":"succinate","cpd00020":"pyruvate",
          "cpd00141":"propionate","cpd00221":"D-lactate"}
# 3. the headline methylotrophic route -- must come only from another member
ROUTE={"cpd00116":"methanol","cpd00441":"trimethylamine","cpd00098":"choline",
       "cpd00540":"dimethylamine","cpd00443":"methylamine"}
# 4. step 21's dependency layer
DEPEND={"cpd00166":"cobalamin","cpd03424":"vitamin B12"}
AA={"cpd00023":"L-glutamate","cpd00033":"glycine","cpd00035":"L-alanine",
    "cpd00039":"L-lysine","cpd00041":"L-aspartate","cpd00051":"L-arginine",
    "cpd00053":"L-glutamine","cpd00054":"L-serine","cpd00060":"L-methionine",
    "cpd00065":"L-tryptophan","cpd00066":"L-phenylalanine","cpd00069":"L-tyrosine",
    "cpd00084":"L-cysteine","cpd00107":"L-leucine","cpd00119":"L-histidine",
    "cpd00129":"L-proline","cpd00132":"L-asparagine","cpd00156":"L-valine",
    "cpd00161":"L-threonine","cpd00322":"L-isoleucine"}
# 5. sugars that are not sediment carbon
NOTSED={"cpd00082":"D-fructose"}
GROUPS=[("OXIDANT -- anoxic lagoon",OXIDANTS),("cross-feeding currency",CURRENCY),
        ("methylotrophic route",ROUTE),("B12 dependency layer",DEPEND),
        ("free amino acid",AA),("non-sediment sugar",NOTSED)]
removed=[]
kept=dict(base)
for label,grp in GROUPS:
    for cid,nm in grp.items():
        if cid in kept:
            removed.append((cid,nm,label)); kept.pop(cid)
ADD={"cpd00048":("Sulfate","10"),"cpd00239":("H2S","10"),"cpd00011":("CO2","10"),
     "cpd00242":("HCO3","10"),"cpd10516":("Fe3+","10"),
     "cpd00027":("D-Glucose","10"),"cpd00154":("D-Xylose","10")}
print(f"  ...    gapseq base: {len(base)}   removed: {len(removed)}   added: {len(ADD)}")
print("  ...    removals, with the reason each one is removed:")
for cid,nm,label in sorted(removed,key=lambda x:(x[2],x[0])):
    print(f"           {cid:<10} {nm:<18} {label}")
final=dict(kept); final.update(ADD)

# ---- THE ASSERTION. Not an intention -- a check that can fail. ----
bad=[c for c in final if c in OXIDANTS]
if bad:
    print(f"  [FAIL] ANOXIC ASSERTION FAILED: {bad} present. Refusing to write.")
    sys.exit(1)
print("  [OK]   ANOXIC ASSERTION PASSED: no O2, nitrate, nitrite, NO or N2O in the medium")
for grp_label,grp in (("currency",CURRENCY),("route",ROUTE),("B12",DEPEND)):
    b=[c for c in final if c in grp]
    if b: print(f"  [FAIL] {grp_label} compound(s) leaked back in: {b}"); sys.exit(1)
print("  [OK]   no cross-feeding currency, methylotrophic substrate or B12 in the medium")

def write(name,items):
    p=os.path.join(out,f"{name}-medium.csv")
    with open(p,"w") as fh:
        fh.write("compounds,name,maxFlux\n")
        for cid,(nm,mf) in sorted(items.items()): fh.write(f"{cid},{nm},{mf}\n")
    print(f"  [OK]   {p}  ({len(items)} compounds)")
write("ETO_SED",final)
write("ETO_SED_NOSO4",{k:v for k,v in final.items() if k!="cpd00048"})
p=os.path.join(out,"etoliko_media.tsv")
with open(p,"w") as fh:
    fh.write("medium\tdescription\tcompound\tname\n")
    for nm_,items in (("ETO_SED",final),
                      ("ETO_SED_NOSO4",{k:v for k,v in final.items() if k!="cpd00048"})):
        for cid,(n2,_) in sorted(items.items()):
            fh.write(f"{nm_}\tEtoliko in-silico sediment, anoxic\t{cid}\t{n2}\n")
print(f"  [OK]   {p}")
print("\n  ...    final ETO_SED composition:")
for cid,(nm,mf) in sorted(final.items()): print(f"           {cid:<10} {nm}")
PY

say "B. re-check every medium file on disk for oxidants"
note "including any written by 41i or 41j, which may still contain O2"
for f in "$OUT"/*-medium.csv "$OUT"/etoliko_media.tsv; do
  [ -f "$f" ] || continue
  if grep -qE 'cpd00007|cpd00209|cpd00075' "$f"; then
    no "$(basename "$f") CONTAINS AN OXIDANT:"
    grep -nE 'cpd00007|cpd00209|cpd00075' "$f" | sed 's/^/           /'
  else
    ok "$(basename "$f") clean"
  fi
done

say "C. the fill command -- feasibility is testable only after this"
note "gapseq fill takes the medium. A model filled on ETO_SED can grow on ETO_SED;"
note "one filled on gapseq's default cannot, which is exactly what 41j section D saw."
cat <<'EOS'
           # one genome first, to prove the chain end-to-end:
           GS=$WS2/envs/gapseq_local
           M=$WS2/41_env_checks/ETO_SED-medium.csv
           F=<one MAG .fa from 07_mags/bins>
           $GS/gapseq find      -O -p all -b 200 -m Bacteria -K 3 $F
           $GS/gapseq find-transport -O -K 3 $F
           $GS/gapseq draft     -r <F>-all-Reactions.tbl -t <F>-Transporter.tbl -c $F -u 200 -l 100
           $GS/gapseq fill      -m <F>-draft.RDS -n $M -c <F>-rxnWeights.RDS -g <F>-rxnXgenes.RDS
           # then: does it grow on the medium it was filled on? That is the real test.
EOS
note "archaeal MAGs need -m Archaea in 'find', not Bacteria. Eight of the 29."
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
