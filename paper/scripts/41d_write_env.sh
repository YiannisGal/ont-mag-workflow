#!/usr/bin/env bash
# 41d -- settle the exchange-reaction convention, then WRITE the environment file.
#
# 41c section B called system python3, which has no reframed. reframed lives in the
# SMETANA venv. The open question it could not answer: gapseq writes R_EX_cpd00013_e0,
# and SMETANA's media loader builds an id from the mediadb 'compound' column by some
# internal rule. If that rule yields R_EX_cpd00013_e, every compound fails to load,
# SMETANA runs, and reports no cross-feeding. Silent. So: read the rule from the
# installed source, then verify empirically against the model.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; mkdir -p "$OUT"
CULT="${CULT:-$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis}"
M="$CULT/gapseq_project/YGCA01_Clean.xml"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41d  $(date -Is)"

PY=""
for p in "$WS2/envs/smetana/bin/python" "$WS2/envs/smetana/bin/python3"; do
  [ -x "$p" ] && { PY="$p"; break; }; done
[ -z "$PY" ] && { no "smetana venv python not found under $WS2/envs/smetana/bin"; exit 1; }
ok "interpreter: $PY"
"$PY" -c 'import reframed,smetana;print("  ...    reframed",reframed.__version__ if hasattr(reframed,"__version__") else "?","| smetana ok")' 2>&1 | tail -2

say "A. what id does SMETANA build from a mediadb compound?"
SP=$("$PY" -c 'import smetana,os;print(os.path.dirname(smetana.__file__))' 2>/dev/null)
note "smetana package at $SP"
grep -rn "fmt_func\|R_EX_\|EX_{\|_e'" "$SP" 2>/dev/null | head -12 | sed 's/^/           /'
RF=$("$PY" -c 'import reframed,os;print(os.path.dirname(reframed.__file__))' 2>/dev/null)
grep -rn "def from_compounds" -A12 "$RF/core/environment.py" 2>/dev/null | head -20 | sed 's/^/           /'

say "B. empirical: what does the model actually expose, and does it match?"
"$PY" - "$M" "$OUT" <<'PY'
import sys, os
m,out=sys.argv[1],sys.argv[2]
from reframed import load_cbmodel
from reframed.core.environment import Environment
mod=None
for fl in ('fbc2','cobra'):
    try:
        mod=load_cbmodel(m,flavor=fl); print(f"  [OK]   loaded with flavor={fl!r}"); break
    except Exception as e: print(f"  ...    flavor={fl!r} failed: {type(e).__name__}")
if mod is None:
    print("  [FAIL] reframed cannot load a gapseq model at all -- stop here"); sys.exit(0)
ex=sorted(r for r in mod.reactions if 'EX_' in r)
print(f"  ...    {len(mod.reactions)} reactions, {len(ex)} exchanges. First five:")
for r in ex[:5]: print(f"           {r}")
probe=['cpd00013','cpd00009','cpd00048','cpd00001']
try:
    env=Environment.from_compounds(probe)
    built=sorted(env.keys()); print("  ...    Environment.from_compounds() builds:")
    for b in built: print(f"           {b}")
    hit=[b for b in built if b in mod.reactions]
    if len(hit)==len(built):
        print(f"  [OK]   all {len(built)} generated ids exist in the model -- plain cpd id works")
        print("  CONVENTION: mediadb compound column = bare cpd id")
        open(f"{out}/EXCHANGE_CONVENTION.txt","w").write("bare\n")
    else:
        print(f"  [FAIL] {len(built)-len(hit)} of {len(built)} generated ids are NOT in the model.")
        # find what suffix WOULD work
        import re
        sample=[r for r in ex if 'cpd00013' in r or 'cpd00009' in r]
        print("  ...    the model spells them:", ", ".join(sample[:3]) or "(absent)")
        cand=set()
        for r in ex:
            mm=re.match(r'^(R_)?EX_(cpd\d{5})(.*)$',r)
            if mm: cand.add(mm.group(3))
        print(f"  ...    observed suffixes on model exchanges: {sorted(cand)}")
        for suf in sorted(cand):
            test=[f"cpd00013{suf}"]
            try:
                e2=Environment.from_compounds(test); k=list(e2.keys())[0]
                mark="MATCH" if k in mod.reactions else "no"
                print(f"           compound='cpd00013{suf}' -> {k}   {mark}")
                if k in mod.reactions:
                    open(f"{out}/EXCHANGE_CONVENTION.txt","w").write(suf+"\n")
            except Exception as e2e: print(f"           compound='cpd00013{suf}' -> error {e2e}")
except Exception as e:
    print(f"  [FAIL] Environment.from_compounds raised {type(e).__name__}: {e}")
PY

say "C. write the environment file"
CONV="$OUT/EXCHANGE_CONVENTION.txt"
if [ ! -s "$CONV" ]; then
  no "convention undetermined -- NOT writing the environment file"
  no "a mediadb that does not match the models yields 'no cross-feeding' silently"
else
  SUF=$(cat "$CONV"); [ "$SUF" = "bare" ] && SUF=""
  ok "compound column suffix: '${SUF:-<none>}'"
  python3 - "$OUT" "$SUF" <<'PY'
import sys,os
out,suf=sys.argv[1],(sys.argv[2] if len(sys.argv)>2 else "")
if suf=="bare": suf=""
# VERIFIED ids only (41c matched each against the model's own name).
# The medium is a WHITELIST: methanol, TMA, choline, B12 and amino acids are
# excluded by absence, so their ids never have to be correct.
BASE=[("cpd00001","h2o"),("cpd00067","h+"),("cpd00009","phosphate"),
      ("cpd00013","ammonium (nh3)"),("cpd00048","sulfate"),("cpd00239","sulfide h2s"),
      ("cpd00011","co2"),("cpd00205","k+"),("cpd00971","na+"),("cpd00254","mg"),
      ("cpd00063","ca2+"),("cpd00099","cl-"),("cpd10515","fe2+"),("cpd00030","mn2+"),
      ("cpd00034","zn2+"),("cpd00058","cu2+"),("cpd00149","co2+"),("cpd00244","ni2+"),
      ("cpd11574","molybdate")]
CARB=[("cpd00027","d-glucose"),("cpd00154","xylose")]        # verified
UNV =[("cpd00158","cellobiose"),("cpd00280","galacturonate")] # id NOT yet verified
media={
 "ETO_SED":      BASE+CARB,                       # run 1, primary
 "ETO_SED_RICH": BASE+CARB+[("cpd00166","B12?"),], # run 2 -- see note below
 "ETO_SED_NOSO4":[x for x in BASE if x[0]!="cpd00048"]+CARB,  # run 3, the test
}
p=os.path.join(out,"etoliko_media.tsv")
with open(p,"w") as fh:
    fh.write("medium\tdescription\tcompound\tname\n")
    for med,items in media.items():
        d={"ETO_SED":"in-silico Etoliko sediment: anoxic, sulfate non-limiting, plant-polymer carbon",
           "ETO_SED_RICH":"sensitivity: sediment + B12 (id UNVERIFIED - do not run until checked)",
           "ETO_SED_NOSO4":"contrast: sediment WITHOUT sulfate"}[med]
        for cid,nm in items: fh.write(f"{med}\t{d}\t{cid}{suf}\t{nm}\n")
print(f"  [OK]   wrote {p}")
for med,items in media.items(): print(f"           {med:<16} {len(items)} compounds")
print("  ...    NOT supplied, deliberately: acetate, H2, formate (cross-feeding currencies);")
print("         methanol, TMA, choline (the headline methylotrophic route);")
print("         B12 and free amino acids (step 21's dependency layer).")
print("  [WARN] ETO_SED_RICH uses cpd00166 for B12, which 41c could NOT verify.")
print("         Verify it against a real MAG model before running run 2.")
print(f"  ...    also unverified and therefore OMITTED from every medium: "
      + ", ".join(f"{c} ({n})" for c,n in UNV))
print("         Carbon is currently glucose + xylose only. Adding cellobiose and")
print("         galacturonate needs their ids confirmed -- they are the cellulose and")
print("         pectin arms, so this is a real gap, not a cosmetic one.")
PY
fi

say "D. curated-HMM resources for the mcr question"
note "product strings gave mcr=0 for 8 obligate methanogens. Bakta is a BACTERIAL"
note "annotator; this is the second archaeal under-annotation after 39c's secretome."
for d in "$WS2"/*kofam* "$WS2"/*/*kofam* "$HOME"/kofam* ${SCRATCH_ROOT}/*kofam*; do
  [ -e "$d" ] && echo "           found: $d"; done 2>/dev/null
command -v hmmsearch >/dev/null && ok "hmmsearch present: $(command -v hmmsearch)" || no "hmmsearch absent"
command -v exec_annotation >/dev/null && ok "kofamscan present" || note "kofamscan not on PATH"
find "$WS2" -maxdepth 3 \( -name '*.hmm' -o -name 'ko_list*' -o -name 'profiles' -type d \) 2>/dev/null | head -8 | sed 's/^/           found: /'
note "mcrA is KEGG K00399 (mcrA), K00401 (mcrB), K00402 (mcrG); TIGR03256 for mcrA."
note "8 archaeal proteomes x 3 profiles is seconds of hmmsearch."
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
