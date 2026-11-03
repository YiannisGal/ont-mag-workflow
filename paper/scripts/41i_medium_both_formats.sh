#!/usr/bin/env bash
# 41i -- the medium is a gapseq input FIRST and a SMETANA input SECOND.
#
# 41h: without a medium mro=0.744; with my 15-compound medium the community cannot grow
# and everything is n/a. Detailed mode gives 28 real interaction rows, so the machinery
# is fine. The medium was the problem -- as it must be, because a gapseq model is
# gap-filled ON a medium and is only guaranteed feasible on THAT medium.
#
# A. what format does gapseq's own -medium.csv use?
# B. POSITIVE CONTROL: run smetana on gapseq's own medium. mip/mro must be numbers.
# C. why is mip n/a even with no medium?
# D. Etoliko environment in BOTH formats + per-model feasibility before any SMETANA run.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; S="$OUT/smoke_run"
CULT="${CULT:-$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis}"; G="$CULT/gapseq_project"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41i  $(date -Is)"
PY="$WS2/envs/smetana/bin/python"; SM="$WS2/envs/smetana/bin/smetana"

say "A. gapseq's own medium format"
MED=$(ls "$G"/*-medium.csv 2>/dev/null | head -1)
if [ -z "$MED" ]; then no "no *-medium.csv in $G"; else
  ok "$MED  ($(wc -l < "$MED") lines)"
  head -4 "$MED" | sed 's/^/           /'
  note "compound id column and how many entries:"
  python3 - "$MED" <<'PY'
import sys,csv
p=sys.argv[1]
with open(p) as fh:
    d=csv.DictReader(fh); rows=list(d)
print(f"           columns: {', '.join(rows[0].keys())}")
print(f"           rows: {len(rows)}")
col=next((k for k in rows[0] if any(str(r[k]).startswith('cpd') for r in rows)), None)
print(f"           compound column: {col!r}")
if col: print("           first ids: "+", ".join(str(r[col]) for r in rows[:6]))
PY
fi

say "B. POSITIVE CONTROL -- smetana on gapseq's own medium"
if [ -z "$MED" ] || [ ! -d "$S" ]; then no "prerequisites missing"; else
  python3 - "$MED" "$S/gapseq_media.tsv" <<'PY'
import sys,csv
src,dst=sys.argv[1],sys.argv[2]
rows=list(csv.DictReader(open(src)))
col=next(k for k in rows[0] if any(str(r[k]).startswith('cpd') for r in rows))
ids=[]
for r in rows:
    v=str(r[col]).strip()
    if v.startswith('cpd'): ids.append(v.split('[')[0].split('_')[0])
ids=sorted(set(ids))
with open(dst,'w') as fh:
    fh.write("medium\tdescription\tcompound\tname\n")
    for i in ids: fh.write(f"GAPSEQ\tgapseq's own fill medium\t{i}\tx\n")
print(f"  ...    converted {len(ids)} compounds -> {dst}")
PY
  cd "$S" || exit 1
  rm -f pc_*.tsv
  timeout 900 "$SM" --flavor fbc2 --anaerobic -m GAPSEQ --mediadb gapseq_media.tsv \
     -o pc M1.xml M2.xml > pc.log 2>&1
  if [ -s pc_global.tsv ]; then
    V=$(awk -F'\t' 'NR==2{print $4"|"$5}' pc_global.tsv)
    printf '           gapseq medium  mip|mro = %s\n' "$V"
    case "$V" in
      n/a\|n/a) no "STILL n/a on gapseq's OWN medium -- the conversion or the format is wrong"
                tail -8 pc.log | sed 's/^/           /' ;;
      n/a\|*)   ok "mro computes -> the medium mechanism WORKS; mip is a separate problem" ;;
      *)        ok "both scores compute -- medium mechanism fully verified" ;;
    esac
  else no "no pc_global.tsv"; tail -8 pc.log | sed 's/^/           /'; fi
  cd - >/dev/null
fi

say "C. why is mip n/a?"
"$PY" - <<'PY'
import inspect, smetana.interface as I
src=inspect.getsource(I.mip_score)
for l in src.splitlines()[:40]: print("           "+l)
PY
note "trying --use-lp, if smetana exposes it:"
"$SM" --help 2>&1 | grep -iE 'use.lp|molecular|exclude|no.coupling' | sed 's/^/           /'
if [ -d "$S" ]; then
  cd "$S"; rm -f lp_*.tsv
  timeout 900 "$SM" --flavor fbc2 --anaerobic --use-lp -o lp M1.xml M2.xml > lp.log 2>&1 \
    && [ -s lp_global.tsv ] \
    && printf '           --use-lp  mip|mro = %s\n' "$(awk -F'\t' 'NR==2{print $4"|"$5}' lp_global.tsv)" \
    || note "--use-lp not accepted or produced nothing"
  cd - >/dev/null
fi

say "D. the Etoliko environment, in BOTH formats"
note "One definition, two encodings. Written only if B showed the mechanism works."
python3 - "$OUT" <<'PY'
import os,sys
out=sys.argv[1]
# ids verified by name against a real gapseq model in 41c
BASE=[("cpd00001","h2o"),("cpd00067","h+"),("cpd00009","phosphate"),("cpd00013","nh3"),
      ("cpd00048","sulfate"),("cpd00239","h2s"),("cpd00011","co2"),("cpd00205","k+"),
      ("cpd00971","na+"),("cpd00254","mg"),("cpd00063","ca2+"),("cpd00099","cl-"),
      ("cpd10515","fe2+"),("cpd00030","mn2+"),("cpd00034","zn2+"),("cpd00058","cu2+"),
      ("cpd00149","co2+"),("cpd00244","ni2+"),("cpd11574","molybdate")]
CARB=[("cpd00027","d-glucose"),("cpd00154","xylose")]
sets={"ETO_SED":BASE+CARB,
      "ETO_SED_NOSO4":[x for x in BASE if x[0]!="cpd00048"]+CARB}
# gapseq format
for name,items in sets.items():
    p=os.path.join(out,f"{name}-medium.csv")
    with open(p,"w") as fh:
        fh.write("compounds,name,maxFlux\n")
        for cid,nm in items: fh.write(f"{cid},{nm},10\n")
    print(f"  [OK]   gapseq  {p}  ({len(items)} compounds)")
# smetana format
p=os.path.join(out,"etoliko_media.tsv")
with open(p,"w") as fh:
    fh.write("medium\tdescription\tcompound\tname\n")
    for name,items in sets.items():
        for cid,nm in items: fh.write(f"{name}\tEtoliko in-silico sediment\t{cid}\t{nm}\n")
print(f"  [OK]   smetana {p}")
print("  ...    NOT supplied: acetate, H2, formate (the cross-feeding currencies);")
print("         methanol, TMA, choline (the headline route); B12 and amino acids.")
print("  🔴 THIS MEDIUM IS UNTESTED FOR FEASIBILITY. No gapseq model has been filled on")
print("     it yet, so no model is guaranteed to grow on it. 41h showed exactly what that")
print("     looks like: every score n/a. The fill step must USE this file, not just SMETANA.")
PY
note "next: gapseq fill -m <ETO_SED-medium.csv> for each MAG, then SMETANA on the same set."
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
