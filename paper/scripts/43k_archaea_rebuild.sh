#!/usr/bin/env bash
# 43k_archaea_rebuild.sh -- build MTP4, and make both archaea methanogens again.
#
# RUN IN tmux. gapseq find on a 4.7 Mb genome is hours. The script is RESUMABLE: every
# stage skips if its output already exists, so re-running after a drop costs nothing.
#     tmux new -s 43k
#     bash scripts/43k_archaea_rebuild.sh 2>&1 | tee logs/43k.log
#
# WHY
#   43j resolved the crosswalk:
#       ETO_LIQ_MAG_003 -> LIQ_das_vamb_003   99.84% complete, 4.72 Mb, 16S+23S+5S
#       ETO_LIQ_MAG_011 -> LIQ_das_vamb_017   90.48% complete, 2.91 Mb, no rRNA
#   MTP4 -- the Methanosarcinaceae at 0.6170 sediment breadth that Part 1's coexistence
#   argument rests on -- is LIQ_das_vamb_003, and it is NOT in the in-situ ten. The ten
#   were selected in step 41, before MTP4 was identified in 41b. We modelled
#   LIQ_das_vamb_017 (Methanofollis) throughout, and it is the worse genome of the two.
#
#   43j also showed the ETO_SED-filled LIQ_das_vamb_017 takes up glucose, cellobiose and
#   galacturonate at maximum rate and secretes H2 67.7, acetate 50.3, formaldehyde 27.1.
#   It is a sugar fermenter. Gapfilling on a medium the organism cannot inhabit does not
#   produce a model of that organism.
#
# WHAT THIS DOES
#   1. builds ETO_METH-medium.csv in gapseq's own format, BY SUBTRACTION from the
#      existing ETO_SED-medium.csv, so the format cannot drift from what gapseq accepts
#   2. builds LIQ_das_vamb_003 (MTP4) from genome to model, --biomass Archaea
#   3. refills LIQ_das_vamb_017 on ETO_METH from its existing draft
#   4. GATES on three assertions that can each fail:
#         grows on ETO_METH  /  methane exchange carries flux  /  NO sugar uptake
#      The third is the one 43b needed and did not have.
#
# It does not touch the nine bacterial models, workspace_1, or any shared conda env.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43k_archaea"; mkdir -p "$OUT"
export TMPDIR="${SCRATCH_ROOT}/tmp43k"; mkdir -p "$TMPDIR"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"
MDB="$WS2/41_env_checks/etoliko_media.tsv"
SEDCSV="$WS2/41_env_checks/ETO_SED-medium.csv"
METHCSV="$WS2/41_env_checks/ETO_METH-medium.csv"
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017
echo "43k  $(date -Is)"; echo "out: $OUT"; echo "TMPDIR: $TMPDIR (${#TMPDIR} chars)"

# ---------------------------------------------------------------- A. locate
say "A. the paths, PINNED -- discovery cost more than it saved"
# gapseq is a local checkout with the executable at the REPO ROOT, not under bin/.
# My 43k v1 locator required */bin/* and found nothing. Both paths are now fixed,
# recovered from logs/42*.log, and only CHECKED here rather than searched for.
GAPSEQ="$WS2/envs/gapseq_local/gapseq"
ENVDIR="$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis/gapseq_project/.pixi/envs/default"
if [ ! -x "$GAPSEQ" ]; then bad "not executable: $GAPSEQ"; echo "STEP_STATUS: FAILED"; exit 1; fi
if [ ! -d "$ENVDIR" ]; then bad "no pixi env at $ENVDIR"; echo "STEP_STATUS: FAILED"; exit 1; fi
case "$ENVDIR" in
  *damaged*) bad "refusing to use a .damaged env -- its references were overwritten by -U"
             echo "STEP_STATUS: FAILED"; exit 1;;
esac
ok "gapseq : $GAPSEQ"
ok "deps   : $ENVDIR"
note "that env is READ ONLY -- binaries only, never written (the -U incident)"
export PATH="$ENVDIR/bin:$PATH"
note "dependency check, AFTER PATH is set (42's mistake was checking before):"
MISS=0
for t in Rscript blastp hmmsearch bedtools exonerate barrnap; do
  if command -v "$t" >/dev/null 2>&1; then printf '     [OK]   %-10s %s\n' "$t" "$(command -v "$t")"
  else printf '     [--]   %-10s MISSING\n' "$t"; MISS=$((MISS+1)); fi
done
[ "$MISS" -gt 0 ] && note "$MISS dependency/ies missing; gapseq may still run if they are optional"
"$GAPSEQ" -v 2>&1 | head -3 | sed 's/^/     /' || true

note "genome FASTA for $MTP4:"
FA=""
while IFS= read -r c; do [ -s "$c" ] && { FA="$c"; break; }; done < <(
  find "$WS2" -type f \( -name "${MTP4}.fa" -o -name "${MTP4}.fna" -o -name "${MTP4}.fasta" \) 2>/dev/null
  find "$WS2" -type f \( -name "LIQ_vamb_003.fa" -o -name "LIQ_vamb_003.fna" \) 2>/dev/null )
[ -n "$FA" ] || { bad "no genome FASTA for $MTP4"; echo "STEP_STATUS: FAILED"; exit 1; }
ok "FASTA: $FA  ($(grep -c '^>' "$FA") contigs, $(grep -v '^>' "$FA" | tr -d '\n' | wc -c) bp)"

note "existing draft for $HYD (to refill without rebuilding):"
HYDRDS=$(find "$WS2/42_gems" -name "*${HYD}*draft*.RDS" -o -name "*${HYD}*.RDS" 2>/dev/null | head -1)
[ -n "$HYDRDS" ] && ok "draft RDS: $HYDRDS" || no "no draft RDS for $HYD -- it will be rebuilt too"

note "taxonomy of $MTP4, from whatever table holds it:"
grep -rhiE "vamb_003|ETO_LIQ_MAG_003" "$WS2"/*/*gtdb*.tsv "$WS2"/*/*taxonom*.tsv "$WS2"/*/*classif*.tsv 2>/dev/null | head -4 | cut -c1-220 | sed 's/^/     /' || true

# ---------------------------------------------------------------- B. medium
say "B. ETO_METH-medium.csv, built BY SUBTRACTION from the file gapseq already accepts"
if [ -s "$METHCSV" ]; then
  ok "already present: $METHCSV ($(($(wc -l < "$METHCSV")-1)) compounds)"
else
  [ -s "$SEDCSV" ] || { bad "no $SEDCSV to derive from"; echo "STEP_STATUS: FAILED"; exit 1; }
  "$PY" - "$SEDCSV" "$METHCSV" <<'PYEOF'
import sys, csv
src, dst = sys.argv[1], sys.argv[2]
SUGARS = {'cpd00027':'D-Glucose','cpd00154':'Xylose','cpd00158':'CELB','cpd00280':'D-Galacturonate'}
ADD    = {'cpd11640':'H2','cpd00011':'CO2'}
rows=[];
with open(src, newline='') as fh:
    rd = csv.DictReader(fh); fields = rd.fieldnames; rows = list(rd)
print("   source columns: %s" % fields)
idcol = next((c for c in fields if c.lower() in ('compounds','compound','id','cpd')), fields[0])
namecol = next((c for c in fields if c.lower() in ('name','description')), None)
kept, removed = [], []
for r in rows:
    (removed if r[idcol].strip() in SUGARS else kept).append(r)
print("   removed %d sugar row(s):" % len(removed))
for r in removed: print("     %-12s %s" % (r[idcol], r.get(namecol,'')))
have = {r[idcol].strip() for r in kept}
added=[]
for cpd, nm in ADD.items():
    if cpd in have: continue
    tmpl = dict(kept[0])                       # copy every column so the format matches exactly
    tmpl[idcol] = cpd
    if namecol: tmpl[namecol] = nm
    kept.append(tmpl); added.append((cpd, nm))
print("   added %d row(s): %s" % (len(added), ', '.join('%s %s' % a for a in added)))
with open(dst,'w',newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
    for r in kept: w.writerow(r)
print("   wrote %s with %d compounds" % (dst, len(kept)))
ids = {r[idcol].strip() for r in kept}
assert not (ids & set(SUGARS)),      "ASSERTION FAILED: a sugar survived into ETO_METH"
assert 'cpd00007' not in ids,        "ASSERTION FAILED: O2 in an anoxic medium"
assert 'cpd11640' in ids,            "ASSERTION FAILED: no H2 in a methanogenic medium"
print("   assertions passed: no sugars, no O2, H2 present")
PYEOF
  [ $? -eq 0 ] || { bad "medium build failed its own assertions"; echo "STEP_STATUS: FAILED"; exit 1; }
fi
note "and the same medium in the SMETANA media db (added by 43j):"
awk -F'\t' 'NR==1{next} $1=="ETO_METH"{n++} END{printf "     ETO_METH rows in media db: %d\n", n+0}' "$MDB"

# ---------------------------------------------------------------- C. flags
say "C. verify gapseq's flags PER SUBCOMMAND before passing any (Rule 68)"
for sub in find find-transport draft fill; do
  echo "   --- gapseq $sub ---"
  H=$("$GAPSEQ" $sub -h 2>&1 | head -40)
  if echo "$H" | grep -qiE 'error|traceback|unexpected'; then
    no "  -h did not return help text for '$sub'; flags will be used conservatively"
  fi
  echo "$H" | grep -oE '^\s+-[A-Za-z]\b[^ ]*' | tr -d ' ' | tr '\n' ' ' | sed 's/^/     flags: /'; echo
done

# ---------------------------------------------------------------- D. build MTP4
say "D. build $MTP4 (MTP4) from genome to model"
MD="$OUT/$MTP4"; mkdir -p "$MD"
cd "$MD"
B="$MTP4"
if [ -s "$MD/$B-all-Pathways.tbl" ] || [ -s "$MD/$B-Pathways.tbl" ]; then
  ok "pathways already found -- skipping gapseq find"
else
  note "gapseq find  (this is the long one; started $(date -Is))"
  "$GAPSEQ" find -p all -b 200 -m Archaea -t Archaea "$FA" 2>&1 | tail -5 | sed 's/^/     /' || \
  "$GAPSEQ" find -p all -b 200 "$FA" 2>&1 | tail -5 | sed 's/^/     /'
fi
if ls "$MD"/*-Transporter.tbl >/dev/null 2>&1; then
  ok "transporters already found -- skipping"
else
  note "gapseq find-transport"
  "$GAPSEQ" find-transport "$FA" 2>&1 | tail -5 | sed 's/^/     /'
fi
if ls "$MD"/*-draft.RDS >/dev/null 2>&1; then
  ok "draft already present -- skipping"
else
  note "gapseq draft --biomass Archaea   (42d: MAGs without 16S break Gram inference)"
  RXN=$(ls "$MD"/*-all-Reactions.tbl "$MD"/*-Reactions.tbl 2>/dev/null | head -1)
  TR=$(ls "$MD"/*-Transporter.tbl 2>/dev/null | head -1)
  "$GAPSEQ" draft -r "$RXN" -t "$TR" -p "$(ls "$MD"/*-all-Pathways.tbl "$MD"/*-Pathways.tbl 2>/dev/null | head -1)" \
      -c "$FA" -u 200 -l 100 --biomass Archaea 2>&1 | tail -8 | sed 's/^/     /'
fi
DRAFT=$(ls "$MD"/*-draft.RDS 2>/dev/null | head -1)
[ -n "$DRAFT" ] || { bad "no draft RDS produced for $MTP4"; echo "STEP_STATUS: FAILED"; exit 1; }
ok "draft: $(basename "$DRAFT")"
if ls "$MD"/*.xml >/dev/null 2>&1 && ! ls "$MD"/*draft*.xml >/dev/null 2>&1; then
  ok "filled model already present -- skipping fill"
else
  note "gapseq fill -n ETO_METH"
  WT=$(ls "$MD"/*-rxnWeights.RDS 2>/dev/null | head -1)
  XG=$(ls "$MD"/*-rxnXgenes.RDS 2>/dev/null | head -1)
  "$GAPSEQ" fill -m "$DRAFT" -n "$METHCSV" -c "$WT" -g "$XG" 2>&1 | tail -8 | sed 's/^/     /'
fi
MOD=$(ls "$MD"/*.xml 2>/dev/null | grep -v draft | head -1)
[ -n "$MOD" ] && ok "model: $(basename "$MOD")" || { bad "no filled model for $MTP4"; echo "STEP_STATUS: FAILED"; exit 1; }

# ---------------------------------------------------------------- E. refill hyd
say "E. refill $HYD on ETO_METH"
HD="$OUT/$HYD"; mkdir -p "$HD"
if ls "$HD"/*.xml 2>/dev/null | grep -qv draft; then
  ok "already refilled -- skipping"
elif [ -n "$HYDRDS" ]; then
  cd "$HD"
  SRC=$(dirname "$HYDRDS")
  WT=$(ls "$SRC"/*-rxnWeights.RDS 2>/dev/null | head -1)
  XG=$(ls "$SRC"/*-rxnXgenes.RDS 2>/dev/null | head -1)
  note "gapseq fill -m $(basename "$HYDRDS") -n ETO_METH"
  "$GAPSEQ" fill -m "$HYDRDS" -n "$METHCSV" ${WT:+-c "$WT"} ${XG:+-g "$XG"} 2>&1 | tail -8 | sed 's/^/     /'
else
  no "no draft for $HYD -- skipping; the community will run with MTP4 only"
fi

# ---------------------------------------------------------------- F. gate
say "F. GATE -- three assertions, each able to fail"
mkdir -p "$OUT/models"
for d in "$MD" "$HD"; do
  m=$(ls "$d"/*.xml 2>/dev/null | grep -v draft | head -1)
  [ -n "$m" ] && cp "$m" "$OUT/models/$(basename "$d").xml"
done
ls -l "$OUT/models" | sed 's/^/     /'
"$PY" - "$OUT/models" "$METHCSV" "$SEDCSV" <<'PYEOF'
import sys, os, csv, glob, warnings, re
from reframed import load_cbmodel, Environment, FBA
mdir, methcsv, sedcsv = sys.argv[1:4]
def cpds(path):
    out=[]
    with open(path, newline='') as fh:
        rd=csv.DictReader(fh)
        idc=next((c for c in rd.fieldnames if c.lower() in ('compounds','compound','id','cpd')), rd.fieldnames[0])
        for r in rd: out.append(r[idc].strip())
    return out
meth, sed = cpds(methcsv), cpds(sedcsv)
SUGAR={'cpd00027','cpd00154','cpd00158','cpd00280'}
fails=0
for p in sorted(glob.glob(os.path.join(mdir,'*.xml'))):
    n=os.path.basename(p).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m=load_cbmodel(p, flavor='fbc2')
    print("\n   %s  (%d reactions)" % (n, len(m.reactions)))
    mm=m.copy()
    Environment.from_compounds(meth, fmt_func=lambda x:"R_EX_{}_e".format(x)).apply(mm,inplace=True,warning=False)
    s=FBA(mm); g=s.fobj if s and s.fobj else 0.0
    print("     1. grows on ETO_METH        : %.4f  %s" % (g, "OK" if g>1e-6 else "FAIL"))
    if g<=1e-6: fails+=1; continue
    ch4=[r for r,rx in mm.reactions.items()
         if r in set(mm.get_exchange_reactions())
         and any(re.match(r'^methane', (mm.metabolites[k].name or ''), re.I) for k in rx.stoichiometry)]
    f=max((s.values.get(r,0.0) for r in ch4), default=0.0)
    print("     2. methane exchange flux    : %.4f  %s" % (f, "OK" if f>1e-6 else "FAIL -- not methanogenic"))
    if f<=1e-6: fails+=1
    up=[(r,v) for r,v in s.values.items()
        if v<-1e-6 and r in set(mm.get_exchange_reactions())
        and any(k.replace('M_','').rsplit('_',1)[0] in SUGAR for k in mm.reactions[r].stoichiometry)]
    print("     3. sugar uptake on ETO_METH : %s  %s" % (up or 'none', "OK" if not up else "FAIL"))
    if up: fails+=1
    # and the check 43b never had: does it eat sugar when sugar is offered?
    t=m.copy()
    Environment.from_compounds(sed, fmt_func=lambda x:"R_EX_{}_e".format(x)).apply(t,inplace=True,warning=False)
    s2=FBA(t); g2=s2.fobj if s2 and s2.fobj else 0.0
    up2=[]
    if s2 and s2.values:
        up2=[(r,round(v,2)) for r,v in s2.values.items()
             if v<-1e-6 and r in set(t.get_exchange_reactions())
             and any(k.replace('M_','').rsplit('_',1)[0] in SUGAR for k in t.reactions[r].stoichiometry)]
    print("     (on ETO_SED: growth %.4f, sugar uptake %s)" % (g2, up2 or 'none'))
print()
print("GATE_STATUS: %s  (%d assertion failure(s))" % ("PASSED" if fails==0 else "FAILED", fails))
sys.exit(0 if fails==0 else 3)
PYEOF
GRC=$?
if [ "$GRC" -ne 0 ]; then
  bad "the rebuilt archaea did not pass. Do NOT run 43l on these."
  note "a model that fails assertion 2 is not a methanogen; one that fails 3 is a"
  note "sugar fermenter wearing an archaeal name. Both were true of the 43b model."
  echo; echo "STEP_STATUS: FAILED  $(date -Is)"; exit 1
fi
ok "both archaea grow on ETO_METH, make methane, and eat no sugar"
say "DONE"
echo "models in $OUT/models -- 43l assembles the 11-member community from these"
echo "STEP_STATUS: COMPLETE  $(date -Is)"
