#!/usr/bin/env bash
# 43j_which_archaeon.sh -- is the modelled archaeon the one the argument rests on?
#                          and is the ETO_SED-filled methanogen a methanogen at all?
#
# 43i settled two things and broke one model.
#
#   SETTLED: the archaeal reconstruction is NOT the problem. LIQ_das_vamb_017 grows at
#   114.0 with unlimited uptake and 0 of its 58 biomass precursors are unreachable.
#   The Bakta annotation gap did not cripple gapseq's model. That hypothesis is dead.
#
#   SETTLED: the substrate supply is real. On ETO_SED the fermenters export H2 (up to
#   169) and formate (up to 111). A hydrogenotroph has partners here.
#   But NOBODY exports methanol or trimethylamine -- every value 0.000 -- and
#   methylamine and dimethylamine are absent from every model.
#
#   BROKEN: the ETO_SED-filled LIQ_das_vamb_017 grows at 0.8429 on a medium with no H2,
#   no formate, no acetate, no methanol and no methylamine. Adding any of them changes
#   nothing. ETO_SED's only carbon is glucose, xylose, cellobiose and galacturonate, so
#   43b's gapfilling gave a Methanofollis SUGAR CATABOLISM. That model is biologically
#   invalid, it was inside "10/10 grow on ETO_SED", and the 0-interaction result was
#   computed with it. Rule 74, sharpened: a gapfilling medium that an organism could not
#   inhabit does not produce a model of that organism.
#
# TWO QUESTIONS HERE, in order of consequence.
#
#   A. WHICH ARCHAEON IS MODELLED? Part 1's coexistence argument rests on MTP4
#      (ETO_LIQ_MAG_003, Methanosarcinaceae, 0.6170 sediment breadth), a METHYLOTROPH.
#      The in-situ ten contain LIQ_das_vamb_017, a Methanofollis -- a HYDROGENOTROPH.
#      If those are different organisms then the community analysis has been modelling
#      an archaeon the argument does not depend on. This must be resolved from the
#      crosswalk, not from the resemblance of two names ending in 003.
#
#   B. IS IT DOING METHANOGENESIS? Section C reports every non-zero exchange flux for
#      the ETO_SED-filled model on ETO_SED: what enters, what leaves, and whether
#      methane is produced at all.
#
# Fast. No 3-hour step. Section F only PROPOSES the fix and prints the command.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43j"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"
GAP="$WS2/43_smetana/models"; SED="$WS2/43_smetana/models_sed"
MDB="$WS2/41_env_checks/etoliko_media.tsv"; ARC=LIQ_das_vamb_017
export PYTHONWARNINGS="ignore::UserWarning"
echo "43j  $(date -Is)"; echo "out: $OUT"

# ---------------------------------------------------------------- A. crosswalk
say "A. WHICH ARCHAEA EXIST, and which one is MTP4?"
note "searching every crosswalk / taxonomy table on disk rather than assuming a path."
mapfile -t XW < <(find "$WS2" -maxdepth 3 \( -name '*CROSSWALK*' -o -name '*crosswalk*' \
  -o -name '*gtdb*.tsv' -o -name '*taxonomy*.tsv' -o -name '*classif*.tsv' \) -type f 2>/dev/null | head -20)
if [ "${#XW[@]}" -eq 0 ]; then no "no crosswalk-like file found under $WS2"; fi
for f in "${XW[@]}"; do
  echo "   --- $f ($(wc -l < "$f") lines) ---"
  head -1 "$f" | sed 's/^/     /'
done
echo
note "every row mentioning Archaea, Methano*, or MTP:"
for f in "${XW[@]}"; do
  if grep -qiE 'archae|methano|MTP' "$f" 2>/dev/null; then
    echo "   --- $f ---"
    grep -iE 'archae|methano|MTP' "$f" | head -20 | sed 's/^/     /'
  fi
done
echo
note "and specifically: where does MTP4 / ETO_LIQ_MAG_003 map, and where does $ARC?"
for f in "${XW[@]}"; do
  for k in MTP4 ETO_LIQ_MAG_003 "$ARC" vamb_017; do
    if grep -q "$k" "$f" 2>/dev/null; then
      echo "   [$k] in $(basename "$f"):"
      grep "$k" "$f" | head -6 | sed 's/^/     /'
    fi
  done
done

# ---------------------------------------------------------------- B. genus
say "B. what does the model itself say it is?"
note "gapseq writes the source genome into the SBML. This reads it, so the identity"
note "does not depend on my reading of a filename."
"$PY" - "$GAP" "$SED" "$ARC" <<'PYEOF'
import sys, os, warnings, re
from reframed import load_cbmodel
gap, sed, arc = sys.argv[1:4]
for tag, d in (('ETO_FILL-filled', gap), ('ETO_SED-filled', sed)):
    p = os.path.join(d, arc + '.xml')
    if not os.path.exists(p): continue
    txt = open(p, encoding='utf-8', errors='replace').read(20000)
    ids = set(re.findall(r'(?:id|name)="([^"]{4,60})"', txt))
    interesting = [i for i in ids if re.search(r'methano|archae|MAG|ETO|bin', i, re.I)]
    print("   [%s] header identifiers: %s" % (tag, sorted(interesting)[:8] or 'none obvious'))
PYEOF

# ---------------------------------------------------------------- C. is it a methanogen
say "C. what is the ETO_SED-filled model actually eating, and does it make methane?"
"$PY" - "$SED" "$GAP" "$ARC" "$MDB" "$OUT" <<'PYEOF'
import sys, os, csv, json, warnings, re
from reframed import load_cbmodel, Environment, FBA
sed, gap, arc, mdb, out = sys.argv[1:6]
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['medium']=='ETO_SED': med.append(row['compound'].strip())

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    m = load_cbmodel(os.path.join(sed, arc + '.xml'), flavor='fbc2')

# find methane by NAME, not by remembered id
meth = {mid: met.name for mid, met in m.metabolites.items()
        if re.match(r'^methane', (met.name or ''), re.I)}
print("   methane metabolites in the model: %s" % (meth or 'NONE'))
mcr = {r: rx.name for r, rx in m.reactions.items()
       if re.search(r'methyl-?coenzyme|coenzyme m reductase|\bmcr\b', (rx.name or ''), re.I)}
print("   methyl-CoM reductase reactions  : %s" % (mcr or 'NONE FOUND BY NAME'))

mm = m.copy()
Environment.from_compounds(med, fmt_func=lambda x: "R_EX_{}_e".format(x)).apply(mm, inplace=True, warning=False)
s = FBA(mm)
g = s.fobj if s and s.fobj else 0.0
print("\n   growth on ETO_SED = %.4f" % g)
if s is None or not s.values:
    print("   no flux distribution returned"); sys.exit(0)
exch = set(mm.get_exchange_reactions())
rows = [(r, v) for r, v in s.values.items() if r in exch and abs(v) > 1e-6]
def nm(r):
    st = mm.reactions[r].stoichiometry
    mid = list(st)[0] if st else None
    return mm.metabolites[mid].name if mid in mm.metabolites else r
print("\n   CONSUMED (uptake, flux < 0):")
for r, v in sorted([x for x in rows if x[1] < 0], key=lambda x: x[1]):
    print("     %-26s %10.3f   %s" % (r, v, nm(r)))
print("\n   PRODUCED (secretion, flux > 0):")
for r, v in sorted([x for x in rows if x[1] > 0], key=lambda x: -x[1]):
    print("     %-26s %10.3f   %s" % (r, v, nm(r)))
ch4 = [(r, v) for r, v in rows if any(k in mm.reactions[r].stoichiometry for k in meth)]
print("\n   methane exchange flux: %s" % (ch4 or 'ZERO -- it is not making methane'))
sugars = [ (r,v) for r,v in rows if v < 0 and re.search(r'glucose|xylose|cellobiose|galacturonate', nm(r) or '', re.I)]
if sugars:
    print("\n   >>> IT IS EATING SUGAR:")
    for r,v in sugars: print("       %-26s %10.3f  %s" % (r,v,nm(r)))
    print("   >>> a Methanofollis does not ferment sugars. This model is an artefact of")
    print("       gapfilling on a medium the organism could not inhabit.")
json.dump({'growth': g, 'consumed': [[r,v] for r,v in rows if v<0],
           'produced': [[r,v] for r,v in rows if v>0]},
          open(os.path.join(out,'archaeon_fluxes.json'),'w'), indent=2)
PYEOF

# ---------------------------------------------------------------- D. donations
say "D. what does the archaeon DONATE in the 43h table? (62 rows, 0 as receiver)"
DET="$WS2/43_smetana/run_43h/gapped/gap_detailed.tsv"
if [ -s "$DET" ]; then
  tail -n +2 "$DET" | awk -F'\t' -v a="$ARC" '$4==a && $9>0{printf "     %-16s to %-24s smetana=%s\n",$5,$3,$9}' \
    | sort -k5,5 | head -25
  echo "   distinct compounds donated: $(tail -n +2 "$DET" | awk -F'\t' -v a="$ARC" '$4==a && $9>0{print $5}' | sort -u | wc -l)"
else no "no $DET"; fi

# ---------------------------------------------------------------- E. methanogenic medium
say "E. ETO_METH -- a medium a methanogen could actually inhabit"
note "H2 + CO2 + N, P, S and the trace metals from ETO_SED. NO SUGARS. This is the"
note "medium the archaeon should have been gapfilled on: gapfilling is only meaningful"
note "on an environment the organism can live in."
"$PY" - "$MDB" "$GAP" "$ARC" "$OUT" <<'PYEOF'
import sys, os, csv, warnings
from reframed import load_cbmodel, Environment, FBA
mdb, gap, arc, out = sys.argv[1:5]
rows=[]; hdr=None
with open(mdb) as fh:
    rd=csv.reader(fh, delimiter='\t'); hdr=next(rd)
    for r in rd: rows.append(r)
SUGARS={'cpd00027','cpd00154','cpd00158','cpd00280'}
base=[r for r in rows if r[0]=='ETO_SED' and r[2] not in SUGARS]
add=[('cpd11640','H2'),('cpd00011','CO2')]
have={r[2] for r in base}
meth=[[ 'ETO_METH','methanogenic: H2/CO2, no sugars', r[2], r[3]] for r in base]
for c,n in add:
    if c not in have: meth.append(['ETO_METH','methanogenic: H2/CO2, no sugars',c,n])
print("   ETO_METH: %d compounds (%d from ETO_SED minus %d sugars, plus H2/CO2)"
      % (len(meth), len(base), len(SUGARS)))
print("   removed: cpd00027 D-Glucose, cpd00154 Xylose, cpd00158 CELB, cpd00280 D-Galacturonate")
if any(r[0]=='ETO_METH' for r in rows):
    print("   already in media db -- not appending")
else:
    with open(mdb,'a') as fh:
        w=csv.writer(fh, delimiter='\t')
        for r in meth: w.writerow(r)
    print("   appended to %s" % mdb)
cpds=[r[2] for r in meth]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    m = load_cbmodel(os.path.join(gap, arc+'.xml'), flavor='fbc2')
mm=m.copy()
Environment.from_compounds(cpds, fmt_func=lambda x: "R_EX_{}_e".format(x)).apply(mm, inplace=True, warning=False)
s=FBA(mm); g=s.fobj if s and s.fobj else 0.0
print("\n   ETO_FILL-filled %s on ETO_METH: %.4f" % (arc, g))
if g < 1e-6:
    print("   >>> as expected: it needs gapfilling on THIS medium, not on ETO_SED.")
    print("   >>> that is the point of section F.")
PYEOF

say "F. the proposed fix -- NOT run here"
note "refill the archaeon on ETO_METH so it stays a methanogen, keep the bacteria as"
note "they are, and re-run the community. The gapseq command is printed rather than"
note "executed: refilling is a rebuild, and the last time I ran one on my own judgement"
note "it produced a sugar-fermenting Methanofollis."
cat <<EOF
     # from the gapseq env used in 42/43b, for LIQ_das_vamb_017 only:
     gapseq fill -m <draft.RDS> -n <ETO_METH medium csv> -c <transporter rxn weights> \\
                 -b <biomass> -g <pathways>
     # then re-verify:  does it grow on ETO_METH, and does methane exchange carry flux?
EOF
note "That is written as 43k once section A tells us WHICH archaeon should be in the"
note "community. Refilling the wrong one is wasted effort."
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
