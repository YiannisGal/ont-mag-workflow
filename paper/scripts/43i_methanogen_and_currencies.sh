#!/usr/bin/env bash
# 43i_methanogen_and_currencies.sh -- can the community feed the methanogen?
#
# RUN IN tmux. Section G is optional and slow; set SKIP_RERUN=1 to stop before it.
#
# WHAT 43h ESTABLISHED, and one thing it got wrong
#   models_sed  : filled on ETO_SED, 10/10 grow alone, 0 interactions. Prototrophs.
#   models      : filled on ETO_FILL, 0/10 grow alone, many interactions. Auxotrophs.
#   Neither is the truth; both restate the gapfilling instruction. That dependence is a
#   methods statement, not something to resolve by preferring the prettier table.
#
#   MY SECTION B VERDICT WAS WRONG. It read R_Community_Growth = 0 and announced the
#   community could not grow. R_Community_Growth is an all-members objective; section C
#   then showed EIGHT of ten reach 1.1657 inside that same community on that same medium
#   while all ten are dead alone. That is cross-feeding, measured. The aggregate is zero
#   only because two members are dead and one zero kills the product. I tested the
#   aggregate and reported it as a property of the community. Rule 70, again.
#
#   The two that stay dead: LIQ_das_vamb_017 (Methanofollis) and LIQ_dassub_semibin_017.
#   That is also the honest explanation of mip=n/a: all members at min_growth is
#   infeasible. And mro=n/a for this set is correct, since MRO needs each organism to
#   grow alone. Neither is a bug.
#
# THE QUESTION HERE
#   ETO_SED has no H2, no formate, no acetate, no methanol, no methylamines. Every
#   methanogenic substrate is absent by construction, because those are FERMENTATION
#   PRODUCTS -- the community is supposed to make them. So either the fermenters are not
#   exporting them, or the archaeon cannot take them up. With Bakta missing 95% of
#   archaeal markers and this being the smallest model at 880 reactions, an
#   under-annotated archaeon is a live possibility and must be distinguished from a
#   real biological absence.
#
# IDs are NOT taken from memory. Section A finds compounds BY NAME in the models
# themselves and reports the IDs it found, the way 41c had to.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/43_smetana/run_43i"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; CLI="$WS2/envs/smetana/bin/smetana"
MDB="$WS2/41_env_checks/etoliko_media.tsv"
GAP="$WS2/43_smetana/models"; SED="$WS2/43_smetana/models_sed"
ARC=LIQ_das_vamb_017
SKIP_RERUN="${SKIP_RERUN:-0}"
export PYTHONWARNINGS="ignore::UserWarning"
echo "43i  $(date -Is)"; echo "out: $OUT"

# ---------------------------------------------------------------- 0. salvage
say "0. the detailed table from 43h, read from disk"
note "the live output was scrambled by the SSH drops; this re-reads the file."
DET="$WS2/43_smetana/run_43h/gapped/gap_detailed.tsv"
if [ -s "$DET" ]; then
  N=$(($(wc -l < "$DET") - 1)); echo "   rows: $N"
  echo "   --- donors (n compounds donated, smetana>0) ---"
  tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$4]++} END{for(x in c) printf "     %-26s %d\n",x,c[x]}' | sort -k2,2nr
  echo "   --- receivers ---"
  tail -n +2 "$DET" | awk -F'\t' '$9>0{c[$3]++} END{for(x in c) printf "     %-26s %d\n",x,c[x]}' | sort -k2,2nr
  echo "   --- is the archaeon anywhere in this table? ---"
  A=$(tail -n +2 "$DET" | awk -F'\t' -v a="$ARC" '$3==a && $9>0' | wc -l)
  B=$(tail -n +2 "$DET" | awk -F'\t' -v a="$ARC" '$4==a && $9>0' | wc -l)
  echo "     $ARC as RECEIVER: $A row(s);  as DONOR: $B row(s)"
  [ "$A" -gt 0 ] && { echo "     compounds it receives:"; tail -n +2 "$DET" | awk -F'\t' -v a="$ARC" '$3==a && $9>0{printf "       %-16s from %-24s smetana=%s\n",$5,$4,$9}' | head -20; }
else
  no "no gap_detailed.tsv at $DET"
fi

# ---------------------------------------------------------------- A. currencies
say "A. the currencies, found BY NAME in the models (not by remembered ID)"
"$PY" - "$GAP" "$OUT" <<'PYEOF'
import sys, os, glob, re, json, warnings
from reframed import load_cbmodel
gap, out = sys.argv[1], sys.argv[2]
WANT = {
 'H2'            : r'^h2-[ce]0?$|^hydrogen',
 'Formate'       : r'^formate',
 'Acetate'       : r'^acetate',
 'Methanol'      : r'^methanol',
 'Methylamine'   : r'^methylamine|^monomethylamine',
 'Dimethylamine' : r'^dimethylamine',
 'Trimethylamine': r'^trimethylamine',
 'CO2'           : r'^co2-',
 'H2S'           : r'^h2s-|^hydrogen sulfide',
 'Sulfate'       : r'^sulfate-',
 'Pyruvate'      : r'^pyruvate',
 'L-Lactate'     : r'^l-lactate',
 'Ethanol'       : r'^ethanol',
 'Butyrate'      : r'^butyrate',
 'Propionate'    : r'^propionate',
}
models = {}
for p in sorted(glob.glob(os.path.join(gap,'*.xml'))):
    n = os.path.basename(p).replace('.xml','')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        models[n] = load_cbmodel(p, flavor='fbc2')

# resolve each name to an ID using every model, so a gap in one does not lose it
resolved = {}
for label, pat in WANT.items():
    rx = re.compile(pat, re.I)
    hits = {}
    for n, m in models.items():
        for mid, met in m.metabolites.items():
            if rx.match(met.name or ''):
                cpd = mid.replace('M_','').rsplit('_',1)[0]
                hits.setdefault(cpd, met.name)
    resolved[label] = hits
    print("   %-15s -> %s" % (label, ', '.join('%s (%s)' % (k, v) for k, v in list(hits.items())[:3]) or 'NOT FOUND IN ANY MODEL'))
json.dump({k: list(v) for k, v in resolved.items()}, open(os.path.join(out,'currency_ids.json'),'w'), indent=2)

print()
print("   E = has an exchange reaction;  t = has a transporter;  . = metabolite absent")
labels = [l for l in WANT if resolved[l]]
w = max(len(n) for n in models)
print("   %-*s  %s" % (w, 'model', '  '.join('%-6s' % l[:6] for l in labels)))
table = {}
for n, m in sorted(models.items()):
    exch = set(m.get_exchange_reactions()); cells = []; row = {}
    for l in labels:
        best = '.'
        for cpd in resolved[l]:
            cm, em = 'M_%s_c' % cpd, 'M_%s_e' % cpd
            has_ex = ('R_EX_%s_e' % cpd) in exch
            has_tr = any(cm in r.stoichiometry and em in r.stoichiometry for r in m.reactions.values())
            if has_ex and has_tr: best = 'Et'; break
            if has_ex and best not in ('Et',): best = 'E'
            elif has_tr and best == '.': best = 't'
            elif (cm in m.metabolites or em in m.metabolites) and best == '.': best = '-'
        cells.append(best); row[l] = best
    table[n] = row
    print("   %-*s  %s" % (w, n, '  '.join('%-6s' % c for c in cells)))
json.dump(table, open(os.path.join(out,'currency_table.json'),'w'), indent=2)
PYEOF

# ---------------------------------------------------------------- B. why dead
say "B. why can the methanogen not grow, even in the community?"
note "for each biomass precursor of $ARC, ask whether the model can make it at all"
note "when everything importable is freely available. A precursor that is unreachable"
note "under unlimited uptake is a hole in the RECONSTRUCTION, not a shortage in the medium."
"$PY" - "$GAP" "$ARC" "$OUT" <<'PYEOF'
import sys, os, warnings, json
from reframed import load_cbmodel, Environment, FBA
gap, arc, out = sys.argv[1], sys.argv[2], sys.argv[3]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    m = load_cbmodel(os.path.join(gap, arc + '.xml'), flavor='fbc2')
exch = set(m.get_exchange_reactions())
bio = [r for r in m.reactions if r.lower().startswith('r_bio') and r not in exch]
print("   %s: %d reactions, %d exchanges, biomass=%s" % (arc, len(m.reactions), len(exch), bio))
if not bio: sys.exit(0)
BIO = bio[0]

mm = m.copy()
Environment.complete(mm, max_uptake=1000.0).apply(mm, inplace=True, warning=False)
s = FBA(mm); g = s.fobj if s and s.fobj else 0.0
print("   growth with EVERYTHING it has an exchange for, unlimited: %.6f" % g)
if g < 1e-6:
    print("   >>> it cannot grow even given every compound it can import.")
    print("   >>> the reconstruction is incomplete. No medium and no partner can fix that.")
else:
    print("   >>> it CAN grow when fed. So the block is what is available, not the network.")

subs = [k for k, v in m.reactions[BIO].stoichiometry.items() if v < 0]
print("   biomass precursors: %d" % len(subs))
blocked = []
for cm in subs:
    t = m.copy()
    Environment.complete(t, max_uptake=1000.0).apply(t, inplace=True, warning=False)
    r_id = 'SINK_TEST'
    try:
        t.add_reaction_from_str('%s: %s --> ' % (r_id, cm))
    except Exception:
        continue
    t.set_objective({r_id: 1})
    s = FBA(t); v = s.fobj if s and s.fobj else 0.0
    if v < 1e-9:
        blocked.append((cm, m.metabolites[cm].name))
print("   precursors that CANNOT be produced even with unlimited uptake: %d" % len(blocked))
for cm, nm in blocked[:25]:
    print("     %-20s %s" % (cm, nm))
json.dump([[a, b] for a, b in blocked], open(os.path.join(out, 'archaeon_blocked.json'), 'w'), indent=2)
if blocked:
    print("   >>> these are gaps in the archaeal reconstruction. Given Bakta called 95%")
    print("       of this genome's KOfam-confirmed markers 'hypothetical protein', that")
    print("       is the expected failure mode and it is measurable, not speculative.")
PYEOF

# ---------------------------------------------------------------- C. rescue
say "C. does a methanogenic substrate rescue it? (ETO_SED + one addition at a time)"
note "H2, formate, acetate, methanol, methylamines are FERMENTATION PRODUCTS. Adding"
note "them to the environment asks a different question from adding them to the lagoon:"
note "it asks whether the model can use them at all."
"$PY" - "$GAP" "$SED" "$ARC" "$MDB" "$OUT" <<'PYEOF'
import sys, os, csv, json, warnings
from reframed import load_cbmodel, Environment, FBA
gap, sed, arc, mdb, out = sys.argv[1:6]
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['medium']=='ETO_SED': med.append(row['compound'].strip())
ids = json.load(open(os.path.join(out,'currency_ids.json')))
def load(d, n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_cbmodel(os.path.join(d, n+'.xml'), flavor='fbc2')
for setname, d in (('ETO_FILL-filled', gap), ('ETO_SED-filled', sed)):
    m = load(d, arc)
    exch = set(m.get_exchange_reactions())
    print("\n   [%s] %s" % (setname, arc))
    base = m.copy()
    Environment.from_compounds(med, fmt_func=lambda x: "R_EX_{}_e".format(x)).apply(base, inplace=True, warning=False)
    s = FBA(base); print("     ETO_SED alone                         %.4f" % (s.fobj if s and s.fobj else 0.0))
    for label, cpds in ids.items():
        for cpd in cpds:
            if ('R_EX_%s_e' % cpd) not in exch: continue
            t = m.copy()
            Environment.from_compounds(med+[cpd], fmt_func=lambda x: "R_EX_{}_e".format(x)).apply(t, inplace=True, warning=False)
            s = FBA(t); g = s.fobj if s and s.fobj else 0.0
            flag = "   <-- RESCUED" if g > 1e-6 else ""
            print("     + %-14s (%s)  %.4f%s" % (label, cpd, g, flag))
            break
PYEOF

# ---------------------------------------------------------------- D. exporters
say "D. does anyone EXPORT what the methanogen needs?"
note "an uptake route is useless if no member secretes the compound. This asks each"
note "of the other nine to maximise export of each currency while growing."
"$PY" - "$GAP" "$ARC" "$MDB" "$OUT" <<'PYEOF'
import sys, os, csv, glob, json, warnings
from reframed import load_cbmodel, Environment, FBA
gap, arc, mdb, out = sys.argv[1:5]
med=[]
with open(mdb) as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['medium']=='ETO_SED': med.append(row['compound'].strip())
ids = json.load(open(os.path.join(out,'currency_ids.json')))
want = {}
for label, cpds in ids.items():
    if label in ('H2','Formate','Acetate','Methanol','Methylamine','Dimethylamine','Trimethylamine'):
        for c in cpds: want.setdefault(label, []).append(c)
print("   maximum secretion while growing, on ETO_SED (0 = cannot or will not secrete)")
rows=[]
for p in sorted(glob.glob(os.path.join(gap,'*.xml'))):
    n = os.path.basename(p).replace('.xml','')
    if n == arc: continue
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = load_cbmodel(p, flavor='fbc2')
    exch = set(m.get_exchange_reactions())
    out_row = {}
    for label, cpds in want.items():
        best = 0.0
        for cpd in cpds:
            r = 'R_EX_%s_e' % cpd
            if r not in exch: continue
            t = m.copy()
            Environment.from_compounds(med, fmt_func=lambda x: "R_EX_{}_e".format(x)).apply(t, inplace=True, warning=False)
            t.set_objective({r: 1})
            s = FBA(t); v = s.fobj if s and s.fobj else 0.0
            best = max(best, v)
        out_row[label] = best
    rows.append((n, out_row))
labels = sorted(want)
w = max(len(n) for n,_ in rows)
print("   %-*s  %s" % (w,'model','  '.join('%-9s'%l[:9] for l in labels)))
for n, r in rows:
    print("   %-*s  %s" % (w, n, '  '.join('%-9.3f'%r.get(l,0.0) for l in labels)))
PYEOF

if [ "$SKIP_RERUN" = "1" ]; then note "SKIP_RERUN=1"; say DONE; echo "STEP_STATUS: COMPLETE  $(date -Is)"; exit 0; fi

# ---------------------------------------------------------------- E. medium
say "E. ETO_SED_FERM -- ETO_SED plus the fermentation products, as a stated variant"
note "this is NOT a correction to ETO_SED. It is a second, declared environment that"
note "asks what the community does when fermentation products are present. Both go in"
note "the methods; neither replaces the other."
"$PY" - "$MDB" "$OUT" <<'PYEOF'
import sys, os, csv, json
mdb, out = sys.argv[1], sys.argv[2]
ids = json.load(open(os.path.join(out,'currency_ids.json')))
add=[]
for label in ('H2','Formate','Acetate','Methanol','Methylamine','Dimethylamine','Trimethylamine'):
    for c in ids.get(label,[]): add.append((c,label)); break
rows=[]; hdr=None
with open(mdb) as fh:
    rd = csv.reader(fh, delimiter='\t'); hdr = next(rd)
    for r in rd: rows.append(r)
have = {r[2] for r in rows if r[0]=='ETO_SED'}
new = [r[:] for r in rows if r[0]=='ETO_SED']
for r in new: r[0]='ETO_SED_FERM'
n_add=0
for cpd,label in add:
    if cpd in have: continue
    new.append(['ETO_SED_FERM','ETO_SED plus fermentation products',cpd,label]); n_add+=1
if any(r[0]=='ETO_SED_FERM' for r in rows):
    print("   ETO_SED_FERM already in the media db -- not rewriting")
else:
    with open(mdb,'a') as fh:
        w=csv.writer(fh, delimiter='\t')
        for r in new: w.writerow(r)
    print("   appended ETO_SED_FERM: %d compounds (%d added to ETO_SED's 24)" % (len(new), n_add))
    for c,l in add: print("     %-12s %s" % (c,l))
PYEOF
say "F. SMETANA on ETO_SED_FERM, gap-retaining set, global only (fast)"
TD="$OUT/ferm"; mkdir -p "$TD"; cp "$GAP"/*.xml "$TD/"
( cd "$TD" && "$CLI" -o ferm -m ETO_SED_FERM --mediadb "$MDB" --flavor fbc2 -g "*.xml" 2>&1 ) | sed 's/^/   /' | tail -8
[ -s "$TD/ferm_global.tsv" ] && column -t -s$'\t' "$TD/ferm_global.tsv" | sed 's/^/     /'
note "the detailed run on this medium is the next ~3 h step; not started here."
say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
