#!/usr/bin/env bash
# 41f -- (A) stop hand-building communities and let SMETANA do it, end to end.
#        (B) the mcr answer from 11_kofam, and how badly Bakta missed the archaea.
#
# 41e built Community() itself and got 0 pool exchanges. reframed only creates
# *_pool reactions when extracellular compartments are merged -- flags SMETANA passes
# and I did not. Third time I have tested my reconstruction of a tool instead of the
# tool. So A runs the actual smetana CLI on two real models and looks at the output.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; R="$OUT/relabelled"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41f  $(date -Is)"
PY="$WS2/envs/smetana/bin/python"; SM="$WS2/envs/smetana/bin/smetana"

say "A. reframed Community(): which flags create pool exchanges?"
"$PY" - "$R/YGCA01_Clean.xml" "$R/YGCA01_Contaminant.xml" <<'PY' 2>&1 | grep -v UserWarning | grep -v '^  warn'
import sys, inspect, warnings
warnings.filterwarnings('ignore')
from reframed import load_cbmodel, Community
print("  ...    Community signature:")
for p in inspect.signature(Community.__init__).parameters.values():
    if p.name!='self': print(f"           {p.name} = {p.default!r}")
a,b=sys.argv[1],sys.argv[2]
m1=load_cbmodel(a,flavor='fbc2'); m1.id="M1"
m2=load_cbmodel(b,flavor='fbc2'); m2.id="M2"
for kw in ({}, {"merge_extracellular_compartments":True},
           {"merge_extracellular_compartments":True,"interacting":True},
           {"merge_extracellular_compartments":True,"create_biomass":True}):
    try:
        c=Community("t",[m1,m2],**kw)
        mg=c.merged_model if hasattr(c,'merged_model') else c.merged
        pool=[r for r in mg.reactions if r.endswith('_pool')]
        print(f"  ...    {kw or '(defaults)'} -> {len(pool)} pool exchanges")
        if pool:
            for r in sorted(pool)[:4]: print(f"           {r}")
            ps=set(pool)
            good=[s for s in ("cpd00013","cpd00009","cpd00048","cpd00001")
                  if f"R_EX_M_{s}_e_pool" in ps]
            print(f"           probe compounds addressable: {len(good)}/4")
            if len(good)==4:
                print("  [OK]   relabelled models + these flags = SMETANA can address them")
                break
    except Exception as e:
        print(f"  ...    {kw or '(defaults)'} -> {type(e).__name__}: {str(e)[:70]}")
PY

say "B. end-to-end smoke test: run the real smetana binary"
if [ ! -f "$R/YGCA01_Clean.xml" ]; then no "relabelled models missing -- run 41e first"; else
  printf 'medium\tdescription\tcompound\tname\n' > "$OUT/smoke_media.tsv"
  for c in cpd00001 cpd00009 cpd00013 cpd00048 cpd00030 cpd00034 cpd00063 cpd00099 \
           cpd00205 cpd00254 cpd00971 cpd10515 cpd00067 cpd00011 cpd00027; do
    printf 'SMOKE\tsmoke test\t%s\tx\n' "$c" >> "$OUT/smoke_media.tsv"; done
  ok "smoke medium: 15 compounds"
  printf 'SMOKE\tM1\nSMOKE\tM2\n' > "$OUT/smoke_comm.tsv"
  cd "$R" || exit 1
  cp YGCA01_Clean.xml M1.xml 2>/dev/null; cp YGCA01_Contaminant.xml M2.xml 2>/dev/null
  note "running: smetana --flavor fbc2 --anaerobic -m SMOKE --mediadb ... M1.xml M2.xml"
  timeout 900 "$SM" --flavor fbc2 --anaerobic -m SMOKE --mediadb "$OUT/smoke_media.tsv" \
     -o "$OUT/smoke" M1.xml M2.xml 2>&1 | tail -20 | sed 's/^/           /'
  rc=${PIPESTATUS[0]}
  echo "           exit=$rc"
  F=$(ls "$OUT"/smoke*.tsv 2>/dev/null | head -1)
  if [ -n "$F" ] && [ -s "$F" ]; then
    n=$(( $(wc -l < "$F") - 1 ))
    ok "output $F has $n data row(s)"
    head -3 "$F" | sed 's/^/           /'
    [ "$n" -eq 0 ] && no "ZERO rows -- this is the silent-null failure; addressing is still wrong"
  else
    no "no output file -- smetana did not complete"
  fi
  cd - >/dev/null
fi

say "C. mcr and methylotrophy from 11_kofam -- and the reach of that probe"
python3 - "$WS2" "$OUT" <<'PY'
import sys,os,glob,csv,re,collections
ws2,out=sys.argv[1],sys.argv[2]
K=os.path.join(ws2,'11_kofam')
tgt=os.path.join(K,'kofam_targets.hmm')
KOS={'K00399':'mcrA','K00401':'mcrB','K00402':'mcrG',
     'K04480':'mtaB methanol','K14080':'mtaA','K14081':'mtaC',
     'K14082':'mtmB MMA','K16178':'mtmC','K14083':'mttB TMA','K14084':'mttC',
     'K16177':'mtbB DMA','K00577':'mtrA','K00200':'fwdA'}
# Rule 56: a negative is bounded by the probe. Which KOs were even searched?
searched=set()
if os.path.exists(tgt):
    for line in open(tgt,encoding='utf-8',errors='replace'):
        if line.startswith('NAME'):
            searched.add(line.split()[1])
print(f"  ...    kofam_targets.hmm contains {len(searched)} profiles")
inset=[k for k in KOS if k in searched]; outset=[k for k in KOS if k not in searched]
print(f"  ...    of the {len(KOS)} markers of interest, {len(inset)} were searched:")
print("           searched : "+", ".join(f"{k}({KOS[k]})" for k in inset))
if outset:
    print("           NOT searched: "+", ".join(f"{k}({KOS[k]})" for k in outset))
    print("           -> absence of these is NO evidence. They were never looked for.")
tax=os.path.join(ws2,'07_mags','MAG_TAXONOMY.tsv')
rows=list(csv.DictReader(open(tax),delimiter='\t'))
nk=next(k for k in rows[0] if k.strip().lower()=='working_name')
gk=next(k for k in rows[0] if k.strip().lower()=='genus')
dk=next(k for k in rows[0] if k.strip().lower() in ('domain','kingdom'))
arch=[r for r in rows if 'archae' in (r[dk] or '').lower()]
hyp=collections.Counter(); tot=collections.Counter(); per={}
print()
print(f"           {'genome':<24}{'genus':<20}"+"".join(f"{KOS[k].split()[0]:>9}" for k in inset))
for r in arch:
    g=r[nk]; f=os.path.join(K,'hits',g+'.tbl')
    d=collections.Counter()
    if os.path.exists(f):
        for line in open(f,encoding='utf-8',errors='replace'):
            if line.startswith('#'): continue
            p=line.split()
            if len(p)<6: continue
            ko=p[2]
            if ko in KOS:
                d[ko]+=1
                if 'hypothetical' in line.lower(): hyp[g]+=1
                tot[g]+=1
    per[g]=d
    print(f"           {g:<24}{r[gk]:<20}"+"".join(f"{d[k] if d[k] else '.':>9}" for k in inset))
print()
mcr=[g for g,d in per.items() if d['K00399'] or d['K00401'] or d['K00402']]
print(f"  [{'OK' if len(mcr)>=6 else 'FAIL'}]   mcr in {len(mcr)}/{len(arch)} archaea by KOfam"
      + ("  -- product strings said 0/8" if mcr else ""))
mty={'K04480','K14080','K14081','K14082','K16178','K14083','K14084','K16177'} & set(inset)
mlist=[(g,sum(per[g][k] for k in mty)) for g in per if sum(per[g][k] for k in mty)>0]
print(f"  ...    methylotrophic methyltransferase KOs present in {len(mlist)} archaeon(s):")
for g,n in sorted(mlist,key=lambda x:-x[1]):
    gen=next(r[gk] for r in arch if r[nk]==g)
    print(f"           {g:<24}{gen:<20}{n} hit(s)")
print()
th=sum(hyp.values()); tt=sum(tot.values())
print(f"  🔴 Bakta called {th} of {tt} of these KOfam-confirmed archaeal hits "
      f"'hypothetical protein' ({100*th/tt:.0f}%)" if tt else "  ...    no hits to score")
print("     That is the size of the archaeal annotation gap, measured rather than asserted.")
print("     Any archaeal conclusion in this thesis drawn from Bakta product strings")
print("     needs re-deriving from 11_kofam.")
with open(os.path.join(out,'archaea_kofam_markers.tsv'),'w') as fh:
    fh.write("genome\tgenus\t"+"\t".join(inset)+"\n")
    for r in arch:
        g=r[nk]; fh.write(f"{g}\t{r[gk]}\t"+"\t".join(str(per[g][k]) for k in inset)+"\n")
print(f"  [OK]   wrote {out}/archaea_kofam_markers.tsv")
PY
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
