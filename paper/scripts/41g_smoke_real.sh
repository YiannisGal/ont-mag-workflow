#!/usr/bin/env bash
# 41g -- (A) a smoke test that cannot pass by reading its own inputs.
#        (B) are the methylamine KO counts distinct proteins, or one family hitting
#            three homologous profiles? And do the hits survive ko_list thresholds?
#
# 41f B globbed "$OUT/smoke*.tsv" and matched smoke_comm.tsv -- the communities file it
# had just written itself -- then reported it as smetana output. A check that passes on
# its own setup is not a check. Outputs now go to an EMPTY directory.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; R="$OUT/relabelled"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41g  $(date -Is)"
PY="$WS2/envs/smetana/bin/python"; SM="$WS2/envs/smetana/bin/smetana"

say "A0. versions -- 41f found reframed's Community takes only 3 args"
"$PY" - <<'PY'
import importlib.metadata as md
for p in ("smetana","reframed"):
    try: print(f"  ...    {p} {md.version(p)}")
    except Exception as e: print(f"  ...    {p}: {e}")
PY
note "smetana/interface.py expects R_EX_M_{}_e_pool; reframed 1.6.0 Community() has no"
note "merge flags. If these versions disagree, that is the whole problem."

say "A. smoke test into an EMPTY directory"
S="$OUT/smoke_run"; rm -rf "$S"; mkdir -p "$S"
[ -f "$R/YGCA01_Clean.xml" ] || { no "relabelled models missing"; exit 1; }
cp "$R/YGCA01_Clean.xml" "$S/M1.xml"; cp "$R/YGCA01_Contaminant.xml" "$S/M2.xml"
printf 'medium\tdescription\tcompound\tname\n' > "$S/media.tsv"
for c in cpd00001 cpd00009 cpd00013 cpd00048 cpd00030 cpd00034 cpd00063 cpd00099 \
         cpd00205 cpd00254 cpd00971 cpd10515 cpd00067 cpd00011 cpd00027; do
  printf 'SMOKE\tsmoke\t%s\tx\n' "$c" >> "$S/media.tsv"; done
BEFORE=$(find "$S" -type f | sort)
note "files present before the run:"; echo "$BEFORE" | sed "s#$S/#           #"
cd "$S" || exit 1
note "running smetana (detailed mode, 900 s cap)"
timeout 900 "$SM" --flavor fbc2 --anaerobic -m SMOKE --mediadb media.tsv \
   -o result M1.xml M2.xml > run.log 2>&1
rc=$?
echo "           exit=$rc"
tail -12 run.log | sed 's/^/           /'
AFTER=$(find "$S" -type f -newer media.tsv ! -name run.log ! -name 'M?.xml' ! -name media.tsv | sort)
if [ -z "$AFTER" ]; then
  no "smetana created NO new output file. exit=$rc is not success."
  no "Nothing about a 10-genome run can be trusted until this two-genome run works."
else
  ok "new files created by smetana:"
  echo "$AFTER" | sed "s#$S/#           #"
  for f in $AFTER; do
    n=$(( $(wc -l < "$f") - 1 ))
    printf '           %-28s %s data row(s)\n' "$(basename "$f")" "$n"
    [ "$n" -le 0 ] && no "$(basename "$f") is EMPTY -- the silent-null failure, confirmed"
    [ "$n" -gt 0 ] && head -3 "$f" | sed 's/^/             /'
  done
fi
cd - >/dev/null

say "B. methylamine KOs: distinct proteins or one family, three profiles?"
python3 - "$WS2" "$OUT" <<'PY'
import sys,os,csv,collections
ws2,out=sys.argv[1],sys.argv[2]
K=os.path.join(ws2,'11_kofam','hits')
# per-KO score thresholds. Without them these are raw hmmsearch hits.
kl=None
for c in (os.path.join(ws2,'11_kofam','ko_list'), os.path.join(ws2,'36_gapclose','ko_list')):
    if os.path.exists(c): kl=c; break
TH={}
if kl:
    with open(kl) as fh:
        rd=csv.DictReader(fh,delimiter='\t')
        for r in rd:
            try: TH[r['knum']]=(float(r['threshold']), r.get('score_type',''))
            except Exception: pass
    print(f"  ...    ko_list at {kl}: {len(TH)} thresholds")
else:
    print("  [WARN] no ko_list found -- hits below KEGG's own cutoff cannot be removed.")
    print("         Every count below is therefore an UPPER bound.")
MT=['K04480','K14080','K14082','K14083','K16177']
NAMES={'K04480':'mtaB','K14080':'mtaA','K14082':'mtmB','K14083':'mttB','K16177':'mtbB',
       'K00399':'mcrA','K00401':'mcrB','K00402':'mcrG'}
tax=os.path.join(ws2,'07_mags','MAG_TAXONOMY.tsv')
rows=list(csv.DictReader(open(tax),delimiter='\t'))
nk=next(k for k in rows[0] if k.strip().lower()=='working_name')
gk=next(k for k in rows[0] if k.strip().lower()=='genus')
dk=next(k for k in rows[0] if k.strip().lower() in ('domain','kingdom'))
arch=[r for r in rows if 'archae' in (r[dk] or '').lower()]
print(f"\n           {'genome':<24}{'genus':<20}{'raw hits':>9}{'distinct prot':>15}{'pass thr':>10}")
detail={}
for r in arch:
    g=r[nk]; f=os.path.join(K,g+'.tbl')
    prot=collections.defaultdict(set); raw=0; passed=collections.defaultdict(set)
    if os.path.exists(f):
        for line in open(f,encoding='utf-8',errors='replace'):
            if line.startswith('#'): continue
            p=line.split()
            if len(p)<6: continue
            ko=p[2]
            if ko not in MT: continue
            raw+=1; prot[p[0]].add(ko)
            try: sc=float(p[5])
            except Exception: sc=None
            t=TH.get(ko)
            if t is None or (sc is not None and sc>=t[0]): passed[p[0]].add(ko)
    detail[g]=(prot,passed)
    print(f"           {g:<24}{r[gk]:<20}{raw:>9}{len(prot):>15}{len(passed):>10}")
print("\n  ...    KO co-assignment for the three strong genomes (how many KOs each protein hits):")
for g,(prot,passed) in detail.items():
    if len(prot)<5: continue
    c=collections.Counter(len(v) for v in prot.values())
    gen=next(r[gk] for r in arch if r[nk]==g)
    print(f"           {g}  ({gen})")
    for k in sorted(c): print(f"             {c[k]:>3} protein(s) match {k} KO(s)")
    multi=sum(v for k,v in c.items() if k>1)
    print(f"             -> {multi}/{len(prot)} proteins are assigned to MORE THAN ONE KO."
          + ("  Counts are NOT copy numbers." if multi else ""))
print("\n  ...    per-KO distinct-protein counts (the number that may be quoted):")
hdr="           "+f"{'genome':<24}"+"".join(f"{NAMES[k]:>8}" for k in MT)
print(hdr)
for r in arch:
    g=r[nk]; prot,passed=detail[g]
    line=f"           {g:<24}"
    for ko in MT:
        n=sum(1 for v in passed.values() if ko in v)
        line+=f"{n if n else '.':>8}"
    print(line)
print("\n  ...    'passed thr' uses ko_list where available; where it is not, the raw hit")
print("         stands and the number is an upper bound. State which applies in Methods.")
PY
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
