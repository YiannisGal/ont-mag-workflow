#!/usr/bin/env bash
# 41b_fixes.sh -- repairs to 41, plus the MTP4 question 41 raised.
#
# 41's section A read cobrar's bundled e_coli_core.xml example out of the DAMAGED
# pixi env. The locator had only negative tests. A now requires POSITIVE evidence of
# gapseq provenance and prints every candidate it rejected and why.
#
# A. gapseq model locator, strict          (41's namespace probe, re-armed)
# B. the 8th archaeon: MTP4 -- what is it, and is it a methylotroph?
# C. 07_mags layout discovery, then ANI    (41 D found <2 FASTAs)
# D. aromatic screen, de-duplicated and ID-mapped (41 E counted 59 tables for 29 MAGs
#    and compared ETO_* ids against das_* names, so its in-situ line was meaningless)
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"; export TMPDIR="${TMPDIR41:-${SCRATCH_ROOT}/tmp41}"
mkdir -p "$OUT" "$TMPDIR" || { echo "FATAL"; exit 1; }
CULT="${CULT:-$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis}"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [FAIL] %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
echo "41b_fixes  $(date -Is)"

say "A. gapseq model locator -- positive provenance required"
note "rejecting: .pixi envs, R libraries, extdata, anything marked damaged"
mapfile -t ALLXML < <(find "$CULT" -name '*.xml' 2>/dev/null | sort)
note "${#ALLXML[@]} .xml files under \$CULT in total"
GOOD=(); 
for x in "${ALLXML[@]}"; do
  base="${x%.xml}"; why=""
  case "$x" in
    *.pixi*)      why="inside a pixi env" ;;
    */library/*)  why="inside an R library" ;;
    *extdata*)    why="package example data" ;;
    *damaged*)    why="inside a directory marked damaged" ;;
    *-draft.xml)  why="is a draft, not a filled model" ;;
  esac
  # POSITIVE test: gapseq writes <name>-draft.xml beside its filled <name>.xml
  if [ -z "$why" ] && [ ! -f "${base}-draft.xml" ]; then
    why="no sibling -draft.xml -- not a gapseq fill product"
  fi
  if [ -n "$why" ]; then printf '           reject  %-58s  %s\n' "$(basename "$x")" "$why"
  else GOOD+=("$x"); printf '           ACCEPT  %s\n' "$x"; fi
done
if [ "${#GOOD[@]}" -eq 0 ]; then
  no "no file passes the gapseq provenance test"
  no "the namespace is UNDETERMINED -- do not write the environment file"
elif [ "${#GOOD[@]}" -gt 2 ]; then
  no "${#GOOD[@]} candidates; expected 2 (YGCA01_Clean, YGCA01_Contaminant). Not choosing."
else
  ok "${#GOOD[@]} genuine gapseq model(s)"
  python3 - "${GOOD[0]}" "$OUT" <<'PY'
import re,sys,collections
p,out=sys.argv[1],sys.argv[2]; t=open(p,encoding='utf-8',errors='replace').read()
rx=sorted({m[1] for m in re.findall(r'id="(R_)?(EX_[A-Za-z0-9_\-\[\]]+)"',t)})
print(f"  ...    {len(rx)} exchange reactions"); [print(f"           {r}") for r in rx[:6]]
if len(rx)<50:
    print(f"  [FAIL] only {len(rx)} exchanges -- a genome-scale anaerobe model has hundreds.")
    print("         This looks like a toy/core model. Refusing."); sys.exit(0)
c=collections.Counter()
for r in rx:
    b=r[3:]
    c['ModelSEED cpd#####' if re.match(r'^cpd\d{5}',b) else
      'KEGG C#####' if re.match(r'^C\d{5}',b) else 'BiGG-like/other']+=1
print("  ...    namespace tally:"); [print(f"           {v:>6}  {k}") for k,v in c.most_common()]
top,n=c.most_common(1)[0]
print(f"  [{'OK' if n/len(rx)>=.9 else 'FAIL'}]   namespace {top} at {n/len(rx):.0%}")
nm=dict(re.findall(r'<species[^>]*\bid="([^"]+)"[^>]*\bname="([^"]*)"',t))
low={k:(v or '').strip().lower() for k,v in nm.items()}
open(f"{out}/gapseq_metabolite_names.tsv","w").write(
    "metabolite_id\tname\n"+"".join(f"{k}\t{v}\n" for k,v in sorted(nm.items())))
print(f"  [OK]   wrote {out}/gapseq_metabolite_names.tsv  ({len(nm)} species)")
WANT={"ammonium":["ammonium","nh4"],"phosphate":["phosphate"],"sulfate":["sulfate"],
 "sulfide":["h2s","hydrogen sulfide","sulfide"],"bicarbonate":["bicarbonate","hco3","co2"],
 "water":["h2o"],"proton":["h+"],"potassium":["k+"],"sodium":["na+"],"magnesium":["mg2+"],
 "calcium":["ca2+"],"chloride":["cl-"],"iron2":["fe2+"],"cellobiose":["cellobiose"],
 "xylose":["xylose","d-xylose"],"galacturonate":["d-galacturonate","galacturonate"],
 "acetate":["acetate"],"h2":["h2","hydrogen"],"X-methanol":["methanol"],
 "X-TMA":["trimethylamine"],"X-choline":["choline"],"X-B12":["cobalamin"]}
miss=[];amb=[]
for lab,keys in WANT.items():
    ex=sorted({k for k,v in low.items() if v in keys})
    sb=sorted({k for k,v in low.items() if any(q in v for q in keys)}) if not ex else []
    hit,how=(ex,"EXACT") if ex else (sb,"SUBSTR"); req=not lab.startswith("X-")
    if not hit:
        if req: miss.append(lab)
        print(f"           [--   ] {lab:<16}   0"); continue
    if how=="SUBSTR" and len(hit)>1 and req: amb.append(f"{lab} -> {', '.join(hit[:4])}")
    print(f"           [{'OK   ' if how=='EXACT' else 'AMBIG'}] {lab:<16} {len(hit):>3} {how}  {hit[0]}  {low[hit[0]]!r}")
if miss: print(f"  [FAIL] missing: {', '.join(miss)}")
if amb:  print("  [FAIL] ambiguous:"); [print(f"           {a}") for a in amb]
print("  [FAIL] DO NOT write the environment file." if (miss or amb)
      else "  [OK]   every compound resolves uniquely -- environment file can be written")
PY
fi

say "B. the 8th archaeon -- MTP4"
TAX="$WS2/07_mags/MAG_TAXONOMY.tsv"
if [ ! -f "$TAX" ]; then no "MAG_TAXONOMY.tsv missing"; else
python3 - "$TAX" "$WS2" "$OUT" <<'PY'
import sys,csv,os,glob,re,collections
tax,ws2,out=sys.argv[1],sys.argv[2],sys.argv[3]
rows=list(csv.DictReader(open(tax),delimiter='\t'))
nk=next(k for k in rows[0] if k.strip().lower()=='working_name')
dk=next((k for k in rows[0] if k.strip().lower() in ('domain','kingdom')),None)
arch=[r for r in rows if dk and 'archae' in (r[dk] or '').lower()]
print(f"  ...    {len(arch)} archaea; full lineage of each:")
lin=[k for k in rows[0] if k.strip().lower() in
     ('domain','phylum','class','order','family','genus','species')]
for r in arch:
    print(f"           {r[nk]}")
    print("             "+" | ".join(f"{r[k]}" for k in lin if r.get(k)))
# methanogenesis / methylotrophy markers, by product string across each archaeon
MK={'mcrA  methyl-CoM reductase':[r'methyl-coenzyme m reductase'],
    'mtaB  methanol:CoM':[r'methanol.*methyltransferase',r'\bmtaB\b'],
    'mttB  TMA:CoM':[r'trimethylamine.*methyltransferase',r'\bmttB\b'],
    'mtmB  MMA:CoM':[r'monomethylamine.*methyltransferase',r'\bmtmB\b'],
    'mtbB  DMA:CoM':[r'dimethylamine.*methyltransferase',r'\bmtbB\b'],
    'fwd/fmd  CO2 red.':[r'formylmethanofuran dehydrogenase'],
    'hdrABC':[r'heterodisulfide reductase']}
print("  ...    methanogenesis markers by product string (SCREEN, not curated HMM):")
hdr="           "+f"{'genome':<24}"+"".join(f"{k.split()[0]:>10}" for k in MK)
print(hdr)
for r in arch:
    g=r[nk]; d=os.path.join(ws2,'09_annotation',g)
    f=os.path.join(d,g+'.tsv')
    if not os.path.exists(f):
        cand=[x for x in glob.glob(os.path.join(d,'*.tsv')) if 'hypothetical' not in x]
        f=cand[0] if cand else None
    if not f: print(f"           {g:<24}  (no annotation table)"); continue
    txt=open(f,encoding='utf-8',errors='replace').read().lower()
    line=f"           {g:<24}"
    for k,pats in MK.items():
        n=sum(len(re.findall(p,txt)) for p in pats)
        line+=f"{n if n else '.':>10}"
    print(line)
print("  ...    mcrA present => a methanogen. mtaB/mttB/mtmB/mtbB present => methylotrophic.")
print("  ...    If MTP4 carries mcrA + any mt*B, the in-situ methylotrophic CONSUMER")
print("         problem in Part 1 changes: MTP4 sits at sed_breadth 0.6170, not 0.08-0.18.")
PY
fi

say "C. 07_mags layout, then ANI"
note "what is actually in 07_mags:"
find "$WS2/07_mags" -maxdepth 3 -type d 2>/dev/null | sed "s#$WS2/#           #" | head -12
note "genome-file extensions present, with counts and depth:"
find "$WS2/07_mags" -maxdepth 4 -type f \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \) 2>/dev/null \
 | awk -F/ '{n=split($0,a,"."); ext=a[n]; d=NF; print ext"  depth="d}' | sort | uniq -c | sed 's/^/           /'
NEXP=$(( $(wc -l < "$WS2/07_mags/MAG_TAXONOMY.tsv") - 1 ))
note "MAG_TAXONOMY lists $NEXP genomes -- the genome directory must hold that many"
BD=""; BEST=0
while IFS= read -r cand; do
  c=$(find "$cand" -maxdepth 1 -type f \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \) \
      ! -name '*all_mags*' 2>/dev/null | wc -l)
  [ "$c" -gt 0 ] && printf '           %-40s %s file(s)\n' "${cand#$WS2/}" "$c"
  [ "$c" -gt "$BEST" ] && { BEST=$c; BD="$cand"; }
done < <(find "$WS2/07_mags" -maxdepth 3 -type d 2>/dev/null | sort)
if [ "$BEST" -eq 0 ]; then BD=""
elif [ "$BEST" -ne "$NEXP" ]; then
  no "richest directory ${BD#$WS2/} holds $BEST genome files but taxonomy lists $NEXP"
  no "refusing -- a partial genome set gives a partial ANI matrix that looks complete"
  BD=""
else
  ok "using ${BD#$WS2/}: $BEST files, exactly matching the taxonomy count"
fi
if [ -z "$BD" ]; then no "genome directory not resolved -- ANI not attempted"
elif ! command -v skani >/dev/null 2>&1; then no "skani not on PATH"
else
  find "$BD" -maxdepth 1 -type f \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \) \
    ! -name '*all_mags*' | sort > "$TMPDIR/g.txt"
  ok "$(wc -l < "$TMPDIR/g.txt") genomes -> skani"
  timeout 1800 skani dist --ql "$TMPDIR/g.txt" --rl "$TMPDIR/g.txt" -t 4 \
    -o "$OUT/mag_x_mag_ani.tsv" 2>"$OUT/skani.err"; rc=$?
  if [ $rc -ne 0 ]; then no "skani exit $rc"; else
    ok "$(( $(wc -l < "$OUT/mag_x_mag_ani.tsv") - 1 )) pairs -> $OUT/mag_x_mag_ani.tsv"
    echo "  ...    SED x LIQ pairs at >=95% ANI (the Chakrabartyella pair ~99.75% is the control):"
    awk -F'\t' 'NR>1&&$3>=95{r=$1;q=$2;gsub(/.*\//,"",r);gsub(/.*\//,"",q);
      sub(/\.[^.]*$/,"",r);sub(/\.[^.]*$/,"",q); if(r!=q)printf "           %6.2f%%  %-26s %s\n",$3,r,q}' \
      "$OUT/mag_x_mag_ani.tsv" | sort -rn -k1 | head -20
  fi
fi

say "D. aromatic screen, repaired"
note "41 screened 59 tables for 29 genomes (Bakta writes a .hypotheticals table too)"
python3 - "$WS2" "$OUT" <<'PY'
import sys,os,glob,csv,re,collections
ws2,out=sys.argv[1],sys.argv[2]
tax=os.path.join(ws2,'07_mags','MAG_TAXONOMY.tsv')
rows=list(csv.DictReader(open(tax),delimiter='\t'))
nk=next(k for k in rows[0] if k.strip().lower()=='working_name')
names=[r[nk] for r in rows]
print(f"  ...    {len(names)} genomes from MAG_TAXONOMY")
used=[];skipped=0
for g in names:
    f=os.path.join(ws2,'09_annotation',g,g+'.tsv')
    if os.path.exists(f): used.append((g,f))
    else: skipped+=1
print(f"  ...    exact <dir>/<dir>.tsv matched for {len(used)}; {skipped} missing")
allt=glob.glob(os.path.join(ws2,'09_annotation','*','*.tsv'))
print(f"  ...    a bare glob would have matched {len(allt)} -- that was 41's error")
PAT=re.compile(r'benzoyl-coa|benzoate|4-hydroxybenzo|phenylphosphate|catechol|resorcinol|'
               r'cyclohexa-1,5-diene|phenol\b|bamb|bcra',re.I)
per={}
for g,f in used:
    seen=set()
    with open(f,encoding='utf-8',errors='replace') as fh:
        for line in fh:
            p=line.rstrip('\n').split('\t')
            if len(p)<2: continue
            if PAT.search(line): seen.add(p[0])   # de-duplicate BY LOCUS
    if seen: per[g]=len(seen)
tot=sum(per.values())
print(f"  ...    {tot} distinct loci across {len(per)} genomes (41 said 68 across 19, not deduped)")
insitu=os.path.join(out,'in_situ_10.txt')
ids=[l.strip() for l in open(insitu)] if os.path.exists(insitu) else []
# 41 compared ETO_* ids against das_* names. Map through the crosswalk.
xw=None
for c in glob.glob(os.path.join(ws2,'**','MAG_CROSSWALK.tsv'),recursive=True): xw=c;break
m={}
if xw:
    xr=list(csv.DictReader(open(xw),delimiter='\t'))
    cols=list(xr[0].keys())
    ec=next((k for k in cols if any(str(r[k]).startswith('ETO_') for r in xr)),None)
    wc=next((k for k in cols if any(str(r[k]) in names for r in xr)),None)
    print(f"  ...    crosswalk {os.path.basename(xw)}: ETO column={ec!r} working column={wc!r}")
    if ec and wc: m={r[ec]:r[wc] for r in xr}
if not m:
    print("  [FAIL] no ETO->working_name mapping; cannot restrict to the in-situ set")
else:
    tgt={m.get(i,i) for i in ids}
    print(f"  ...    in-situ 10 as working names: {', '.join(sorted(tgt))}")
    hits={g:n for g,n in per.items() if g in tgt}
    if hits:
        print(f"  [OK]   {len(hits)} of the in-situ 10 carry aromatic-degradation hits:")
        for g,n in sorted(hits.items(),key=lambda x:-x[1]): print(f"           {n:>3}  {g}")
    else:
        print("  [OK]   none of the in-situ 10 carry a hit (now a real negative, not an ID mismatch)")
print("  ...    top genomes overall:")
for g,n in sorted(per.items(),key=lambda x:-x[1])[:12]: print(f"           {n:>3}  {g}")
open(os.path.join(out,'aromatic_screen_dedup.tsv'),'w').write(
    "genome\tdistinct_loci\n"+"".join(f"{g}\t{n}\n" for g,n in sorted(per.items(),key=lambda x:-x[1])))
print(f"  [OK]   wrote {out}/aromatic_screen_dedup.tsv")
PY

say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
