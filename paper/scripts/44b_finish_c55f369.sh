#!/usr/bin/env bash
# 44b_finish.sh -- the two sections 44 could not run, and the ANI names resolved.
#
# 44 SETTLED, and these are not re-opened:
#   - the in-situ ten are EXACTLY the top ten by sediment breadth, cut at 0.8936.
#     Nothing outside the ten is at or above it. MTP4 (ETO_LIQ_MAG_003) is at 0.6170,
#     below a clean gap (0.8936 -> 0.6912 -> 0.6573 -> 0.6170). The exclusion was a
#     threshold decision, not an oversight. Part 1 stays at ten genomes.
#   - ETO_LIQ_MAG_011, the Methanofollis that receives formate at 0.60, IS in the ten
#     (10th, 0.8936). The syntrophy involves a genuine in-situ community member.
#   - cross_abundance.tsv is clean: 0 values written as '.nnn'. The loss was in the paste.
#   - Methanosarcina: TWO named (ETO_LIQ_MAG_010 0.0825, ETO_LIQ_MAG_018 0.1796, both
#     culture-only). MTP4's GTDB genus is literally 'MTP4', an unnamed Methanosarcinaceae
#     genus -- so "three" counted Methanosarcinaceae and "two" counted Methanosarcina.
#     Both numbers were right about different things.
#
# WHY 44's SECTIONS D AND F FAILED -- both my file-finding, neither a data problem:
#   D  find_one matched gtdbtk.bac120.summary.tsv first. Archaea are in ar53. Searching
#      for one taxonomy file when a tool writes two per domain is the same
#      equal-treatment failure the notebook has flagged since Volume I.
#   F  I globbed "$KOF/*.tbl". 43m found LIQ_das_vamb_003.tbl by RECURSIVE search, so the
#      tables are in a subdirectory. 0 .tbl files parsed, and the section reported that
#      honestly rather than printing zeros.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/44_loose_ends"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; command -v "$PY" >/dev/null 2>&1 || PY=python3
KOF="$WS2/11_kofam"; KOLIST="$KOF/ko_list"
XW="$WS2/07_mags/MAG_CROSSWALK.tsv"; XAB="$WS2/26_crossmap/cross_abundance.tsv"
ANI="$WS2/41_env_checks/mag_x_mag_ani.tsv"; TEN="$WS2/41_env_checks/in_situ_10.txt"
echo "44b  $(date -Is)"

# ---------------------------------------------------------------- A. taxonomy
say "A. the archaeal census -- BOTH gtdbtk summaries this time"
mapfile -t SUMS < <(find "$WS2/08_taxonomy" -name 'gtdbtk.*.summary.tsv' 2>/dev/null | sort)
for f in "${SUMS[@]}"; do printf '  [OK]   %-14s %s rows\n' "$(basename "$f")" "$(($(wc -l < "$f")-1))"; done
[ "${#SUMS[@]}" -eq 0 ] && { no "no gtdbtk summaries under 08_taxonomy"; }
"$PY" - "$XAB" "${SUMS[@]:-}" <<'PYEOF'
import sys, csv, re, os
xab=sys.argv[1]; sums=[a for a in sys.argv[2:] if a and os.path.exists(a)]
tax={}
for f in sums:
    with open(f) as fh:
        rd=csv.DictReader(fh, delimiter='\t')
        c=rd.fieldnames or []
        n=next((x for x in c if re.search(r'user_genome|genome|bin', x, re.I)), c[0])
        t=next((x for x in c if re.search(r'classification|taxonom', x, re.I)), None)
        for r in rd:
            if t: tax[r[n]]=(r[t], os.path.basename(f))
print("   genomes classified across all summaries: %d" % len(tax))
arch={k:v for k,v in tax.items() if 'd__Archaea' in v[0]}
bact={k:v for k,v in tax.items() if 'd__Bacteria' in v[0]}
print("   Archaea %d   Bacteria %d" % (len(arch), len(bact)))
print()
print("   %-26s %-22s %-24s %s" % ('genome','family','genus','species'))
for k,(t,src) in sorted(arch.items()):
    f=re.search(r'f__([^;]*)',t); g=re.search(r'g__([^;]*)',t); s=re.search(r's__([^;]*)',t)
    print("   %-26s %-22s %-24s %s" % (k, f.group(1) if f else '-', g.group(1) if g else '-', s.group(1) if s else '-'))
print()
for key in ('Methanosarcinaceae','Methanosarcina','Methanofollis','Methanomicrobiales','Methanobacteri','Methanogenium'):
    hits=[k for k,(t,_) in tax.items() if key in t]
    print("   %-22s %2d  %s" % (key, len(hits), sorted(hits)))
print()
print("   NOTE: 'Methanosarcina' and 'Methanosarcinaceae' are different counts. MTP4's")
print("   genus is unnamed, so it is Methanosarcinaceae but not Methanosarcina.")
PYEOF

# ---------------------------------------------------------------- B. ANI names
say "B. the ANI pairs, resolved to MAG ids"
"$PY" - "$ANI" "$XW" "$XAB" "$TEN" <<'PYEOF'
import sys, os, csv, re
ani, xw, xab, tenf = sys.argv[1:5]
bin2mag={}; genus={}
if os.path.exists(xw):
    with open(xw) as fh:
        rd=csv.DictReader(fh, delimiter='\t')
        for r in rd:
            bin2mag[r.get('working_name','')]=r.get('mag_id','')
if os.path.exists(xab):
    with open(xab) as fh:
        for r in csv.DictReader(fh, delimiter='\t'): genus[r['mag_id']]=r.get('genus','-')
ten=set()
if os.path.exists(tenf): ten={l.strip() for l in open(tenf) if l.strip()}

def lab(path):
    b=os.path.basename(path).replace('.fa','').replace('.fna','')
    m=bin2mag.get(b,'?')
    return "%s (%s, %s)%s" % (b, m, genus.get(m,'-'), ' [IN THE TEN]' if m in ten else '')

rows=[]
with open(ani) as fh:
    rd=csv.DictReader(fh, delimiter='\t')
    for r in rd:
        try: v=float(r['ANI'])
        except Exception: continue
        rows.append((r['Ref_file'], r['Query_file'], v))
print("   %d comparisons" % len(rows))
hi=[x for x in rows if x[2] >= 95 and os.path.basename(x[0]) != os.path.basename(x[1])]
print("   pairs >= 95%% (same species): %d" % len(hi))
for a,b,v in sorted(hi, key=lambda x:-x[2]):
    print("     %6.2f  %s" % (v, lab(a)))
    print("             %s" % lab(b))

# the specific question: Sedimentibacter SED x LIQ
sed_mags={m for m,g in genus.items() if g=='Sedimentibacter'}
print()
print("   Sedimentibacter MAGs: %s" % sorted(sed_mags))
mag2bin={v:k for k,v in bin2mag.items()}
bins={mag2bin.get(m) for m in sed_mags if mag2bin.get(m)}
print("   their bins: %s" % sorted(b for b in bins if b))
best=None
for a,b,v in rows:
    ba, bb = os.path.basename(a).replace('.fa',''), os.path.basename(b).replace('.fa','')
    if ba in bins and bb in bins and ba != bb:
        if best is None or v > best[2]: best=(ba,bb,v)
if best:
    print("   SED x LIQ Sedimentibacter ANI = %.2f%%   -> %s"
          % (best[2], "SAME species" if best[2] >= 95 else "DIFFERENT species"))
else:
    print("   no Sedimentibacter-vs-Sedimentibacter row found in the ANI table")
print()
print("   in-situ ten: %d genomes, %d pair(s) >= 95%% among them"
      % (len(ten), sum(1 for a,b,v in hi if bin2mag.get(os.path.basename(a).replace('.fa','')) in ten
                                        and bin2mag.get(os.path.basename(b).replace('.fa','')) in ten)))
PYEOF

# ---------------------------------------------------------------- C. kofam
say "C. step 21 re-derived -- KOfam tables found recursively this time"
mapfile -t TBL < <(find "$KOF" -type f \( -name '*.tbl' -o -name '*.txt' -o -name '*.tsv' \) 2>/dev/null | grep -viE 'ko_list|profiles' | head -100)
echo "   candidate result files: ${#TBL[@]}"
for f in "${TBL[@]:0:5}"; do printf '     %s  (%s lines)\n' "${f#$KOF/}" "$(wc -l < "$f")"; done
[ "${#TBL[@]}" -gt 5 ] && echo "     ... and $(( ${#TBL[@]} - 5 )) more"
if [ "${#TBL[@]}" -eq 0 ]; then no "nothing under $KOF looks like KOfam output"; else
"$PY" - "$KOLIST" "$OUT" "${TBL[@]}" <<'PYEOF'
import sys, os, re, json
kolist, out = sys.argv[1], sys.argv[2]
files = sys.argv[3:]
thr={}; defn={}
with open(kolist) as fh:
    for line in fh:
        p=line.rstrip('\n').split('\t')
        if len(p)<3 or not p[0].startswith('K'): continue
        try: t=float(p[1])
        except ValueError: t=None
        thr[p[0]]=t; defn[p[0]]=p[-1]
def kos(pat):
    rx=re.compile(pat, re.I)
    return {k for k,d in defn.items() if rx.search(d or '')}
PROD=kos(r'cob\(|cobalamin|cobyrinate|cobyrinic|cobinamide|precorrin|cobaltochelatase|adenosylcobinamide|cobalt-precorrin')
CONS=kos(r'methionine synthase|methylmalonyl-CoA mutase|ribonucleoside-triphosphate reductase|glutamate mutase|ethanolamine ammonia-lyase|diol dehydratase')
TRAN=kos(r'\bbtu[BFCD]\b|vitamin B12 transport|cobalamin transport')
print("   biosynthesis %d KOs   B12-dependent %d   transport %d" % (len(PROD),len(CONS),len(TRAN)))

def parse(f):
    best={}; n=0
    for line in open(f, errors='replace'):
        if line.startswith('#'): continue
        p=line.split()
        if len(p)<6: continue
        ko=next((x for x in p[:4] if re.fullmatch(r'K\d{5}', x)), None)
        if not ko: continue
        gene=p[0] if p[0]!=ko else p[1]
        sc=None
        for x in p[3:8]:
            try:
                v=float(x)
                if 0 < v < 100000 and '.' in x: sc=v; break
            except ValueError: continue
        if sc is None: continue
        t=thr.get(ko)
        if t is not None and sc < t: continue
        n+=1
        if gene not in best or sc > best[gene][1]: best[gene]=(ko,sc)
    hits={}
    for g,(ko,sc) in best.items(): hits[ko]=hits.get(ko,0)+1
    return hits, n

rows=[]
for f in files:
    g=os.path.basename(f)
    for ext in ('.tbl','.txt','.tsv'): g=g.replace(ext,'')
    h,n=parse(f)
    if not h: continue
    rows.append(dict(genome=g, nko=len(h), prod=sorted(PROD&set(h)),
                     cons=sorted(CONS&set(h)), tran=sorted(TRAN&set(h))))
if not rows:
    print("   parsed 0 genomes -- the file format is not what this parser expects")
    print("   first 3 lines of the first candidate:")
    for i,l in enumerate(open(files[0], errors='replace')):
        if i>=3: break
        print("     %s" % l.rstrip()[:200])
    sys.exit(0)
w=max(len(r['genome']) for r in rows)
print()
print("   %-*s %6s %9s %9s %9s" % (w,'genome','KOs','biosynth','B12-dep','transport'))
for r in sorted(rows, key=lambda r:-len(r['prod'])):
    print("   %-*s %6d %9d %9d %9d" % (w,r['genome'],r['nko'],len(r['prod']),len(r['cons']),len(r['tran'])))
P=[r for r in rows if len(r['prod'])>=10]; C=[r for r in rows if r['cons']]; T=[r for r in rows if r['tran']]
print()
print("   genomes with >=10 biosynthesis KOs : %d" % len(P))
print("   genomes with a B12-DEPENDENT enzyme: %d" % len(C))
print("   genomes with a B12 transporter     : %d" % len(T))
print()
print("   step 21 said 15 producers, 0 consumers, from Bakta product strings.")
print("   >>> producers here: %d    consumers (B12-dependent enzyme): %d" % (len(P), len(C)))
for r in C[:10]:
    print("       %-26s %s" % (r['genome'], [ '%s %s' % (k,(defn.get(k) or '')[:34]) for k in r['cons'] ][:3]))
print()
print("   The '>=10 biosynthesis KOs' cut is MINE and is reported with the call.")
print("   Per-genome counts are above so the cut can be moved.")
json.dump(rows, open(os.path.join(out,'b12_kofam.json'),'w'), indent=2)
PYEOF
fi
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
