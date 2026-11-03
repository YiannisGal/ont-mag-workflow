#!/usr/bin/env bash
# 44c_b12_properly.sh -- fix the parser, then decide whether B12 can be answered at all.
#
# 44b SETTLED and these are closed:
#   - Sedimentibacter SED x LIQ ANI = 84.18%  -> DIFFERENT species
#     (SED_das_semibin_002 / ETO_SED_MAG_001  vs  LIQ_dassub_semibin_017 / ETO_LIQ_MAG_025)
#   - the only pair >=95% anywhere is Chakrabartyella at 99.75%
#     (SED_das_vamb_004 / ETO_SED_MAG_002  vs  LIQ_das_semibin_004 / ETO_LIQ_MAG_012),
#     both inside the ten. "2 pairs" was one pair counted in both directions.
#     >>> THE IN-SITU TEN ARE NINE SPECIES.
#   - 8 archaea, 21 bacteria, 29 total.
#
# TWO BUGS IN 44b, both mine:
#
#   1. THE Methanosarcina COUNT. I matched 'Methanosarcina' as a SUBSTRING, which also
#      matches 'Methanosarcinaceae'. By genus: Methanosarcina = LIQ_das_metabat_030 and
#      LIQ_das_vamb_005 -> TWO. LIQ_das_vamb_003's genus is 'MTP4'. Methanosarcinaceae
#      (family) = THREE. The script printed 3 and then explained in its own note why 3
#      was wrong. This is the h2 -> h2o error of 41c again.
#
#   2. THE SCORE PARSER. hmmsearch --tblout is:
#        0 target  1 acc  2 query(KO)  3 acc  4 E-value  5 score  6 bias  ...
#      I scanned p[3:8] for the first float containing a dot and called it the score, so
#      for "1.6e-226  750.1" it took 1.6e-226 AS THE SCORE, compared it to a threshold of
#      ~200, and discarded a real hit. That is why every genome showed 1-13 KOs.
#
# AND A PRIOR QUESTION THE NUMBERS RAISE. The files are in 11_kofam/hits/ at 100-417
# lines. A genome-wide KOfam run on a bacterial genome yields hundreds to thousands of
# hits. If these are a TARGETED marker search, then "0 B12 producers" is a statement
# about which KOs were searched, not about the genomes -- Rule 35. Section B decides
# that from the data and REFUSES to report B12 numbers if the input is not genome-wide.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/44_loose_ends"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; command -v "$PY" >/dev/null 2>&1 || PY=python3
KOF="$WS2/11_kofam"; KOLIST="$KOF/ko_list"
echo "44c  $(date -Is)"

# ---------------------------------------------------------------- A. format
say "A. what format are these files, exactly?"
F1=$(find "$KOF" -type f -name '*.tbl' 2>/dev/null | head -1)
[ -n "$F1" ] || { bad "no .tbl under $KOF"; exit 1; }
echo "   $F1"
echo "   --- first 6 non-comment lines, with column indices ---"
grep -v '^#' "$F1" | head -6 | "$PY" -c "
import sys
for line in sys.stdin:
    p=line.split()
    print('     ' + '  '.join('%d:%s' % (i,x) for i,x in enumerate(p[:9])))
"
echo "   --- comment lines (hmmer writes its column names here) ---"
grep '^#' "$F1" | head -4 | cut -c1-190 | sed 's/^/     /'
note "how many .tbl files, and how many lines each:"
find "$KOF" -type f -name '*.tbl' | wc -l | sed 's/^/     files: /'
find "$KOF" -type f -name '*.tbl' -exec wc -l {} + | sort -n | tail -4 | sed 's/^/     /'

# ---------------------------------------------------------------- B. gate
say "B. is this a GENOME-WIDE KOfam run, or a targeted marker search?"
note "parsing with score = column 5 (0-based), the hmmsearch --tblout position."
note "a genome-wide run on a bacterial MAG gives hundreds of distinct KOs."
"$PY" - "$KOLIST" "$KOF" "$OUT" <<'PYEOF'
import sys, os, re, glob, json, statistics
kolist, kof, out = sys.argv[1:4]
thr={}; defn={}
with open(kolist) as fh:
    for line in fh:
        p=line.rstrip('\n').split('\t')
        if len(p)<3 or not p[0].startswith('K'): continue
        try: t=float(p[1])
        except ValueError: t=None
        thr[p[0]]=t; defn[p[0]]=p[-1]
print("   ko_list: %d KOs, %d with a numeric threshold" % (len(thr), sum(1 for v in thr.values() if v is not None)))

def parse(f):
    """hmmsearch --tblout: 0 target 1 acc 2 query 3 acc 4 E-value 5 score 6 bias"""
    best={}; raw=0; nothr=0
    for line in open(f, errors='replace'):
        if line.startswith('#'): continue
        p=line.split()
        if len(p) < 6: continue
        gene, ko = p[0], p[2]
        if not re.fullmatch(r'K\d{5}', ko):
            ko = next((x for x in p[:4] if re.fullmatch(r'K\d{5}', x)), None)
            if not ko: continue
        try: score=float(p[5])
        except ValueError: continue
        raw += 1
        t=thr.get(ko)
        if t is None: nothr += 1
        elif score < t: continue
        if gene not in best or score > best[gene][1]: best[gene]=(ko, score)
    kos={}
    for g,(ko,sc) in best.items(): kos[ko]=kos.get(ko,0)+1
    return kos, raw, nothr

rows=[]
for f in sorted(glob.glob(os.path.join(kof,'**','*.tbl'), recursive=True)):
    g=os.path.basename(f).replace('.tbl','')
    kos, raw, nothr = parse(f)
    rows.append((g, len(kos), raw, nothr))
if not rows: print("   no files parsed"); sys.exit(1)
w=max(len(r[0]) for r in rows)
print()
print("   %-*s %10s %10s %10s" % (w,'genome','distinct KO','raw rows','no thrshld'))
for g,n,raw,nothr in sorted(rows, key=lambda r:-r[1]):
    print("   %-*s %10d %10d %10d" % (w,g,n,raw,nothr))
med=statistics.median(r[1] for r in rows)
print()
print("   median distinct KOs per genome: %.0f" % med)
json.dump([{'genome':g,'kos':n,'raw':r,'nothr':t} for g,n,r,t in rows],
          open(os.path.join(out,'kofam_inventory.json'),'w'), indent=2)
if med < 200:
    print()
    print("   >>> NOT GENOME-WIDE. %d KOs per genome is a targeted marker search." % med)
    print("   >>> A B12 count from this input would be a statement about which KOs were")
    print("       searched, not about the genomes. Rule 35. Section C will NOT run.")
    print("GENOMEWIDE: NO")
    sys.exit(3)
print("   >>> genome-wide. Section C can proceed.")
print("GENOMEWIDE: YES")
PYEOF
GW=$?

# ---------------------------------------------------------------- C. B12
if [ "$GW" -eq 0 ]; then
  say "C. B12, re-derived genome-wide"
  "$PY" - "$KOLIST" "$KOF" "$OUT" <<'PYEOF'
import sys, os, re, glob, json
kolist, kof, out = sys.argv[1:4]
thr={}; defn={}
with open(kolist) as fh:
    for line in fh:
        p=line.rstrip('\n').split('\t')
        if len(p)<3 or not p[0].startswith('K'): continue
        try: t=float(p[1])
        except ValueError: t=None
        thr[p[0]]=t; defn[p[0]]=p[-1]
def kos_matching(pat):
    rx=re.compile(pat, re.I)
    return {k for k,d in defn.items() if rx.search(d or '')}
PROD=kos_matching(r'cob\(|cobalamin|cobyrinate|cobyrinic|cobinamide|precorrin|cobaltochelatase|adenosylcobinamide|cobalt-precorrin')
CONS=kos_matching(r'methionine synthase|methylmalonyl-CoA mutase|ribonucleoside-triphosphate reductase|glutamate mutase|ethanolamine ammonia-lyase|diol dehydratase')
TRAN=kos_matching(r'\bbtu[BFCD]\b|vitamin B12 transport|cobalamin transport')
print("   biosynthesis %d  B12-dependent %d  transport %d KOs" % (len(PROD),len(CONS),len(TRAN)))
def parse(f):
    best={}
    for line in open(f, errors='replace'):
        if line.startswith('#'): continue
        p=line.split()
        if len(p)<6: continue
        gene, ko = p[0], p[2]
        if not re.fullmatch(r'K\d{5}', ko): continue
        try: score=float(p[5])
        except ValueError: continue
        t=thr.get(ko)
        if t is not None and score < t: continue
        if gene not in best or score > best[gene][1]: best[gene]=(ko,score)
    s=set()
    for g,(ko,sc) in best.items(): s.add(ko)
    return s
rows=[]
for f in sorted(glob.glob(os.path.join(kof,'**','*.tbl'), recursive=True)):
    g=os.path.basename(f).replace('.tbl','')
    h=parse(f)
    rows.append(dict(genome=g, nko=len(h), prod=sorted(PROD&h), cons=sorted(CONS&h), tran=sorted(TRAN&h)))
w=max(len(r['genome']) for r in rows)
print()
print("   %-*s %6s %9s %9s %9s" % (w,'genome','KOs','biosynth','B12-dep','transport'))
for r in sorted(rows, key=lambda r:-len(r['prod'])):
    print("   %-*s %6d %9d %9d %9d" % (w,r['genome'],r['nko'],len(r['prod']),len(r['cons']),len(r['tran'])))
for cut in (5,10,15,20):
    P=[r for r in rows if len(r['prod'])>=cut]
    print("   producers at >=%2d biosynthesis KOs : %d" % (cut, len(P)))
C=[r for r in rows if r['cons']]; T=[r for r in rows if r['tran']]
print("   genomes with a B12-DEPENDENT enzyme: %d" % len(C))
print("   genomes with a B12 transporter     : %d" % len(T))
print()
print("   step 21 said 15 producers, 0 consumers (Bakta product strings).")
print("   The cut is reported WITH the call -- four cuts printed above, pick and state one.")
json.dump(rows, open(os.path.join(out,'b12_kofam.json'),'w'), indent=2)
PYEOF
else
  say "C. NOT RUN"
  note "the KOfam input is not genome-wide, so B12 cannot be answered from it."
  note "What would be needed: a genome-wide KOfam/hmmsearch of all 29 MAG proteomes"
  note "against the full KOfam profile set, then this same thresholded, deduplicated"
  note "parse. That is a compute step measured in hours, not a lookup."
  note "Until then step 21's '15 producers, 0 consumers' stands as UNVERIFIED, and the"
  note "Volume II caveat (the enrichment medium supplies B12 at 0.05 mg/L) is the"
  note "operative statement for the cultures."
fi

# ---------------------------------------------------------------- D. taxonomy
say "D. the archaeal counts, by GENUS not by substring"
"$PY" - "$WS2" <<'PYEOF'
import sys, os, csv, re, glob
ws2=sys.argv[1]
rows={}
for f in glob.glob(os.path.join(ws2,'08_taxonomy','**','gtdbtk.*.summary.tsv'), recursive=True):
    with open(f) as fh:
        rd=csv.DictReader(fh, delimiter='\t')
        c=rd.fieldnames or []
        n=next((x for x in c if re.search(r'user_genome|genome|bin',x,re.I)), c[0])
        t=next((x for x in c if re.search(r'classification|taxonom',x,re.I)), None)
        for r in rd:
            if t: rows[r[n]]=r[t]
def rank(tax, letter):
    m=re.search(r'%s__([^;]*)' % letter, tax)
    return (m.group(1) or '').strip() if m else ''
from collections import Counter
gen=Counter(); fam=Counter()
for k,t in rows.items():
    if 'd__Archaea' not in t: continue
    gen[rank(t,'g') or '(unnamed)'] += 1
    fam[rank(t,'f') or '(unnamed)'] += 1
print("   archaeal GENERA (exact):")
for g,n in gen.most_common(): print("     %-24s %d" % (g,n))
print("   archaeal FAMILIES (exact):")
for f,n in fam.most_common(): print("     %-24s %d" % (f,n))
print()
print("   >>> Methanosarcina (genus)      : %d" % gen.get('Methanosarcina',0))
print("   >>> Methanosarcinaceae (family) : %d" % fam.get('Methanosarcinaceae',0))
print("   >>> the difference is MTP4, whose genus is unnamed. Two counts, both correct,")
print("       of two different things. State which one the sentence means.")
PYEOF
say "DONE"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
