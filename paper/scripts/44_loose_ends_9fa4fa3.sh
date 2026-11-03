#!/usr/bin/env bash
# 44_loose_ends.sh -- the four open numbers, and the in-situ set decision.
#
# All read-only. Nothing is rebuilt, nothing is written outside 44_loose_ends/.
# Minutes, not hours. No tmux needed.
#
# WHAT THIS SETTLES
#   B  the in-situ set: what threshold produced "the ten", where MTP4 falls against it,
#      and therefore whether Part 1's genome set is ten or eleven
#   C  Sedimentibacter SED x LIQ ANI -- Part 1's species count
#   D  the Methanosarcina count: two in cross_abundance.tsv, three in an earlier draft
#   E  cross_abundance.tsv re-read from the file, not from the paste that lost leading 0.
#   F  step 21's "15 B12 producers, 0 consumers", re-derived from KOfam at curated
#      thresholds instead of from Bakta product strings
#
# Every section locates its own inputs and says what it found. A section that cannot find
# its input says so and does not guess -- the alternative is a confident wrong number,
# which is what put three of these four on the open list to begin with.
set -uo pipefail
WS2="${WS2:-${WS2}}"
OUT="$WS2/44_loose_ends"; mkdir -p "$OUT"
say(){ printf '\n=== %s ===\n' "$*"; }; ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }; note(){ printf '  ...    %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
PY="$WS2/envs/smetana/bin/python"; command -v "$PY" >/dev/null 2>&1 || PY=python3
MTP4=LIQ_das_vamb_003; HYD=LIQ_das_vamb_017
echo "44  $(date -Is)"; echo "out: $OUT"

# ---------------------------------------------------------------- A. inputs
say "A. locate the inputs, and say which were found"
find_one(){ # pattern... -> first existing match
  for pat in "$@"; do
    f=$(find "$WS2" -maxdepth 4 -type f -name "$pat" 2>/dev/null | head -1)
    [ -n "$f" ] && { echo "$f"; return 0; }
  done; return 1
}
TEN=$(find_one 'in_situ_10.txt' 'in_situ*.txt') || true
XAB=$(find_one 'cross_abundance.tsv') || true
GTDB=$(find_one '*gtdbtk*summary.tsv' '*gtdb*.tsv' '*taxonomy*.tsv' '*classif*.tsv') || true
ANI=$(find_one '*ani*.tsv' '*ANI*.tsv' '*fastani*' '*skani*') || true
XW=$(find_one 'MAG_CROSSWALK.tsv') || true
KOF="$WS2/11_kofam"; KOLIST=$(find "$KOF" -maxdepth 2 -name ko_list 2>/dev/null | head -1) || true
for v in TEN XAB GTDB ANI XW KOLIST; do
  eval "p=\${$v:-}"
  if [ -n "$p" ] && [ -s "$p" ]; then printf '  [OK]   %-8s %s\n' "$v" "$p"
  else printf '  [--]   %-8s not found\n' "$v"; fi
done
[ -d "$KOF" ] && ok "kofam dir $KOF ($(ls "$KOF"/*.tbl 2>/dev/null | wc -l) .tbl files)" || no "no $KOF"

# ---------------------------------------------------------------- B. in-situ set
say "B. the in-situ set -- what the threshold was, and where MTP4 falls"
if [ -n "${TEN:-}" ]; then
  echo "   the ten as recorded:"; sed 's/^/     /' "$TEN"
fi
if [ -n "${XAB:-}" ]; then
  note "sediment breadth for all MAGs, from $(basename "$XAB"):"
  "$PY" - "$XAB" "${TEN:-}" "$MTP4" "$HYD" "$OUT" <<'PYEOF'
import sys, os, csv, re
xab, tenf, mtp4, hyd, out = sys.argv[1:6]
ten=set()
if tenf and os.path.exists(tenf):
    ten={l.strip() for l in open(tenf) if l.strip()}
rows=[]
with open(xab) as fh:
    sample=fh.readline(); fh.seek(0)
    delim='\t' if '\t' in sample else ','
    rd=csv.DictReader(fh, delimiter=delim)
    cols=rd.fieldnames or []
    print("     columns: %s" % cols)
    namec=next((c for c in cols if re.search(r'mag|bin|genome|name|id', c, re.I)), cols[0])
    bcol =next((c for c in cols if re.search(r'sed.*breadth|breadth.*sed', c, re.I)), None)
    if bcol is None:
        bcol=next((c for c in cols if re.search(r'breadth', c, re.I)), None)
    print("     name column: %r   breadth column: %r" % (namec, bcol))
    for r in rd:
        try: b=float(r[bcol]) if bcol and r[bcol] not in ('','NA') else None
        except ValueError: b=None
        rows.append((r[namec], b, r))
rows=[r for r in rows if r[1] is not None]
rows.sort(key=lambda x:-x[1])
print()
print("     %-26s %8s  %s" % ('MAG','sed_breadth','in the ten?'))
for n,b,_ in rows:
    tag=''
    if n in ten: tag='  <-- IN THE TEN'
    if mtp4 in n or n in mtp4: tag+='   *** MTP4 ***'
    if hyd  in n or n in hyd:  tag+='   (the modelled hydrogenotroph)'
    print("     %-26s %8.4f%s" % (n, b, tag))
if ten:
    inb=[b for n,b,_ in rows if n in ten]
    outb=[(n,b) for n,b,_ in rows if n not in ten]
    if inb:
        print()
        print("     lowest breadth INSIDE the ten : %.4f" % min(inb))
        above=[(n,b) for n,b in outb if b >= min(inb)]
        print("     MAGs OUTSIDE the ten at or above that breadth: %d" % len(above))
        for n,b in sorted(above, key=lambda x:-x[1]):
            print("       %-26s %.4f%s" % (n, b, '   *** MTP4 ***' if mtp4 in n else ''))
        if not above:
            print("       none -- the ten are exactly the top 10 by sediment breadth,")
            print("       and MTP4 is therefore genuinely below the cut, not an oversight.")
PYEOF
else no "no cross_abundance.tsv -- section B cannot run"; fi

# ---------------------------------------------------------------- C. ANI
say "C. Sedimentibacter SED x LIQ, and every pair at or above 95% ANI"
if [ -n "${ANI:-}" ]; then
  note "from $(basename "$ANI"), header:"; head -1 "$ANI" | sed 's/^/     /'
  note "all pairs >= 95% (same species by the usual cutoff):"
  "$PY" - "$ANI" <<'PYEOF'
import sys, csv, re
p=sys.argv[1]
with open(p) as fh:
    first=fh.readline(); fh.seek(0)
    delim='\t' if '\t' in first else ','
    rd=csv.reader(fh, delimiter=delim)
    rows=list(rd)
hdr=rows[0] if rows and not re.match(r'^[\d.]+$', rows[0][-1] or '0') else None
data=rows[1:] if hdr else rows
hits=[]
for r in data:
    nums=[(i,float(x)) for i,x in enumerate(r) if re.fullmatch(r'\d+\.?\d*', (x or '').strip())]
    ani=[v for i,v in nums if 50 <= v <= 100]
    if not ani: continue
    v=max(ani)
    names=[x for x in r if not re.fullmatch(r'[\d.]+', (x or '').strip())][:2]
    if v >= 95 and len(names) >= 2 and names[0] != names[1]:
        hits.append((names[0], names[1], v))
seen=set(); uniq=[]
for a,b,v in hits:
    k=tuple(sorted([a,b]))
    if k in seen: continue
    seen.add(k); uniq.append((a,b,v))
for a,b,v in sorted(uniq, key=lambda x:-x[2]):
    print("     %-30s %-30s %6.2f" % (a,b,v))
print("     pairs >= 95%%: %d" % len(uniq))
print()
print("     rows mentioning Sedimentibacter-relevant bins:")
for r in data:
    line='\t'.join(r)
    if re.search(r'semibin_002|metabat_018', line):
        print("       %s" % line[:180])
PYEOF
else
  no "no ANI table found."
  note "step 40 recorded 'one pair >=95% in all 29 (Chakrabartyella, 99.75%)'."
  note "If that run's output is elsewhere, point WS2 at it; otherwise ANI must be re-run"
  note "and that is a compute step, not a lookup, so it gets its own script."
fi

# ---------------------------------------------------------------- D. taxonomy
say "D. the archaeal census -- how many Methanosarcina, exactly?"
if [ -n "${GTDB:-}" ]; then
  note "from $(basename "$GTDB"):"
  "$PY" - "$GTDB" "${XW:-}" <<'PYEOF'
import sys, os, csv, re
g=sys.argv[1]; xw=sys.argv[2] if len(sys.argv)>2 else ''
with open(g) as fh:
    first=fh.readline(); fh.seek(0)
    rd=csv.DictReader(fh, delimiter='\t' if '\t' in first else ',')
    cols=rd.fieldnames or []
    namec=next((c for c in cols if re.search(r'user_genome|genome|bin|mag|name|id', c, re.I)), cols[0])
    taxc =next((c for c in cols if re.search(r'classification|taxonom|lineage', c, re.I)), None)
    rows=list(rd)
print("     %d rows; name=%r taxonomy=%r" % (len(rows), namec, taxc))
if not taxc: sys.exit(0)
arch=[r for r in rows if 'Archaea' in (r[taxc] or '')]
print("     archaeal genomes: %d" % len(arch))
for r in arch:
    t=r[taxc]
    fam=re.search(r'f__([^;]*)', t); gen=re.search(r'g__([^;]*)', t); sp=re.search(r's__([^;]*)', t)
    print("       %-26s f=%-24s g=%-22s s=%s" % (r[namec],
          fam.group(1) if fam else '-', gen.group(1) if gen else '-', sp.group(1) if sp else '-'))
for key in ('Methanosarcina','Methanosarcinaceae','Methanofollis','Methanomicrobi'):
    n=[r[namec] for r in rows if key in (r[taxc] or '')]
    print("     %-20s %d  %s" % (key, len(n), n))
PYEOF
else no "no GTDB summary found -- the Methanosarcina count stays open"; fi

# ---------------------------------------------------------------- E. abundance
say "E. cross_abundance.tsv, re-read from the file"
if [ -n "${XAB:-}" ]; then
  note "the pasted copy lost leading '0.' on 18 rows. These are the file's own values."
  "$PY" - "$XAB" <<'PYEOF'
import sys, csv, re
p=sys.argv[1]
with open(p) as fh:
    first=fh.readline(); fh.seek(0)
    rd=csv.reader(fh, delimiter='\t' if '\t' in first else ',')
    rows=list(rd)
hdr=rows[0]; data=rows[1:]
w=[max(len(str(r[i])) if i<len(r) else 0 for r in rows) for i in range(len(hdr))]
print("     " + "  ".join(str(h).ljust(w[i]) for i,h in enumerate(hdr)))
for r in data:
    print("     " + "  ".join(str(r[i] if i<len(r) else '').ljust(w[i]) for i in range(len(hdr))))
print()
sus=[]
for r in data:
    for i,x in enumerate(r):
        if re.fullmatch(r'\.\d+', (x or '').strip()): sus.append((r[0], hdr[i], x))
print("     values written as '.nnn' rather than '0.nnn': %d %s" % (len(sus), sus[:5]))
PYEOF
else no "no cross_abundance.tsv"; fi

# ---------------------------------------------------------------- F. B12
say "F. step 21 re-derived -- B12 producers and consumers, from KOfam"
note "KO numbers are read out of ko_list BY DESCRIPTION, never from memory (Rule 1)."
note "hits are deduplicated to one best per gene and thresholded at ko_list's own value."
if [ -n "${KOLIST:-}" ] && [ -d "$KOF" ]; then
  "$PY" - "$KOLIST" "$KOF" "$OUT" <<'PYEOF'
import sys, os, glob, re, csv, json
kolist, kof, out = sys.argv[1:4]
# --- 1. thresholds and definitions
thr={}; defn={}
with open(kolist) as fh:
    for line in fh:
        p=line.rstrip('\n').split('\t')
        if len(p) < 3 or not p[0].startswith('K'): continue
        try: t=float(p[1])
        except ValueError: t=None
        thr[p[0]]=t; defn[p[0]]=p[-1]
print("     ko_list: %d KOs" % len(thr))

def kos(pat):
    rx=re.compile(pat, re.I)
    return {k for k,d in defn.items() if rx.search(d or '')}

PROD = kos(r'cob\(|cobalamin|cobyrinate|cobyrinic|cobinamide|precorrin|cobaltochelatase|adenosylcobinamide|cobalt-precorrin')
CONS = kos(r'methionine synthase|methylmalonyl-CoA mutase|ribonucleoside-triphosphate reductase|glutamate mutase|ethanolamine ammonia-lyase|diol dehydratase')
TRAN = kos(r'\bbtu[BFCD]\b|vitamin B12 transport|cobalamin transport')
for nm, s in (('de novo / salvage biosynthesis', PROD), ('B12-DEPENDENT enzymes', CONS), ('B12 transporters', TRAN)):
    print("     %-32s %3d KOs  e.g. %s" % (nm, len(s), sorted(s)[:5]))
print()

# --- 2. per-genome hits, deduplicated and thresholded
def parse(tbl):
    best={}
    with open(tbl) as fh:
        for line in fh:
            if line.startswith('#'): continue
            p=line.split()
            if len(p) < 6: continue
            gene, ko = p[0], p[2]
            if not ko.startswith('K'): continue
            try: score=float(p[5])
            except ValueError: continue
            t=thr.get(ko)
            if t is None or score < t: continue          # ko_list's own threshold
            if gene not in best or score > best[gene][1]: best[gene]=(ko, score)
    hit={}
    for gene,(ko,sc) in best.items(): hit.setdefault(ko,0); hit[ko]+=1
    return hit

rows=[]
for tbl in sorted(glob.glob(os.path.join(kof,'*.tbl'))):
    g=os.path.basename(tbl).replace('.tbl','')
    h=parse(tbl)
    rows.append(dict(genome=g,
                     prod=sorted(PROD & set(h)), cons=sorted(CONS & set(h)),
                     tran=sorted(TRAN & set(h)), nko=len(h)))
if not rows:
    print("     no .tbl files parsed"); sys.exit(0)
w=max(len(r['genome']) for r in rows)
print("     %-*s %6s %9s %9s %9s" % (w,'genome','KOs','biosynth','B12-dep','transport'))
for r in sorted(rows, key=lambda r:-len(r['prod'])):
    print("     %-*s %6d %9d %9d %9d" % (w, r['genome'], r['nko'], len(r['prod']), len(r['cons']), len(r['tran'])))
P=[r for r in rows if len(r['prod']) >= 10]
C=[r for r in rows if r['cons']]
T=[r for r in rows if r['tran']]
print()
print("     genomes with >=10 biosynthesis KOs (a de novo pathway) : %d" % len(P))
print("     genomes with any B12-DEPENDENT enzyme                  : %d" % len(C))
print("     genomes with a B12 transporter                         : %d" % len(T))
print()
print("     step 21 said: 15 producers, 0 consumers.")
print("     >>> 'consumer' was defined by Bakta product strings. Defined by KOfam as")
print("         'carries a B12-dependent enzyme', the count is %d." % len(C))
if C:
    print("     genomes carrying B12-dependent enzymes:")
    for r in C[:12]:
        print("       %-26s %s" % (r['genome'], [ '%s (%s)' % (k, (defn.get(k) or '')[:38]) for k in r['cons'] ][:3]))
json.dump(rows, open(os.path.join(out,'b12_kofam.json'),'w'), indent=2)
print()
print("     NOTE the threshold: '>=10 biosynthesis KOs' is mine, and it is a choice.")
print("     The per-genome counts above are the data; the cut is reported with the call.")
PYEOF
else no "no ko_list or no 11_kofam -- section F cannot run"; fi

say "DONE"; echo "written to $OUT"; echo "STEP_STATUS: COMPLETE  $(date -Is)"
