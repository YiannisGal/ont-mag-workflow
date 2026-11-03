#!/usr/bin/env bash
# 64_mtbB_mtaB.sh -- the other two methylamine/methanol markers, done right.
#
#   bash scripts/run.sh scripts/64_mtbB_mtaB_XXXXXXX.sh
#
# WHAT IS STILL BROKEN. mtbB (K16178) at 4.9x and mtaB (K04480) at 2.2x the
# sediment mcrA count were declared unreportable in the findings report and have
# never been re-measured. Both were counted on panels that lacked their own
# siblings, which is precisely what put mtmB at 1,978 instead of 18.
#
# WHAT 60-63 TAUGHT, BUILT IN HERE
#   1. Every sibling is present AS A DECOY that can win the read.
#   2. THE TARGET SET IS RESTRICTED TO ARCHAEA AT THE QUERY. My mtmB fetch used
#      a protein name that NCBI applies to the whole MttB superfamily, and seven
#      bacterial proteins went straight into the target set and held 67 sediment
#      reads. Adding "Archaea"[Organism] to every target query prevents that at
#      source instead of correcting it afterwards.
#   3. The same names restricted to BACTERIA are fetched deliberately, as a
#      DECOYbact bucket, so superfamily reads have somewhere correct to go.
#   4. Section B prints the ORGANISM of every reference, grouped by label. Rule
#      107: a protein-name query returns what the database calls that name, not
#      what you mean by it -- so look at what arrived.
#   5. mtmB is in the panel as a CROSS-CHECK. It must come back at 18 sediment
#      reads. If it does not, this run is not comparable to 63 and nothing here
#      should be quoted. That is a positive control for the whole procedure, not
#      just for the search.
#   6. mcrA from our own genomes is the ceiling, measured in the same run.
#   7. 31's archaeal control, reused verbatim.
#
# CALIBRATION, from our own catalogue rather than asserted. Copy counts above
# curated thresholds and read+gene spans (34_unbinned/length_bias.tsv):
#   mcrA  9 copies, 2,914 bp      mtbB  ? copies, 2,278 bp
#   mtmB 11 copies, 1,864 bp      mtaB  ? copies, 2,593 bp
# Section D computes the expected read ratio per marker from the copy counts this
# run finds, and compares it with what was measured. A marker landing near its
# own expectation is reportable; one far above it is still measuring a family.
#
# RUNTIME 20-60 min. Read-only except 64_mtbB_mtaB/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/64_mtbB_mtaB"; mkdir -p "$OUT"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"
exec > >(tee "$OUT/64.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; bad(){ printf '  [FAIL] %s\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }

echo "64 -- mtbB and mtaB, with the family in the panel"
echo "date : $(date -Is)"
echo "out  : $OUT"

SED_IN="$WS2/00_input/sediment_reads_clean.fastq.gz"
LIQ_IN="$WS2/00_input/liquid_reads.fastq.gz"
PROT="$WS2/06_binqc/checkm2/protein_files"
KOLIST="$WS2/36_gapclose/ko_list"
Q31="$WS2/31_readlevel/marker_queries.faa"
THREADS="${THREADS:-8}"; MINID=50; MINLEN=30
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CURL="curl -s --max-time 60"

for f in "$SED_IN" "$LIQ_IN" "$Q31"; do
  [ -s "$f" ] || { bad "missing $f"; echo "STEP_STATUS: FAILED"; exit 1; }
done
have diamond || { bad "diamond not on PATH"; echo "STEP_STATUS: FAILED"; exit 1; }

NSED=$(awk '$1=="N"{print $2; exit}' "$WS2/34_unbinned/read_len_hist.tsv" 2>/dev/null)
NSED="${NSED:-6827767}"
echo "  sediment denominator: $NSED (the clean library, as 31/34/61/63 used)"

# ================================================================ A
say "A. OUR OWN PROTEINS, BY KO, ABOVE THRESHOLD"
cat > "$OUT/panel_kos.tsv" <<'KOS'
K16178	mtbB	TARGET
K04480	mtaB	TARGET
K16176	mtmB	CROSSCHECK
K14083	DECOYmttB	DECOY
K16177	DECOYmtmC	DECOY
K00399	mcrA	CEILING
KOS
awk -F'\t' '{printf "  %-8s %-12s %s\n",$1,$2,$3}' "$OUT/panel_kos.tsv"

: > "$OUT/own.faa"
python3 - "$WS2" "$OUT" "$PROT" "$KOLIST" <<'PY'
import sys, os, glob
WS2, OUT, PROT, KOLIST = sys.argv[1:5]
thr = {}
if os.path.exists(KOLIST):
    with open(KOLIST, errors='replace') as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) > 1:
                try: thr[f[0]] = float(f[1])
                except ValueError: pass
panel = {}
for line in open(os.path.join(OUT, 'panel_kos.tsv')):
    ko, name, role = line.rstrip('\n').split('\t')
    panel[ko] = name
want, seen = {}, set()
for t in sorted(glob.glob(os.path.join(WS2, '45_kofam_full', '*.tbl')) +
                glob.glob(os.path.join(WS2, '45_kofam_full', 'chunks', '*.tbl'))):
    for line in open(t, errors='replace'):
        if line.startswith('#'): continue
        f = line.split()
        if len(f) < 6: continue
        ko = f[2]
        if ko not in panel: continue
        try: sc = float(f[5])
        except ValueError: continue
        cut = thr.get(ko)
        if cut is not None and sc < cut: continue
        tgt = f[0]
        if '|' not in tgt: continue
        b, prot = tgt.split('|', 1)
        if (b, prot, ko) in seen: continue
        seen.add((b, prot, ko))
        want.setdefault(b, {})[prot] = panel[ko]

def rdfa(p):
    n, b = None, []
    for line in open(p, errors='replace'):
        if line.startswith('>'):
            if n: yield n, ''.join(b)
            n = line[1:].split()[0]; b = []
        else: b.append(line.strip())
    if n: yield n, ''.join(b)

per = {}
with open(os.path.join(OUT, 'own.faa'), 'w') as fh:
    for b, d in want.items():
        fa = os.path.join(PROT, b + '.faa')
        if not os.path.exists(fa): continue
        for nm, sq in rdfa(fa):
            if nm in d and len(sq) >= 120:
                fh.write('>%s|OWN|%s|%s\n%s\n' % (d[nm], b, nm, sq))
                per[d[nm]] = per.get(d[nm], 0) + 1
with open(os.path.join(OUT, 'own_copies.tsv'), 'w') as fh:
    for k in sorted(per):
        print('    %-12s %d proteins above threshold' % (k, per[k]))
        fh.write('%s\t%d\n' % (k, per[k]))
PY

# ================================================================ B
say "B. REFERENCES -- TARGETS RESTRICTED TO ARCHAEA AT THE QUERY"
fetch(){   # $1 out  $2 label  $3 retmax  $4 term
  local out="$1" lab="$2" n="$3" term="$4" ids tmp
  ids=$($CURL --get --data-urlencode "term=$term" "$EUT/esearch.fcgi?db=protein&retmax=$n" 2>/dev/null \
        | grep -oP '(?<=<Id>)[0-9]+(?=</Id>)' | tr '\n' ',' | sed 's/,$//')
  if [ -z "$ids" ]; then printf '    %-12s %-62s 0\n' "$lab" "${term:0:62}"; return 1; fi
  tmp="$TMPDIR/.f64"; $CURL "$EUT/efetch.fcgi?db=protein&id=${ids}&rettype=fasta&retmode=text" > "$tmp" 2>/dev/null
  # keep the organism in the header so section B2 can show what actually arrived
  awk -v lab="$lab" '/^>/{ acc=substr($1,2); org=""
                            if (match($0,/\[[^]]*\]$/)) org=substr($0,RSTART+1,RLENGTH-2)
                            gsub(/[^A-Za-z0-9]+/,"_",org); if(org=="") org="unknown"
                            if(acc in s){sk=1;next} s[acc]=1; sk=0
                            print ">" lab "|REF|" org "|" acc; next } { if(!sk) print }' "$tmp" >> "$out"
  printf '    %-12s %-62s %s\n' "$lab" "${term:0:62}" "$(printf '%s' "$ids" | tr ',' '\n' | wc -l)"
  sleep 0.4
}

if ! $CURL -o /dev/null "$EUT/einfo.fcgi"; then
  bad "NCBI unreachable with curl"; echo "STEP_STATUS: FAILED"; exit 1
fi
ok "NCBI reachable"

if [ ! -s "$OUT/refs.faa" ]; then
  : > "$OUT/refs.faa"
  sub "targets and archaeal siblings"
  fetch "$OUT/refs.faa" mtbB      50 '"dimethylamine methyltransferase"[Protein Name] AND "Archaea"[Organism]'
  fetch "$OUT/refs.faa" mtaB      50 '"methanol--corrinoid protein Co-methyltransferase"[Protein Name] AND "Archaea"[Organism]'
  fetch "$OUT/refs.faa" mtaB      30 '"methanol methyltransferase"[Protein Name] AND "Archaea"[Organism]'
  fetch "$OUT/refs.faa" mtmB      40 '"monomethylamine methyltransferase"[Protein Name] AND "Archaea"[Organism]'
  fetch "$OUT/refs.faa" DECOYmttB 40 '"trimethylamine methyltransferase"[Protein Name] AND "Archaea"[Organism]'
  fetch "$OUT/refs.faa" DECOYmtmC 40 '"monomethylamine corrinoid protein"[Protein Name] AND "Archaea"[Organism]'
  sub "the same names restricted to BACTERIA, on purpose, as a decoy"
  fetch "$OUT/refs.faa" DECOYbact 50 '"methylamine--corrinoid protein Co-methyltransferase"[Protein Name] AND "Bacteria"[Organism]'
  fetch "$OUT/refs.faa" DECOYbact 40 '"trimethylamine methyltransferase family protein"[Protein Name] AND "Bacteria"[Organism]'
  fetch "$OUT/refs.faa" DECOYbact 30 '"dimethylamine methyltransferase"[Protein Name] AND "Bacteria"[Organism]'
fi

sub "B2. WHAT ACTUALLY ARRIVED -- organisms per label (rule 107)"
grep '^>' "$OUT/refs.faa" | awk -F'|' '{print $1"\t"$3}' | sed 's/^>//' | sort -u \
  | awk -F'\t' '{c[$1]=c[$1]" "$2; n[$1]++}
      END{ for(k in n){ printf "\n    %-12s %d organism(s):\n", k, n[k]
                        split(c[k],a," "); for(i=1;i<=length(a);i++) if(a[i]!="") printf "      %s\n", a[i] } }' \
  | head -80

echo
echo "    CHECK: every organism under mtbB, mtaB, mtmB, DECOYmttB and DECOYmtmC"
echo "    should be an archaeon. Anything that is not means the Organism filter"
echo "    did not do its job and the count built on it is not a gene count."

sub "31's archaeal control, verbatim"
awk '/^>/{ p=(index($0,"CONTROL_arc")==2); if(p) print; next } { if(p) print }' "$Q31" > "$OUT/ctrl.faa"
NC=$(grep -c '^>' "$OUT/ctrl.faa" | head -1); NC=${NC:-0}
[ "$NC" -gt 0 ] || { bad "no control"; echo "STEP_STATUS: FAILED"; exit 1; }
ok "$NC CONTROL_arc sequences"

cat "$OUT/own.faa" "$OUT/refs.faa" "$OUT/ctrl.faa" > "$OUT/panel.faa"
ok "panel: $(grep -c '^>' "$OUT/panel.faa" | head -1) sequences"
grep '^>' "$OUT/panel.faa" | awk -F'|' '{c[substr($1,2)]++} END{for(k in c) printf "    %-12s %d\n",k,c[k]}' | sort

# ================================================================ C
say "C. READ SEARCH -- 31's parameters"
diamond makedb --in "$OUT/panel.faa" -d "$OUT/panel" >/dev/null 2>&1 || {
  bad "makedb failed"; echo "STEP_STATUS: FAILED"; exit 1; }
QHASH=$(md5sum "$OUT/panel.faa" | cut -c1-12)
echo "  panel hash: $QHASH"

search(){
  local tag="$1"
  local reads="$2"
  local out="$OUT/${tag}_hits.tsv"
  local stamp="$OUT/${tag}.qhash"
  local rc
  if [ -s "$out" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$QHASH" ]; then
    echo "  [$tag] cached"; return 0; fi
  echo "  [$tag] diamond blastx --frameshift 15 --id $MINID ..."
  diamond blastx -d "$OUT/panel" -q "$reads" --frameshift 15 --id "$MINID" \
      --query-cover 0 --max-target-seqs 5 \
      -f 6 qseqid sseqid pident length evalue bitscore \
      -p "$THREADS" -o "$out.part" 2>"$OUT/${tag}.err"
  rc=$?
  if [ "$rc" -eq 0 ]; then mv "$out.part" "$out"; echo "$QHASH" > "$stamp"
  else rm -f "$out.part"; bad "$tag failed rc=$rc"; return 1; fi
}
search SED "$SED_IN" || exit 1
search LIQ "$LIQ_IN" || exit 1

# ================================================================ D
say "D. RESULT, WITH ITS OWN CALIBRATION"
python3 - "$OUT" "$NSED" "$MINLEN" <<'PY'
import sys, os
from collections import defaultdict
OUT, NSED, MINLEN = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

# read+gene span, 34_unbinned/length_bias.tsv
SPAN = {'mcrA': 2914, 'mtmB': 1864, 'mtbB': 2278, 'mtaB': 2593}
copies = {}
p = os.path.join(OUT, 'own_copies.tsv')
if os.path.exists(p):
    for line in open(p):
        k, v = line.split('\t'); copies[k] = int(v)

nref = defaultdict(int)
for line in open(os.path.join(OUT, 'panel.faa')):
    if line.startswith('>'): nref[line[1:].split('|')[0]] += 1

def load(tag):
    best = {}
    q = os.path.join(OUT, '%s_hits.tsv' % tag)
    if not os.path.exists(q): return best
    for line in open(q):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        try: ln, bs = int(f[3]), float(f[5])
        except ValueError: continue
        if ln < MINLEN: continue
        cur = best.get(f[0])
        if cur is None or bs > cur[1]: best[f[0]] = (f[1], bs)
    return best

for tag in ('SED', 'LIQ'):
    b = load(tag)
    if not b: continue
    cnt = defaultdict(int)
    for _q, (ref, _s) in b.items(): cnt[ref.split('|')[0]] += 1
    print('\n  %s' % tag)
    print('    %-14s %9s %11s %8s' % ('marker', 'reads', 'per 1M', 'n refs'))
    for mk in sorted(cnt, key=lambda k: -cnt[k]):
        pm = ('%11.2f' % (1e6 * cnt[mk] / NSED)) if tag == 'SED' else '%11s' % '-'
        print('    %-14s %9d %s %8d' % (mk, cnt[mk], pm, nref.get(mk, 0)))

    if tag == 'SED':
        x = cnt.get('mtmB', 0)
        print('\n    CROSS-CHECK: mtmB = %d (step 63 measured 18)' % x)
        if abs(x - 18) <= 3:
            print('      -> agrees. This run is comparable to 63.')
        else:
            print('      -> DOES NOT AGREE. The panel differs from 63 in a way that')
            print('         moves the counts, so nothing below is comparable to it.')

        mcra = cnt.get('mcrA', 0)
        print('\n    CEILING, calibrated from our own copy counts')
        print('      %-8s %8s %8s %10s %12s %10s' %
              ('marker', 'reads', 'copies', 'span bp', 'expected r', 'observed r'))
        for mk in ('mtmB', 'mtbB', 'mtaB'):
            c, sp = copies.get(mk), SPAN.get(mk)
            if not mcra or not c or not copies.get('mcrA'):
                print('      %-8s %8d   (no copy count -- cannot calibrate)' % (mk, cnt.get(mk, 0)))
                continue
            exp = (c * sp) / (copies['mcrA'] * SPAN['mcrA'])
            obs = cnt.get(mk, 0) / mcra
            flag = '  OK' if obs <= exp * 2.5 else '  STILL HIGH'
            print('      %-8s %8d %8d %10d %12.2f %10.2f%s'
                  % (mk, cnt.get(mk, 0), c, sp, exp, obs, flag))
        print('\n      expected r = (copies x span) / (mcrA copies x mcrA span).')
        print('      A marker within ~2.5x of its own expectation is reportable;')
        print('      far above it, the panel is still measuring a family.')
        print('\n    DETECTION FLOOR: one read = %.3f per million.' % (1e6 / NSED))
PY

say "WHAT TO WRITE"
cat <<'EOT'
  The cross-check decides whether to read the rest. mtmB must come back at 18.

  For each of mtbB and mtaB: if the observed mtmB-style ratio sits near its own
  expectation, the marker is reportable at that count with the panel stated, and
  the findings report's "unreportable" line can be replaced with a number. If it
  is still far above, the honest sentence stays: that marker measures a corrinoid
  methyltransferase family wider than the gene, and no read-level count for it
  belongs in the chapter.

  Either way section B2 is the part to actually read first. If a bacterium
  appears under mtbB, mtaB or mtmB, the Organism filter failed and the counts
  are not gene counts -- which is the whole reason this script exists.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
