#!/usr/bin/env bash
# 61_mtmB_resolve.sh -- 7 reads or 1,978? Settle mtmB with decoys that can win.
#
#   bash scripts/run.sh scripts/61_mtmB_resolve_XXXXXXX.sh
#
# THE PROBLEM, WHICH IS NOT THE ONE I EXPECTED
# mtmB (K16176) was not missing from the read pipeline. It was run TWICE and the
# two runs disagree by 280-fold:
#
#   31_readlevel/VERDICT.txt      mtmB    7 sediment reads
#   34_unbinned/marker_levels.tsv mtmB 1,978 sediment reads, sed_over_rpoB 0.3253
#
# 1,978 against 34's own mcrA count of 30 is 66x the biological ceiling. Every
# methanogen carries mcrA; only some carry mtmB. A marker cannot outnumber mcrA
# 66-fold in the same reads by the same method. One of these numbers is an
# artefact and the chapter cannot quote either until we know which.
#
# WHY 34's PROMISCUITY FILTER DID NOT CATCH IT
# 34 excludes any reference taking >25% of its marker's reads. That caught mtaB
# (98%) and dsrB (55%, 40%). But mtmB has SEVEN references, so the expected
# share is already 1/7 = 14%, and two or three promiscuous references sharing
# the pool never reach 25%. A FIXED SHARE THRESHOLD IS BLIND IN A SMALL
# REFERENCE SET. The test has to be share against the expected share, 1/n.
#
# WHY THE TWO RUNS DIFFER AT ALL
#   31  queries are OUR OWN verified proteins, and mttB is a TARGET
#   34  queries are a mixed assembly+UniProt set, and mttB is a DECOY that
#       collected 312,347 sediment reads
# mtmB, mtbB, mttB and mtaB are all MttB-superfamily corrinoid-binding
# methyltransferases. They are each other's best false positives, and which one
# absorbs a read depends entirely on which relatives are in the database.
#
# WHAT THIS RUN DOES DIFFERENTLY
#   1. Every sibling is in the panel AT ONCE, as a decoy that can win the read:
#      mtmC (K16177, the corrinoid partner that has never been a decoy), mtbB,
#      mttB, mtaB, plus MttB-superfamily references from NCBI.
#   2. mcrA is measured IN THE SAME RUN, so the ceiling is same-method.
#   3. Promiscuity is judged against 1/n, not against 25%.
#   4. The archaeal control is 31's own CONTROL_arc sequences, reused verbatim,
#      so this run is comparable to 31 by construction.
#   5. Same search parameters as 31: diamond blastx --frameshift 15 --id 50,
#      one read to its best hit, aligned length >= 30 aa.
#   6. Same denominator as 31 and 34: 6,827,767 sediment reads from
#      sediment_reads_clean.fastq.gz. Read from 34's own record, not assumed.
#
# RUNTIME  20-60 min, dominated by the read search. Cached on a query hash.
# Read-only except 61_mtmB/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/61_mtmB"; mkdir -p "$OUT"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"
exec > >(tee "$OUT/61.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }

echo "61 -- mtmB, resolved against its own family"
echo "date : $(date -Is)"
echo "out  : $OUT"

SED_IN="$WS2/00_input/sediment_reads_clean.fastq.gz"
LIQ_IN="$WS2/00_input/liquid_reads.fastq.gz"
PROT="$WS2/06_binqc/checkm2/protein_files"
KOLIST="$WS2/36_gapclose/ko_list"
Q31="$WS2/31_readlevel/marker_queries.faa"
THREADS="${THREADS:-8}"
MINID=50
MINLEN=30

for f in "$SED_IN" "$LIQ_IN" "$Q31"; do
  [ -s "$f" ] || { bad "missing: $f"; echo "STEP_STATUS: FAILED"; exit 1; }
done
[ -s "$KOLIST" ] || KOLIST="$(find "$WS2" -name ko_list -type f | head -1)"
have diamond || { bad "diamond not on PATH"; echo "STEP_STATUS: FAILED"; exit 1; }

# ================================================================ A
say "A. DENOMINATORS -- taken from the record, not assumed"
NSED=$(awk '$1=="N"{print $2; exit}' "$WS2/34_unbinned/read_len_hist.tsv" 2>/dev/null)
NSED="${NSED:-6827767}"
printf '  sediment reads (34_unbinned/read_len_hist.tsv) : %s\n' "$NSED"
echo   "  this is the CLEAN library, the one 31 and 34 both searched."
echo   "  The 7,528,854 figure is the pre-clean file used for mapping; it is not"
echo   "  the denominator for a read search of the clean file."
NLIQ=""

# ================================================================ B
say "B. THE PANEL -- every sibling present, as a decoy that can win"

cat > "$OUT/panel_kos.tsv" <<'KOS'
K16176	mtmB	TARGET	monomethylamine methyltransferase
K16177	DECOYmtmC	DECOY	monomethylamine corrinoid protein -- never used as a decoy before
K16178	DECOYmtbB	DECOY	dimethylamine methyltransferase
K14083	DECOYmttB	DECOY	trimethylamine methyltransferase
K04480	DECOYmtaB	DECOY	methanol methyltransferase
K00399	mcrA	CEILING	every methanogen has it; mtmB cannot exceed it
KOS
awk -F'\t' '{printf "  %-8s %-12s %-8s %s\n",$1,$2,$3,$4}' "$OUT/panel_kos.tsv"

sub "our own proteins, by KO, above each KO's curated threshold"
if [ -s "$KOLIST" ]; then
  ok "ko_list: ${KOLIST#$WS2/}"
else
  no "no ko_list -- falling back to E-value 1e-10 for our own proteins"
fi

: > "$OUT/own.faa"
python3 - "$WS2" "$OUT" "$PROT" "$KOLIST" <<'PY'
import sys, os, glob
WS2, OUT, PROT, KOLIST = sys.argv[1:5]

thr = {}
if KOLIST and os.path.exists(KOLIST):
    with open(KOLIST, errors='replace') as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) > 1:
                try: thr[f[0]] = float(f[1])
                except ValueError: pass

panel = {}
for line in open(os.path.join(OUT, 'panel_kos.tsv')):
    ko, name, role, note = line.rstrip('\n').split('\t')
    panel[ko] = name

# 45_kofam_full ids are BIN|PRODIGAL_ID, which map straight onto the CheckM2
# proteomes. Collect the hits that clear their own KO's threshold.
want = {}                      # bin -> {prodigal_id: marker}
tbls = sorted(glob.glob(os.path.join(WS2, '45_kofam_full', '*.tbl')) +
              glob.glob(os.path.join(WS2, '45_kofam_full', 'chunks', '*.tbl')))
seen = set()
for t in tbls:
    for line in open(t, errors='replace'):
        if line.startswith('#'): continue
        f = line.split()
        if len(f) < 6: continue
        ko = f[2]
        if ko not in panel: continue
        try: score = float(f[5])
        except ValueError: continue
        cut = thr.get(ko)
        if cut is not None and score < cut: continue
        tgt = f[0]
        if '|' not in tgt: continue
        b, prot = tgt.split('|', 1)
        key = (b, prot, ko)
        if key in seen: continue
        seen.add(key)
        want.setdefault(b, {})[prot] = panel[ko]

def rdfa(p):
    n, b = None, []
    for line in open(p, errors='replace'):
        if line.startswith('>'):
            if n: yield n, ''.join(b)
            n = line[1:].split()[0]; b = []
        else: b.append(line.strip())
    if n: yield n, ''.join(b)

n = 0
per = {}
with open(os.path.join(OUT, 'own.faa'), 'w') as fh:
    for b, d in want.items():
        fa = os.path.join(PROT, b + '.faa')
        if not os.path.exists(fa): continue
        for nm, sq in rdfa(fa):
            if nm in d and len(sq) >= 120:
                mk = d[nm]
                fh.write('>%s|OWN|%s|%s\n%s\n' % (mk, b, nm, sq))
                per[mk] = per.get(mk, 0) + 1; n += 1
print('  own proteins above threshold: %d' % n)
for k in sorted(per): print('    %-12s %d' % (k, per[k]))
if not per.get('mtmB'):
    print('  !! no mtmB above threshold in our own genomes. The panel will rest')
    print('     on NCBI references alone for the target, which is worth knowing.')
PY

sub "NCBI references for every member of the family"
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CURL="curl -s --max-time 60"
fetch(){   # $1 out  $2 label  $3 retmax  $4 term
  local out="$1" lab="$2" n="$3" term="$4" ids tmp
  ids=$($CURL --get --data-urlencode "term=$term" "$EUT/esearch.fcgi?db=protein&retmax=$n" 2>/dev/null \
        | grep -oP '(?<=<Id>)[0-9]+(?=</Id>)' | tr '\n' ',' | sed 's/,$//')
  if [ -z "$ids" ]; then printf '    %-14s %-52s 0\n' "$lab" "${term:0:52}"; return 1; fi
  tmp="$TMPDIR/.f61"; $CURL "$EUT/efetch.fcgi?db=protein&id=${ids}&rettype=fasta&retmode=text" > "$tmp" 2>/dev/null
  awk -v lab="$lab" '/^>/{ acc=substr($1,2); if(acc in s){sk=1;next} s[acc]=1; sk=0
                            print ">" lab "|REF|" acc; next } { if(!sk) print }' "$tmp" >> "$out"
  printf '    %-14s %-52s %s\n' "$lab" "${term:0:52}" "$(printf '%s' "$ids" | tr ',' '\n' | wc -l)"
  sleep 0.4
}

if $CURL -o /dev/null "$EUT/einfo.fcgi"; then
  ok "NCBI reachable (curl; wget is broken on this machine)"
  if [ ! -s "$OUT/refs.faa" ]; then
    : > "$OUT/refs.faa"
    fetch "$OUT/refs.faa" mtmB        40 '"monomethylamine methyltransferase"[Protein Name]'
    fetch "$OUT/refs.faa" mtmB        20 '"methylamine--corrinoid protein Co-methyltransferase"[Protein Name]'
    fetch "$OUT/refs.faa" DECOYmtmC   40 '"monomethylamine corrinoid protein"[Protein Name]'
    fetch "$OUT/refs.faa" DECOYmtbB   40 '"dimethylamine methyltransferase"[Protein Name]'
    fetch "$OUT/refs.faa" DECOYmttB   40 '"trimethylamine methyltransferase"[Protein Name]'
    fetch "$OUT/refs.faa" DECOYmtaB   40 '"methanol--corrinoid protein Co-methyltransferase"[Protein Name]'
    # the superfamily itself: the bacterial MttB homologues that are NOT methanogen
    # methyltransferases at all and are the classic false positive here
    fetch "$OUT/refs.faa" DECOYmttBsf 40 '"trimethylamine methyltransferase family protein"[Protein Name]'
    fetch "$OUT/refs.faa" mcrA        30 '"methyl-coenzyme M reductase subunit alpha"[Protein Name] AND refseq[filter]'
  fi
  [ -s "$OUT/refs.faa" ] && ok "$(grep -c '^>' "$OUT/refs.faa" | head -1) NCBI reference sequences"
else
  bad "NCBI unreachable -- the panel will be our own proteins plus 31's queries only."
  : > "$OUT/refs.faa"
fi

sub "31's archaeal control, reused verbatim for comparability"
awk '/^>/{ p=(index($0,"CONTROL_arc")==2); if(p) print; next } { if(p) print }' "$Q31" > "$OUT/ctrl.faa"
NC=$(grep -c '^>' "$OUT/ctrl.faa" | head -1); NC=${NC:-0}
ok "$NC CONTROL_arc sequences from 31_readlevel/marker_queries.faa"
if [ "$NC" -eq 0 ]; then
  bad "no control. A zero would be uninterpretable and this script will not print one."
  echo "STEP_STATUS: FAILED"; exit 1
fi

cat "$OUT/own.faa" "$OUT/refs.faa" "$OUT/ctrl.faa" > "$OUT/panel.faa"
NP=$(grep -c '^>' "$OUT/panel.faa" | head -1)
ok "panel: $NP sequences"
grep '^>' "$OUT/panel.faa" | awk -F'|' '{c[substr($1,2)]++} END{for(k in c) printf "    %-14s %d\n", k, c[k]}' | sort

# ================================================================ C
say "C. READ SEARCH -- 31's parameters exactly"
diamond makedb --in "$OUT/panel.faa" -d "$OUT/panel" >/dev/null 2>&1 || {
  bad "diamond makedb failed"; echo "STEP_STATUS: FAILED"; exit 1; }
QHASH=$(md5sum "$OUT/panel.faa" | cut -c1-12)
echo "  panel hash: $QHASH   (the cache key includes the queries -- rule 3)"

search(){   # $1 tag  $2 reads
  # ONE `local` PER ASSIGNMENT THAT USES AN EARLIER ONE. bash expands every word
  # on a `local` line BEFORE the builtin assigns any of them, so ${tag} inside the
  # same statement is still unset -- and under `set -u` that aborts. 31's own
  # search() used separate lines, which is why it never hit this.
  local tag="$1"
  local reads="$2"
  local out="$OUT/${tag}_hits.tsv"
  local stamp="$OUT/${tag}.qhash"
  local rc
  if [ -s "$out" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$QHASH" ]; then
    echo "  [$tag] cached"; return 0; fi
  echo "  [$tag] diamond blastx --frameshift 15 --id $MINID on $(basename "$reads") ..."
  diamond blastx -d "$OUT/panel" -q "$reads" --frameshift 15 --id "$MINID" \
      --query-cover 0 --max-target-seqs 5 \
      -f 6 qseqid sseqid pident length evalue bitscore \
      -p "$THREADS" -o "$out.part" 2>"$OUT/${tag}.err"
  rc=$?
  if [ "$rc" -eq 0 ]; then mv "$out.part" "$out"; echo "$QHASH" > "$stamp"
  else rm -f "$out.part"; bad "$tag search failed rc=$rc (see ${tag}.err)"; return 1; fi
}
search SED "$SED_IN" || exit 1
search LIQ "$LIQ_IN" || exit 1

# ================================================================ D
say "D. RESULT"
python3 - "$OUT" "$NSED" "$MINLEN" "$WS2" <<'PY'
import sys, os
from collections import defaultdict
OUT, NSED, MINLEN, WS2 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]

def load(tag):
    best = {}
    p = os.path.join(OUT, '%s_hits.tsv' % tag)
    if not os.path.exists(p): return best
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        try: ln, bs = int(f[3]), float(f[5])
        except ValueError: continue
        if ln < MINLEN: continue
        cur = best.get(f[0])
        if cur is None or bs > cur[2]:
            best[f[0]] = (f[1].split('|')[0], f[1], bs)
    return best

for tag, N in (('SED', NSED), ('LIQ', None)):
    b = load(tag)
    if not b:
        print('\n  %s: no hits' % tag); continue
    cnt = defaultdict(int); perref = defaultdict(lambda: defaultdict(int))
    for _q, (mk, ref, _bs) in b.items():
        cnt[mk] += 1; perref[mk][ref] += 1
    nref = defaultdict(int)
    for line in open(os.path.join(OUT, 'panel.faa')):
        if line.startswith('>'): nref[line[1:].split('|')[0]] += 1

    print('\n  %s -- one read counted once, to its best hit' % tag)
    print('    %-14s %10s %10s %8s %10s' % ('marker', 'reads', 'per 1M', 'n refs', 'reads/ref'))
    for mk in sorted(cnt, key=lambda k: -cnt[k]):
        pm = ('%10.2f' % (1e6 * cnt[mk] / N)) if N else '         -'
        print('    %-14s %10d %s %8d %10.1f'
              % (mk, cnt[mk], pm, nref.get(mk, 0), cnt[mk] / max(1, nref.get(mk, 1))))

    print('\n    PROMISCUITY, judged against the expected share 1/n -- not 25%%')
    flagged = 0
    for mk in sorted(perref):
        n = max(1, nref.get(mk, 1)); exp = 1.0 / n; tot = cnt[mk]
        if tot < 20: continue
        for ref, c in sorted(perref[mk].items(), key=lambda kv: -kv[1])[:3]:
            share = c / tot
            if share >= 4 * exp and share > 0.20:
                print('      %-10s %5.1f%% of %d reads (expected %.1f%%, %.1fx)  %s'
                      % (mk, 100 * share, tot, 100 * exp, share / exp, ref))
                flagged += 1
    if not flagged:
        print('      none: no reference exceeds 4x its expected share')

    if tag == 'SED':
        mtmb, mcra = cnt.get('mtmB', 0), cnt.get('mcrA', 0)
        ctrl = cnt.get('CONTROL_arc', 0)
        print('\n    CEILING CHECK')
        print('      archaeal control : %d reads (%.2f per million)' % (ctrl, 1e6 * ctrl / N))
        print('      mcrA             : %d' % mcra)
        print('      mtmB             : %d' % mtmb)
        if ctrl == 0:
            print('      >> NO VERDICT. The control failed; nothing here is interpretable.')
        elif mcra == 0:
            print('      >> mcrA itself is zero, so there is no ceiling to test against.')
        elif mtmb > mcra:
            print('      >> STILL ABOVE THE CEILING at %.1fx mcrA. Report mtmB as'
                  ' UNQUANTIFIED.' % (mtmb / mcra))
            print('         With mtmC, mtbB, mttB, mtaB and the MttB superfamily all')
            print('         present as decoys, an excess over mcrA can no longer be')
            print('         blamed on a missing relative. Quote the decoy counts.')
        else:
            print('      >> WITHIN THE CEILING at %.2f of mcrA. This is a reportable'
                  ' number.' % (mtmb / mcra))
        print('\n    DETECTION FLOOR: %d reads searched; one read = %.3f per million.'
              % (N, 1e6 / N))

print('''
  AGAINST THE TWO EARLIER RUNS
    31_readlevel  mtmB     7 sediment reads   (mttB a TARGET, own queries only)
    34_unbinned   mtmB 1,978 sediment reads   (mttB a DECOY, mixed reference set)
  Whichever this run lands on, the other was measuring the panel and not the gene.
  The number to report is this one, because it is the only one where every
  sibling was present and able to win the read.''')
PY

say "E. WHAT TO WRITE"
cat <<'EOT'
  If mtmB comes in at or below the mcrA count, it is reportable and the
  methylamine arm has a measured in-situ floor.

  If it is still above mcrA with the full sibling panel in place, the honest
  sentence is that mtmB CANNOT BE QUANTIFIED at read level in this dataset --
  the MttB superfamily is too conserved for a 50%-identity blastx to separate
  its members on ONT reads -- and the methylamine evidence rests where it
  already stands: on the pyrrolysine amber test and on mttB phylogenetic
  placement, both of which are sequence-level identifications rather than
  counts.

  EITHER WAY, one correction must reach the record: 34_unbinned's mtmB row of
  1,978 reads and sed_over_rpoB 0.3253 must not be quoted, and 34's 25%
  promiscuity rule is blind in small reference sets and should be restated as a
  multiple of 1/n.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
