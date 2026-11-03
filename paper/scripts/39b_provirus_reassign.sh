#!/usr/bin/env bash
# =============================================================================
# 39b_provirus_reassign.sh  --  CASHING IN 39A: the MAG-level provirus calls
#
# WHAT 39A ESTABLISHED
#   `18_viruses/all_mags_contigs.fna` holds 2,549 headers and every one of them
#   is of the form
#
#       <bin_id>~<contig_name>          e.g. LIQ_das_metabat_001~LIQ_contig_2376
#
#   So the file was never "renamed somewhere in step 18" in any mysterious way:
#   it was concatenated with a `bin~contig` namespace, which is the standard
#   DAS_Tool / dRep convention. **The MAG is written into the header itself.**
#   Step 18 looked up bare contig names against that namespace, matched nothing,
#   and recorded "0 of 38 proviruses assigned to a MAG" plus "254 ambiguous
#   contig names discarded". That was a name-space failure and nothing else.
#
#   39A named the weakest rule that fits ("a contig is a substring of query",
#   100%). It fits, but it is AMBIGUOUS by construction — `contig_1` is inside
#   `contig_11` — so it must not be used to assign anything. The delimiter rule
#   is exact, and this script uses only that.
#
# WHAT THIS DOES
#   Splits each header on the LAST `~`, verifies the left half is a real bin id
#   and the right half a real contig of that bin, and re-runs step 18's
#   assignment. Anything that does not verify is REPORTED, never guessed.
#
# THREE CONTROLS, because the last provirus count was wrong three times
#   1. ROUND TRIP. Every one of the 2,549 headers must resolve to a bin that
#      exists in the crosswalk AND to a contig that exists in that bin's FASTA.
#      A header that resolves to a bin but not to its contig means the delimiter
#      rule is wrong, not that a contig is missing.
#   2. IDENTITY, NOT ROWS. geNomad writes one row per provirus in EVERY stage
#      table, and step 25 summed eight of them into a sevenfold overcount. This
#      counts DISTINCT provirus identifiers and prints raw rows beside distinct
#      ones so a mismatch is visible instead of silently summed.
#   3. AGREEMENT WITH STEP 18. The assembly-level run produced 68 proviruses
#      across 18 MAGs via `provirus_host_links.tsv`. The MAG-level run is a
#      DIFFERENT run on a different input, so the counts need not match — but
#      the SET OF MAGS should overlap heavily. If a MAG-level provirus lands in
#      a genome the host-links file calls provirus-free, that disagreement is
#      the result, and it is printed rather than averaged away.
#
# WHAT IT DOES NOT DO
#   It does not overwrite `MAG_CATALOGUE.tsv`. `provirus_host_links.tsv` is
#   still authoritative for `n_proviruses` (step 25 settled that). This produces
#   an INDEPENDENT second measurement and says whether the two agree.
#
# OUTPUT  39_loose_ends/provirus/{REASSIGN.txt,contig_map.tsv,mag_proviruses.tsv}
# RUNTIME under a minute
# USAGE   ws2 && bash scripts/39b_provirus_reassign.sh 2>&1 | tee logs/39b.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="39b_provirus"
F="$WS2/07_mags"; V="$WS2/18_viruses"
D="$WS2/39_loose_ends/provirus"; mkdir -p "$D"

echo "######## 39b — RE-ASSIGNING THE MAG-LEVEL PROVIRUSES ########"
echo "  started: $(date -Is)"
rm -f "$D/REASSIGN.txt" "$D/contig_map.tsv" "$D/mag_proviruses.tsv"

AMC="$(find "$WS2" -maxdepth 3 -name 'all_mags_contigs.fna' 2>/dev/null | sort | head -1)"
echo "  all_mags_contigs.fna : ${AMC:-<not found>}"
echo

python3 - "$D" "$F" "$V" "$WS2" "${AMC:-}" <<'PY'
import sys, os, re, glob, csv
from collections import Counter, defaultdict

D, F, V, WS2, AMC = sys.argv[1:6]
W = 100
rep = []
STATUS = ['INCOMPLETE']


def say(s=''):
    print(s)
    rep.append(s)


def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t)
        say('=' * W)


def dump(status):
    STATUS[0] = status
    with open(os.path.join(D, 'REASSIGN.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'REASSIGN.txt')}")
    sys.exit(0)


def heads(path, cap=500000):
    out = []
    try:
        with open(path, errors='replace') as fh:
            for line in fh:
                if line.startswith('>'):
                    out.append(line[1:].rstrip('\n'))
                    if len(out) >= cap:
                        break
    except Exception as e:
        say(f"  !! cannot read {path}: {e}")
    return out


rule('1. THE DELIMITER RULE, VERIFIED BOTH WAYS')
if not AMC or not os.path.exists(AMC):
    say("  all_mags_contigs.fna not found. Nothing to re-assign.")
    dump('NO INPUT')

# ---- the truth: which contigs belong to which bin
bins = {}
fas = []
for sub in ('bins', 'final', 'MAGs', ''):
    hits = []
    for e in ('fa', 'fna', 'fasta'):
        hits += glob.glob(os.path.join(F, sub, '*.' + e) if sub
                          else os.path.join(F, '*.' + e))
    hits = [h for h in hits if 'all_mags' not in os.path.basename(h)
            and not os.path.basename(h).startswith('all_')]
    if hits:
        fas = sorted(set(hits))
        break
for fa in fas:
    b = re.sub(r'\.(fa|fna|fasta)$', '', os.path.basename(fa))
    bins[b] = {h.split()[0] if h.split() else h for h in heads(fa)}
say(f"  bin FASTAs on record : {len(bins)}  "
    f"({sum(len(v) for v in bins.values()):,} contigs)")
if not bins:
    say("  !! no per-MAG FASTAs found — cannot verify anything.")
    dump('NO REFERENCE')

q = heads(AMC)
say(f"  headers in all_mags_contigs.fna : {len(q):,}")


def split_hdr(h):
    """split on the LAST '~': a contig name may not contain one, a bin id may."""
    h = h.split()[0] if h.split() else h
    if '~' not in h:
        return None, None
    b, _, c = h.rpartition('~')
    return b, c


ok, bad_bin, bad_contig, no_delim = [], [], [], []
for h in q:
    b, c = split_hdr(h)
    if b is None:
        no_delim.append(h)
    elif b not in bins:
        bad_bin.append(h)
    elif c not in bins[b]:
        bad_contig.append(h)
    else:
        ok.append((b, c))

say()
say(f"  {'outcome':<48}{'headers':>10}{'%':>10}")
say('-' * W)
for lbl, n in (('resolved: bin exists AND contig is in that bin', len(ok)),
               ("no '~' in the header", len(no_delim)),
               ('bin id not in the catalogue', len(bad_bin)),
               ('contig not found in the named bin', len(bad_contig))):
    say(f"  {lbl:<48}{n:>10,}{100.0 * n / max(1, len(q)):>9.1f}%")
say()

pct = 100.0 * len(ok) / max(1, len(q))
# CONTROL 1 — round trip
if pct < 99.0:
    say(f"  [FAIL] CONTROL 1 (round trip): only {pct:.1f}% of headers resolve to")
    say("  a real contig of a real bin. The delimiter rule is not the whole")
    say("  story and nothing may be assigned from it. Examples:")
    for h in (bad_contig or bad_bin or no_delim)[:6]:
        say(f"    {h[:84]}")
    dump('CONTROL 1 FAILED')
say(f"  [PASS] CONTROL 1 (round trip): {pct:.1f}% of headers resolve to a real")
say("  contig of a real bin. The rule is `header.rpartition('~')`, and it is")
say("  EXACT — no substring matching, no arbitration, nothing discarded as")
say("  ambiguous. Step 18's 254 'ambiguous contig names' were an artefact of")
say("  looking up a bare contig name in a bin~contig namespace.")
say()

c2m = {}
dup = 0
for b, c in ok:
    if c in c2m and c2m[c] != b:
        dup += 1
    c2m.setdefault(c, b)
with open(os.path.join(D, 'contig_map.tsv'), 'w') as fh:
    fh.write('contig\tbin_id\n')
    for c, b in sorted(c2m.items()):
        fh.write(f'{c}\t{b}\n')
say(f"  contig -> bin map written: {len(c2m):,} contigs, {dup} name collisions")
if dup:
    say("  NOTE: a bare contig name occurring in more than one bin is exactly why")
    say("  step 25 discarded ambiguous names. The `bin~contig` header removes the")
    say("  ambiguity, so those collisions do not affect this assignment.")
say()

# =====================================================================
rule('2. THE PROVIRUS CALLS, COUNTED BY IDENTITY')
say("  geNomad writes one row per provirus in EVERY stage table. Step 25 summed")
say("  eight of them and produced a sevenfold overcount that passed its own")
say("  audit, because multiplying every count by seven does not change which")
say("  counts are non-zero. So: ONE file, and distinct identifiers.")
say()
cands = sorted(glob.glob(os.path.join(V, '**', 'all_mags_contigs*virus_summary.tsv'),
                         recursive=True))
cands = [c for c in cands if 'provirus_genes' not in c]
if not cands:
    say("  !! no all_mags_contigs*_virus_summary.tsv under 18_viruses/.")
    say("  The MAG-level geNomad output is not where this expects it; the")
    say("  contig map above is still written and usable.")
    dump('NO GENOMAD TABLE')
src = cands[0]
say(f"  reading ONE table: {os.path.relpath(src, WS2)}")
if len(cands) > 1:
    say(f"  ({len(cands)} matched; the others are NOT summed — that was the bug)")
    for c in cands[1:]:
        say(f"     not used: {os.path.relpath(c, WS2)}")

rows = []
with open(src, errors='replace') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        rows.append(r)
if not rows:
    say("  !! table is empty.")
    dump('EMPTY TABLE')

cols = {k.lower(): k for k in rows[0] if k}
topo = cols.get('topology')
seqc = cols.get('seq_name') or cols.get('contig') or list(rows[0])[0]
prov = [r for r in rows if not topo or 'provirus' in str(r.get(topo, '')).lower()]
say(f"  rows in table            : {len(rows):,}")
say(f"  rows with topology=Provirus : {len(prov):,}"
    + ('' if topo else '   (no topology column — all rows taken)'))
ids = sorted({str(r.get(seqc, '')).strip() for r in prov if str(r.get(seqc, '')).strip()})
say(f"  DISTINCT provirus identifiers : {len(ids):,}")
if len(ids) != len(prov):
    say(f"  !! {len(prov) - len(ids)} duplicate identifier(s) collapsed — that is")
    say("     the sevenfold-overcount signature. Investigate before reporting.")
say()

# a provirus id is <bin>~<contig>|provirus_<start>_<end> or similar
per_mag = Counter()
unres = []
for i in ids:
    base = i.split('|')[0]
    b, c = split_hdr(base)
    if b in bins and c in bins.get(b, ()):
        per_mag[b] += 1
    else:
        unres.append(i)
say(f"  proviruses assigned to a MAG : {sum(per_mag.values()):,} "
    f"across {len(per_mag)} genome(s)")
say(f"  unresolved                   : {len(unres)}")
for i in unres[:6]:
    say(f"    {i[:84]}")
with open(os.path.join(D, 'mag_proviruses.tsv'), 'w') as fh:
    fh.write('bin_id\tn_proviruses\n')
    for b, n in sorted(per_mag.items(), key=lambda kv: (-kv[1], kv[0])):
        fh.write(f'{b}\t{n}\n')
say()
if per_mag:
    say(f"  {'bin_id':<32}{'proviruses':>12}")
    say('-' * W)
    for b, n in sorted(per_mag.items(), key=lambda kv: (-kv[1], kv[0]))[:12]:
        say(f"  {b:<32}{n:>12}")
    if len(per_mag) > 12:
        say(f"  ... and {len(per_mag) - 12} more genomes")
say()

# =====================================================================
rule('3. CONTROL — DOES IT AGREE WITH STEP 18?')
hl = os.path.join(V, 'provirus_host_links.tsv')
if not os.path.exists(hl):
    say("  provirus_host_links.tsv not found — the independent comparison cannot")
    say("  be made. This result therefore stands ALONE and should be reported")
    say("  as a single measurement, not a confirmed one.")
    dump('OK (UNCOMPARED)')

hrows = []
with open(hl, errors='replace') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        hrows.append(r)
hmags = Counter()
for r in hrows:
    for v in r.values():
        v = str(v or '').strip()
        if v in bins:
            hmags[v] += 1
            break
say(f"  provirus_host_links.tsv : {len(hrows)} rows, "
    f"{len(hmags)} genome(s) with at least one provirus")
say(f"  this run (MAG-level)    : {sum(per_mag.values())} proviruses, "
    f"{len(per_mag)} genome(s)")
say()
say("  These are TWO DIFFERENT geNomad runs on two different inputs — the")
say("  assemblies and the concatenated MAG contigs — so the counts need not")
say("  match. The genome SETS should overlap. What matters is disagreement.")
say()
A_, B_ = set(per_mag), set(hmags)
say(f"  {'genomes in both':<40}{len(A_ & B_):>8}")
say(f"  {'MAG-level only':<40}{len(A_ - B_):>8}")
say(f"  {'assembly-level (host links) only':<40}{len(B_ - A_):>8}")
say()
if A_ and B_:
    j = len(A_ & B_) / len(A_ | B_)
    say(f"  Jaccard on the genome sets : {j:.2f}")
    if A_ - B_:
        say("  Genomes the MAG-level run calls provirus-bearing and the host-links")
        say("  file does not — each is a genuine disagreement to resolve, not an")
        say("  addition to the count:")
        for b in sorted(A_ - B_)[:8]:
            say(f"    {b}  ({per_mag[b]} provirus(es))")
    say()
say("  WHAT TO REPORT. `provirus_host_links.tsv` remains authoritative for")
say("  `n_proviruses` in the catalogue — step 25 settled that and this does not")
say("  reopen it. The value here is that the MAG-level run is no longer a")
say("  discarded file: it is a SECOND, independent measurement, and the")
say("  agreement (or the named disagreement) above is what should be quoted")
say("  when the provirus numbers are defended.")
dump('OK')
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "analysis completed" "$RC" || OK=1
check "report written" "$([ -s "$D/REASSIGN.txt" ] && echo 0 || echo 1)" || OK=1
check "contig map written" "$([ -s "$D/contig_map.tsv" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$D/REASSIGN.txt" 2>/dev/null | cut -d' ' -f2-)
echo "  status line: ${ST:-<none>}"
case "${ST:-MISSING}" in
  OK|"OK (UNCOMPARED)") : ;;
  *) echo "  [STOP] ${ST:-no status line}"; OK=1 ;;
esac
[ "$OK" -eq 0 ] && runlog "39b_provirus complete" || runlog "**39b_provirus FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
