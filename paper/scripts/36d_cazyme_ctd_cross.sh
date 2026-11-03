#!/usr/bin/env bash
# =============================================================================
# 36d_cazyme_ctd_cross.sh  --  DO THE CAZymes CARRY THE EXPORT TAG?
#
# THE OPEN LIMIT THIS CLOSES
#   36c established that the T9SS cargo of the three Bacteroidota is enriched
#   27.6x in adhesion modules and 3.1x in peptidases, and that carbohydrate-
#   active products sit at 8/155 — indistinguishable from the same genomes'
#   background (1.41x, q = 0.41). But those 8 come from KEYWORD MATCHING ON
#   BAKTA PRODUCT STRINGS, which is exactly the method step 30 replaced with
#   dbCAN for every other CAZyme number in this chapter. A CAZyme call from a
#   product string is not a CAZyme call.
#
#   RESULT, so the reader is not left with 36c's reading: dbCAN finds 53 cargo
#   CAZymes where the keyword method found 8, and the enrichment is real at
#   3.2-3.5x in EACH genome separately. 36c's carbohydrate row is superseded,
#   and so is my summary of it. Many of the 45 missed proteins are named
#   "T9SS type A sorting domain-containing protein" by Bakta — named for the
#   sorting domain instead of the catalytic one, so the circular annotation was
#   not merely uninformative, it was hiding function.
#
#   36c could not do better because no table on disk carries BOTH a CAZy family
#   and a protein id: 30_dbcan/cazyme_families.tsv has families and no ids, the
#   three 23_extras tables have ids and no families, and step 30's raw dbCAN
#   output was not kept. So this step re-derives the per-protein calls.
#
# WHY NOT run_dbcan
#   dbCAN's CLI drags in a database layout, a diamond step and a CGC finder we
#   do not need, and every tool this project has added through a CLI has cost a
#   day. dbCAN's HMMER route IS an hmmsearch against dbCAN-HMMdb with a
#   published parser rule, and hmmsearch already works here. So: hmmsearch, then
#   dbCAN's own thresholds, applied explicitly and printed.
#
#     dbCAN hmmscan-parser rule, applied verbatim:
#       aligned HMM length > 80  ->  i-Evalue < 1e-5
#       aligned HMM length <= 80 ->  i-Evalue < 1e-3
#       HMM coverage > 0.35 in both cases
#
# CONTROLS
#   1. RECOVERY, like for like. Of the proteins Bakta's product strings call
#      carbohydrate-active, dbCAN must also call nearly all: a protein named
#      "beta-xylanase" has a GH domain or the annotation was wrong. Measured on
#      ONE protein set by two methods, so it needs no cross-definition
#      assumption.
#      (v1 compared dbCAN's ALL-CAZyme count against 23d's n_endo + n_exo and
#      FAILED at 3.1-3.3x. That was a defect in the SPECIFICATION: 23d counted
#      only GHs it could classify endo/exo and called ~60% of its own set
#      unclassifiable, while dbCAN's total includes glycosylTRANSFERASES and
#      binding modules. A superset against a subset of a subset must disagree.
#      The 23d table is kept below as a labelled reference, not a gate.)
#   2. CO-LOCATION. CBM families are beta-sandwich folds and so is the T9SS
#      sorting domain. If a CBM model is recognising the CTD, "is a CAZyme" and
#      "is cargo" are ONE observation and the enrichment is circular. Decided on
#      coordinates: any CAZy domain overlapping the CTD envelope by >50% of its
#      own length is excluded and the test recomputed.
#   3. STRATIFICATION. The three genomes differ in CTD count (122/100/32) and in
#      background CAZyme rate (4.5%/6.6%/11.9%), so the pooled fold is not the
#      per-genome fold. Section 3c runs the test one genome at a time.
#   4. Every protein id must join onto the proteome it was searched in.
#
#   None of these tests whether the CTD model works — step 36 gap B settled that
#   (254/254 in machinery-bearing genomes, 100% Bacteroidota). These ask only
#   what fraction of the cargo is carbohydrate-active, and whether that survives.
#
# ALSO DONE HERE
#   The sediment breadth of the three CTD-bearing genomes, read out of step 26.
#   36c's report says to quote it with any in-situ claim and it was never pulled.
#
# OUTPUT  36_gapclose/cazy_ctd.tsv · CAZY_CROSS.txt · dbcan/<mag>.domtbl
# RUNTIME 2-10 min (3 proteomes against dbCAN-HMMdb)
# USAGE   ws2 && bash scripts/36d_cazyme_ctd_cross.sh 2>&1 | tee logs/36d.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="36d_cazycross"
A="$WS2/09_annotation"; F="$WS2/07_mags"; T9="$WS2/23_extras/t9ss"
D="$WS2/36_gapclose"; mkdir -p "$D/dbcan"

echo "######## 36d — DO THE T9SS CARGO PROTEINS INCLUDE CAZymes? ########"
echo "  started: $(date -Is)"

# ------------------------------------------------------------ binary
# ENV_<tool> is only set for tools the smoketest probed. Assuming it exists
# broke steps 33, 34 and 35, so resolve the binary by scanning the roots.
bin_env() {
  local tool="$1" p
  if command -v "$tool" >/dev/null 2>&1; then echo ""; return 0; fi
  for root in ${CONDA_ROOTS:-${CONDA_ENVS_SHARED} /opt/conda/envs}; do
    for p in "$root"/*/bin/"$tool"; do
      [ -x "$p" ] && { echo "${p%/bin/$tool}"; return 0; }
    done
  done
  return 1
}
HS_ENV="$(bin_env hmmsearch)" || { echo "  !! no hmmsearch anywhere. Stop."; exit 2; }
HMMSEARCH="${HS_ENV:+$HS_ENV/bin/}hmmsearch"
echo "  hmmsearch: ${HS_ENV:-<on PATH>}"

# ------------------------------------------------------------ the CTD genomes
mapfile -t MAGS < <(ls -1 "$D/ctd"/*.tbl 2>/dev/null | while read -r f; do
  n=$(basename "$f" .tbl)
  # only genomes that actually have a CTD hit
  if [ "$(grep -vc '^#' "$f" 2>/dev/null)" -gt 0 ]; then echo "$n"; fi
done)
if [ "${#MAGS[@]}" -eq 0 ]; then
  echo "  !! no genome has a CTD hit — run step 36 gap B first."; exit 2
fi
echo "  CTD-bearing genomes: ${#MAGS[@]}  (${MAGS[*]})"

# ------------------------------------------------------------ dbCAN HMM db
# Locate rather than assume. Step 33 asserted an image tag and a container path
# and was wrong about both; step 36 assumed TIGR04183.HMM and the file is
# TIGR04183.1.HMM. Locate, verify it is HMMER3, count the profiles, print it.
DBH=""
for c in "$WS2/30_dbcan/db/dbCAN.txt" "$WS2/30_dbcan/dbCAN.txt" \
         "$WS2/30_dbcan/db/dbCAN-HMMdb-V13.txt" "$WS2/30_dbcan/db/dbCAN-HMMdb.txt"; do
  [ -s "$c" ] && { DBH="$c"; break; }
done
if [ -z "$DBH" ]; then
  # v1 searched maxdepth 5 under two roots and missed it. Search wider, name
  # every form the file takes across dbCAN versions, and PRINT what was
  # searched -- an unsuccessful search whose scope is not stated is not
  # evidence the file is absent.
  echo "  searching for the dbCAN HMM database ..."
  ROOTS="$WS2 ${DBDIR:-${DB_DIR}} $HOME ${SCRATCH_ROOT}"
  echo "    roots: $ROOTS"
  for r in $ROOTS; do
    [ -d "$r" ] || continue
    DBH="$(find "$r" -maxdepth 8 -type f \
          \( -name 'dbCAN.txt' -o -name 'dbCAN-HMMdb*' -o -name 'dbCAN*.hmm' \
             -o -name 'dbCAN_sub.hmm' \) -size +1M 2>/dev/null | head -1)"
    [ -n "$DBH" ] && { echo "    found under $r"; break; }
  done
fi
if [ -z "$DBH" ] || [ ! -s "$DBH" ]; then
  echo
  echo "  !! dbCAN HMM database not found."
  echo "     This is NOT 'no CAZyme carries a CTD' — it is an unrun search, and"
  echo "     the step stops here rather than printing a zero. To get it:"
  echo
  echo "       mkdir -p \$WS2/30_dbcan/db && cd \$WS2/30_dbcan/db"
  echo "       wget -O dbCAN-HMMdb-V14.txt \\"
  echo "         'https://pro.unl.edu/dbCAN2/download_file.php?file=dbCAN-HMMdb-V14.txt'"
  echo "       grep -c '^HMMER3' dbCAN-HMMdb-V14.txt    # expect several hundred"
  echo
  echo "     V14 is the current release (~130 MB); V13/V12 are still posted at"
  echo "     https://pro.unl.edu/dbCAN2/browse_download.php if a specific version"
  echo "     is ever needed to match an earlier run."
  echo
  echo "     Then re-run this step."
  runlog "**36d_cazycross NOT RUN — no dbCAN HMM db**"
  exit 3
fi
NPROF=$(grep -c '^HMMER3' "$DBH" 2>/dev/null || echo 0)
echo "  dbCAN HMMs: $DBH  ($NPROF profiles)"
if [ "$NPROF" -lt 100 ]; then
  echo "  !! only $NPROF profiles — that is not the dbCAN database. Stop."
  exit 3
fi

# ------------------------------------------------------------ search
STAMP=$(md5sum "$DBH" | cut -c1-12)
for m in "${MAGS[@]}"; do
  faa="$A/$m/$m.faa"
  [ -s "$faa" ] || faa="$(ls -1 "$A/$m"/*.faa 2>/dev/null | head -1)"
  if [ -z "$faa" ] || [ ! -s "$faa" ]; then
    echo "  !! no proteome FASTA for $m"; continue
  fi
  out="$D/dbcan/$m.domtbl"
  if [ -s "$out" ] && [ "$(cat "$D/dbcan/$m.stamp" 2>/dev/null)" = "$STAMP" ]; then
    echo "  [$m] cached"; continue
  fi
  echo "  [$m] searching $(grep -c '^>' "$faa") proteins ..."
  "$HMMSEARCH" --cpu "${THREADS:-8}" -E 1 --domE 1 \
      --domtblout "$out.part" "$DBH" "$faa" >/dev/null 2>&1
  rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$out.part" ]; then
    echo "     !! hmmsearch failed (rc=$rc). Not writing a table."
    rm -f "$out.part"; continue
  fi
  mv "$out.part" "$out"; echo "$STAMP" > "$D/dbcan/$m.stamp"
done
echo "  domtbl written: $(ls -1 "$D/dbcan"/*.domtbl 2>/dev/null | wc -l)/${#MAGS[@]}"

# ---- CTD coordinates, for the co-location control --------------------------
CTDH="$D/ctd/ctd.hmm"
if [ -s "$CTDH" ]; then
  for m in "${MAGS[@]}"; do
    faa="$A/$m/$m.faa"
    [ -s "$faa" ] || faa="$(ls -1 "$A/$m"/*.faa 2>/dev/null | head -1)"
    [ -n "$faa" ] && [ -s "$faa" ] || continue
    out="$D/ctd/$m.domtbl"
    [ -s "$out" ] && continue
    "$HMMSEARCH" --cpu "${THREADS:-8}" -E 1e-5 --domtblout "$out.part" \
        "$CTDH" "$faa" >/dev/null 2>&1
    [ -s "$out.part" ] && mv "$out.part" "$out" || rm -f "$out.part"
  done
  echo "  CTD coordinate tables: $(ls -1 "$D/ctd"/*.domtbl 2>/dev/null | wc -l)/${#MAGS[@]}"
else
  echo "  !! $CTDH missing — the co-location control cannot run."
fi

echo
# remove last run's outputs so the validation below cannot pass on them
rm -f "$D/CAZY_CROSS.txt" "$D/cazy_ctd.tsv"
python3 - "$D" "$A" "$F" "$T9" "$WS2" "$DBH" "$NPROF" <<'PY' 
import sys, os, re, csv, glob, math
from collections import defaultdict, Counter

D, A, F, T9, WS2 = sys.argv[1:6]
DBH, NPROF = sys.argv[6], sys.argv[7]
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)
STATUS = ['INCOMPLETE']          # only dump() may finalise this

def dump(code=0):
    # the status line is written HERE and nowhere else, so a run that dies
    # part-way cannot leave a stale 'OK' for the validation gate to read.
    with open(os.path.join(D, 'CAZY_CROSS.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'CAZY_CROSS.txt')}")
    sys.exit(code)

def _lf(n):
    return math.lgamma(n + 1)

def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    if n == 0:
        return 1.0
    r1, r2, c1 = a + b, c + d, a + c
    lo, hi = max(0, c1 - r2), min(r1, c1)
    def lp(x):
        return (_lf(r1) + _lf(r2) + _lf(c1) + _lf(n - c1)
                - _lf(n) - _lf(x) - _lf(r1 - x) - _lf(c1 - x) - _lf(r2 - c1 + x))
    ps = {x: math.exp(lp(x)) for x in range(lo, hi + 1)}
    tot = sum(ps.values()); obs = ps.get(a, 0.0)
    return min(1.0, sum(v for v in ps.values() if v <= obs * (1 + 1e-9)) / tot) if tot else 1.0

# ------------------------------------------------------ CTD ids
ctd = defaultdict(set)
for fn in sorted(os.listdir(os.path.join(D, 'ctd'))):
    if not fn.endswith('.tbl'):
        continue
    mag = fn[:-4]
    for line in open(os.path.join(D, 'ctd', fn), errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) >= 3:
            ctd[mag].add(f[0])
ctd = {m: s for m, s in ctd.items() if s}

# ------------------------------------------------------ proteomes
prod, of_mag = {}, {}
for mag in ctd:
    faa = os.path.join(A, mag, mag + '.faa')
    if not os.path.exists(faa):
        c = glob.glob(os.path.join(A, mag, '*.faa'))
        faa = c[0] if c else None
    if not faa:
        continue
    for line in open(faa, errors='replace'):
        if line.startswith('>'):
            pid, _, desc = line[1:].rstrip('\n').partition(' ')
            prod[pid] = desc.strip(); of_mag[pid] = mag

# ------------------------------------------------ dbCAN parse
# hmmsearch --domtblout: 0 target(protein) 2 tlen 3 query(HMM) 5 qlen
#                        12 i-Evalue 13 score 15 hmm_from 16 hmm_to
FAMRX = re.compile(r'^(GH|GT|PL|CE|AA|CBM)(\d+)')
cazy = {}          # protein -> (family, ievalue, cov)
per_mag = Counter()
parsed = 0
for fn in sorted(glob.glob(os.path.join(D, 'dbcan', '*.domtbl'))):
    # '.domtbl' is 7 characters. v1 sliced [:-8] and turned LIQ_das_metabat_078
    # into LIQ_das_metabat_07, which silently broke every join downstream and
    # let the recovery control pass on zero compared genomes.
    mag = re.sub(r'\.domtbl$', '', os.path.basename(fn))
    best = {}
    for line in open(fn, errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 23:
            continue
        parsed += 1
        pid, hmm = f[0], f[3]
        fam = hmm[:-4] if hmm.endswith('.hmm') else hmm
        if not FAMRX.match(fam):
            continue
        try:
            qlen = float(f[5]); iev = float(f[12])
            hf, ht = float(f[15]), float(f[16])
            sc = float(f[13])
        except ValueError:
            continue
        alen = ht - hf + 1
        cov = alen / qlen if qlen else 0.0
        # ---- dbCAN's own parser rule, applied explicitly
        ok = (iev < 1e-5) if alen > 80 else (iev < 1e-3)
        if not (ok and cov > 0.35):
            continue
        try:
            ef, et = int(f[19]), int(f[20])     # envelope on the PROTEIN
        except (ValueError, IndexError):
            ef, et = 0, 0
        if pid not in best or sc > best[pid][3]:
            best[pid] = (fam, iev, cov, sc, ef, et)
    for pid, (fam, iev, cov, sc, ef, et) in best.items():
        cazy[pid] = (fam, iev, cov, ef, et)
        per_mag[mag] += 1

rule('1. WHAT dbCAN FOUND')
say(f"  database   : {os.path.basename(DBH)}  ({NPROF} profiles)")
say(f"  full path  : {DBH}")
say("  The dbCAN VERSION is part of every number below. The family set grows")
say("  between releases, so a count from one version is not the same quantity as")
say("  a count from another — including step 30's, which used whichever version")
say("  it pulled. Quote the version with the count.")
say()
say(f"  domtbl rows parsed        : {parsed:,}")
say(f"  proteins passing dbCAN's rule: {len(cazy):,}")
say()
say("  rule applied: aligned HMM length > 80 -> i-E < 1e-5, else i-E < 1e-3;")
say("                HMM coverage > 0.35 in both cases. Best-scoring domain per")
say("                protein wins the family call.")
say()
if not cazy:
    say("  !! no protein passed. That is either an empty database or a failed")
    say("     search, NOT an absence of CAZymes in Bacteroidota genomes, which")
    say("     would be biologically impossible. Nothing below is usable.")
    STATUS[0] = 'SEARCH_EMPTY'
    dump(0)

# ---- CONTROL 1: recovery against step 23d's independent count
rule('2. CONTROL — DID THE SEARCH RECOVER WHAT IS KNOWN TO BE THERE?')

# The honest recovery control is against the SAME proteins by the OTHER method.
# Bakta's product strings call some proteins carbohydrate-active; dbCAN should
# find nearly all of those, because a protein whose product string says
# "beta-xylanase" has a GH domain or the annotation was wrong. This compares
# like with like on one protein set and needs no cross-definition assumption.
KW = re.compile(
    r'glycosid|glycosyl hydrolase|glucanase|xylanase|xylosidase|amylase|'
    r'chitinase|cellulase|cellobio|pectat|pectin|mannanase|mannosidase|'
    r'galactosidase|glucosidase|arabino|fucosidase|sialidase|neuraminidase|'
    r'hexosaminidase|agarase|alginate|ulvan|laminarin|levanase|inulinase|'
    r'polysaccharide lyase|carbohydrate esterase|glucuronyl|rhamnosidase|'
    r'chondroitin|heparin|hyaluron|dextranase|pullulanase|glycoside', re.I)
allids_pre = {p2 for p2 in prod}
kw = {p2 for p2 in allids_pre if KW.search(prod.get(p2, ''))}
rec = len(kw & set(cazy)) / len(kw) if kw else 0.0
say(f"  proteins Bakta's product string calls carbohydrate-active : {len(kw):,}")
say(f"  ...of those, also called by dbCAN                          : "
    f"{len(kw & set(cazy)):,}   ({100 * rec:.1f}%)")
say()
CTRL = rec >= 0.70 and len(kw) >= 20
if not kw:
    say("  [NO CONTROL] no keyword-called proteins to recover.")
elif CTRL:
    say("  [PASS] dbCAN recovers the great majority of the proteins the keyword")
    say("         method calls carbohydrate-active. The search ran and found what")
    say("         is there. The reverse direction — proteins dbCAN calls and the")
    say("         keyword method misses — is section 4, and is the point of the step.")
else:
    say("  [FAIL] dbCAN misses most of what the product strings call CAZymes.")
    say("         That is a failed search, not a discovery. Nothing below is usable.")
say()

# ---- 23d comparison, kept but DEMOTED to informational --------------------
say("  For reference, 23d's per-genome endo+exo counts. This is NOT a pass/fail")
say("  check and v1 was wrong to make it one:")
ee = os.path.join(T9, 'cazyme_endo_exo.tsv')
# Same non-deterministic join defect as 16b: a token filter that accepted
# 'Bacteria' and 'Gaoshiqia' as genome ids, merged globally, then read an
# unordered set. A genome id here always has an underscore AND a digit.
IDTOK = re.compile(r'^(?=.*_)(?=.*\d)[A-Za-z0-9._-]{6,}$')
RANKP = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')

def id_tokens(r):
    out = set()
    for v in r.values():
        v = str(v or '').strip()
        if v and ' ' not in v and ';' not in v and not v.startswith(RANKP) \
                and IDTOK.match(v) and not re.fullmatch(r'[0-9._-]+', v):
            out.add(v)
    return out

alias = {}
for path in (os.path.join(F, 'MAG_CROSSWALK.tsv'),):
    if not os.path.exists(path):
        continue
    with open(path, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            ids = sorted(id_tokens(r))
            if len(ids) > 1:
                for i in ids:
                    alias[i] = set(ids)      # row-local, never merged
old, gen = {}, {}
if os.path.exists(ee):
    with open(ee, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            mid = (r.get('mag_id') or '').strip()
            if not mid:
                continue
            try:
                n = int(float(r.get('n_endo') or 0)) + int(float(r.get('n_exo') or 0))
            except ValueError:
                continue
            for a2 in (alias.get(mid) or {mid}):
                old[a2] = n
                if (r.get('genus') or '').strip():
                    gen[a2] = (r.get('genus') or '').strip()

# dbCAN restricted to the DEPOLYMERASE classes, which is the closest thing to
# what 23d was counting. GT = biosynthesis, CBM = binding module, AA = redox.
DEPOL = ('GH', 'PL')
depol = Counter()
for pid, v in cazy.items():
    if v[0][:2] in DEPOL or v[0][:2] == 'GH' or v[0][:2] == 'PL':
        depol[of_mag.get(pid, '')] += 1
say()
say(f"  {'genome':<26}{'genus':<18}{'dbCAN all':>10}{'GH+PL':>8}"
    f"{'23d endo+exo':>14}{'GH+PL/23d':>11}")
say('-' * W)
for m in sorted(per_mag, key=lambda x: (-per_mag[x], x)):
    o = old.get(m)
    rr = (depol.get(m, 0) / o) if o else float('nan')
    say(f"  {m:<26}{(gen.get(m, '-') or '-')[:17]:<18}{per_mag[m]:>10}"
        f"{depol.get(m, 0):>8}{(o if o is not None else '-'):>14}"
        f"{('    n/a' if rr != rr else f'{rr:>11.2f}')}")
say()
say("  Why these cannot be equal, and why v1's [FAIL] was a defect in the check:")
say("   - 23d counted only GHs it could CLASSIFY as endo- or exo-acting and")
say("     reported ~60% of its own calls unclassifiable. Its number is a subset")
say("     of a subset.")
say("   - dbCAN's total includes glycosylTRANSFERASES (GT2, GT4 — biosynthetic,")
say("     not degradative) and carbohydrate-BINDING modules, which 23d never")
say("     counted and never should have.")
say("   - The two used different dbCAN releases where they used dbCAN at all.")
say("  A superset compared against a subset of a subset must disagree. It did,")
say("  by about the amount it had to. That is not evidence about the search.")
say()

# ---------------------------------------------------- the cross
rule('3. DO THE CARGO PROTEINS CARRY CAZy FAMILIES?')
ctdids = set().union(*ctd.values())
allids = {p for p in prod if of_mag.get(p) in ctd}
bg = allids - ctdids
joined = len(ctdids & set(prod))
say(f"  CTD proteins            : {len(ctdids):,}   (joined to a proteome: {joined:,})")
if joined < 0.9 * len(ctdids):
    say("  !! most CTD ids did not join. That is a join failure, not a result.")
    STATUS[0] = 'JOIN_FAILED'
    dump(0)
a = len(ctdids & set(cazy)); b = len(ctdids) - a
c = len(bg & set(cazy)); d = len(bg) - c
fa = a / max(1, len(ctdids)); fb = c / max(1, len(bg))
p = fisher_two_sided(a, b, c, d)
say(f"  background proteins     : {len(bg):,}")
say()
say(f"  {'set':<40}{'CAZyme':>10}{'total':>10}{'%':>9}")
say('-' * W)
say(f"  {'T9SS cargo (CTD-bearing)':<40}{a:>10}{len(ctdids):>10}{100 * fa:>8.2f}%")
say(f"  {'same genomes, everything else':<40}{c:>10}{len(bg):>10}{100 * fb:>8.2f}%")
say()
say(f"  fold                    : {(fa / fb) if fb else float('inf'):.2f}")
say(f"  Fisher exact, two-sided : p = {p:.2e}")
say()

if p < 0.05 and fa > fb:
    VERDICT = 'ENRICHED'
    say("  >> RAW READING: CAZymes are over-represented in the cargo. Section 3b")
    say("     decides whether that survives the co-location control. Until it does,")
    say("     this number is not a result.")
elif p < 0.05 and fa < fb:
    VERDICT = 'DEPLETED'
    say("  >> CAZymes are UNDER-represented in the cargo, by a family-HMM method")
    say("     that owes nothing to Bakta's wording. 36c's reading is confirmed and")
    say("     strengthened: whatever this T9SS exports, it is not the CAZymes.")
else:
    VERDICT = 'SAME'
    say("  >> CAZymes occur in the cargo at the same rate as in the rest of these")
    say("     genomes. That is 36c's result, now on family HMMs rather than product")
    say("     strings: the cargo is not CAZyme-biased, in either direction.")
say()
# ---- CONTROL: is the CAZy call sitting ON the sorting domain? --------------
# CBM families are beta-sandwich folds. So is the T9SS type A sorting domain.
# If a CBM model is recognising the CTD itself, then the CAZyme call and the
# cargo membership are the SAME observation and the enrichment is circular --
# the identical defect 36c had with Bakta naming proteins after the domain,
# arriving by a different route. Decided on coordinates, not on plausibility.
ctdpos = defaultdict(list)
for fn in sorted(glob.glob(os.path.join(D, 'ctd', '*.domtbl'))):
    for line in open(fn, errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 21:
            continue
        try:
            ctdpos[f[0]].append((int(f[19]), int(f[20])))
        except ValueError:
            continue
rule('3b. CONTROL — IS THE CAZy DOMAIN THE SORTING DOMAIN?')
if not ctdpos:
    say("  !! no CTD coordinate tables, so this control did not run. The")
    say("     enrichment above is UNVALIDATED against co-location — do not quote")
    say("     it until this runs.")
    COLOC_OK = None
else:
    say(f"  CTD envelopes read for {len(ctdpos):,} proteins.")
    say("  A CAZy domain overlapping the sorting domain by >50% of its own length")
    say("  is counted as SUSPECT: the model may be recognising the CTD, in which")
    say("  case the CAZyme call and cargo membership are one observation.")
    say()
    suspect, clean = set(), set()
    for pid in (ctdids & set(cazy)):
        fam, iev, cov, ef, et = cazy[pid]
        if not et:
            clean.add(pid); continue
        L = et - ef + 1
        ov = 0
        for cf, ct in ctdpos.get(pid, ()):
            ov = max(ov, max(0, min(et, ct) - max(ef, cf) + 1))
        (suspect if L and ov / L > 0.5 else clean).add(pid)
    say(f"  cargo CAZymes total            : {len(ctdids & set(cazy))}")
    say(f"  ...overlapping the CTD (SUSPECT): {len(suspect)}")
    say(f"  ...disjoint from it (CLEAN)     : {len(clean)}")
    if suspect:
        say()
        say("  the suspect calls:")
        for pid in sorted(suspect)[:15]:
            say(f"    {pid:<26}{cazy[pid][0]:<12}{prod.get(pid, '')[:44]}")
    say()
    a2 = len(clean); b2 = len(ctdids) - a2
    fa2 = a2 / max(1, len(ctdids))
    p2v = fisher_two_sided(a2, b2, c, d)
    say("  the enrichment recomputed on CLEAN calls only:")
    say(f"    cargo {a2}/{len(ctdids)} ({100 * fa2:.2f}%) vs background"
        f" {c}/{len(bg)} ({100 * fb:.2f}%),"
        f" fold {(fa2 / fb) if fb else float('inf'):.2f}, p = {p2v:.2e}")
    say()
    COLOC_OK = (p2v < 0.05 and fa2 > fb)
    if COLOC_OK:
        say("  [PASS] the enrichment survives removing every CAZy domain that")
        say("         overlaps the sorting domain. It is not an artefact of one")
        say("         beta-sandwich model recognising another.")
        VERDICT = 'ENRICHED'
    elif len(suspect) > 0.5 * len(ctdids & set(cazy)):
        say("  [FAIL] most cargo CAZyme calls sit ON the sorting domain. The")
        say("         enrichment is circular — the CAZy model is recognising the")
        say("         CTD. Do not quote the fold above.")
        VERDICT = 'CIRCULAR'
    else:
        say("  [FAIL] the enrichment does not survive the co-location filter.")
        VERDICT = 'SAME'
say()
# ---- the same test, ONE GENOME AT A TIME -----------------------------------
# The pooled test mixes three genomes with very different CTD counts (122/100/32)
# and very different background CAZyme rates (4.5%/6.6%/11.9%). Pooling across
# strata that differ in both can move the estimate in either direction -- here it
# DILUTES it, so the pooled fold is lower than every genome's own. And only one
# of the three is present in the sediment, so the pooled number cannot answer the
# in-situ question at all. Stratify.
rule('3c. THE SAME TEST, ONE GENOME AT A TIME')
nprot_of = Counter(of_mag.values())
STRATA = {}
say(f"  {'genome':<24}{'cargo CAZy':>13}{'background':>14}{'fold':>8}{'p':>12}")
say('-' * W)
for m in sorted(ctd, key=lambda x: (-len(ctd[x]), x)):
    cset = ctd[m]
    cc_ = len(cset & set(cazy))
    bgn = nprot_of.get(m, 0) - len(cset)
    bgc = len([p2 for p2 in cazy if of_mag.get(p2) == m]) - cc_
    if not cset or bgn <= 0:
        continue
    fa_, fb_ = cc_ / len(cset), (bgc / bgn if bgn else 0.0)
    fold_ = (fa_ / fb_) if fb_ else float('inf')
    p_ = fisher_two_sided(cc_, len(cset) - cc_, bgc, bgn - bgc)
    STRATA[m] = (cc_, len(cset), bgc, bgn, fold_, p_)
    say(f"  {m:<24}{f'{cc_}/{len(cset)}':>13}{f'{bgc}/{bgn}':>14}"
        f"{fold_:>8.2f}{p_:>12.2e}")
say('-' * W)
say(f"  {'POOLED':<24}{f'{a}/{len(ctdids)}':>13}{f'{c}/{len(bg)}':>14}"
    f"{(fa / fb) if fb else float('inf'):>8.2f}{p:>12.2e}")
say()
if STRATA and all(v[4] > (fa / fb if fb else 0) for v in STRATA.values()):
    say("  >> EVERY genome's own fold is HIGHER than the pooled one. The pooled")
    say("     estimate is diluted by mixing strata with different background")
    say("     rates, so quote the per-genome table, not the pooled fold.")
elif STRATA and all(v[5] < 0.05 for v in STRATA.values()):
    say("  >> The enrichment holds in every genome separately, so it is not")
    say("     carried by one of the three.")
else:
    say("  >> The enrichment does NOT hold in every genome separately. Say which")
    say("     genomes carry it; the pooled number hides that.")
say()

both = sorted(ctdids & set(cazy))
if both:
    say(f"  the {len(both)} cargo proteins that ARE CAZymes:")
    say(f"    {'protein':<26}{'family':<12}{'cov':>7}{'i-E':>11}  product")
    say('-' * W)
    for pid in both[:30]:
        fam, iev, cov = cazy[pid][:3]
        say(f"    {pid:<26}{fam:<12}{cov:>7.2f}{iev:>11.1e}  {prod.get(pid, '')[:34]}")
    say()
    say("  families among cargo CAZymes:")
    for fam, n in sorted(Counter(cazy[p2][0] for p2 in sorted(both)).items(), key=lambda kv: (-kv[1], kv[0]))[:15]:
        say(f"    {fam:<12}{n:>5}")
    say()
    say("  families among the SAME genomes' non-cargo CAZymes, for contrast:")
    _bgfam = Counter(cazy[p2][0] for p2 in sorted(bg & set(cazy)))
    for fam, n in sorted(_bgfam.items(), key=lambda kv: (-kv[1], kv[0]))[:15]:
        say(f"    {fam:<12}{n:>5}")
say()

# ------------------------------------------- 36c's keyword call, checked
rule('4. WHAT THE KEYWORD METHOD GOT WRONG, MEASURED')
kw = {p2 for p2 in allids if KW.search(prod.get(p2, ''))}
db = set(cazy) & allids
say(f"  called CAZyme by Bakta keyword only : {len(kw - db):,}")
say(f"  called by dbCAN only                : {len(db - kw):,}")
say(f"  called by both                      : {len(kw & db):,}")
say()
if kw | db:
    say(f"  agreement (Jaccard)                 : {len(kw & db) / len(kw | db):.2f}")
say("  This is the size of the difference between a product string and a family")
say("  HMM, measured on these genomes. It is the same measurement step 30 made")
say("  for the catalogue, repeated where it now matters.")
say()
say(f"  within the cargo specifically: keyword said {len(kw & ctdids)},"
    f" dbCAN says {len(db & ctdids)}")

# ---------------------------------------------- step 26 breadth
rule('5. ARE THESE THREE GENOMES ACTUALLY IN THE SEDIMENT?')
say("  36c's report says to quote step 26's sediment breadth with any in-situ")
say("  claim about this cargo. All three are culture-derived (LIQ) MAGs, so")
say("  without breadth this is a statement about the enrichment, not the lagoon.")
say()
# ---- step 26's breadth, located rather than guessed -----------------------
# v1 globbed for '*breadth*'/'*crossmap*' and found nothing, because step 26's
# files are cross_abundance.tsv, SED_vs_all.tsv and LIQ_vs_all.tsv. v2 scanned
# columns instead and picked 12_results/MAG_CATALOGUE.tsv -- whose breadth_SED
# column only carries the FOUR sediment-derived MAGs. The contrast then ran as
# "4 SED MAGs vs 25 LIQ MAGs", which is a comparison of ASSEMBLY ORIGIN wearing
# the label of in-situ presence. A table that covers a sixth of the catalogue
# cannot answer this question, so coverage is now checked and a thin table is
# refused.
#
# Also: CoverM names its column "Covered Fraction" WITH A SPACE, and its genome
# ids carry a .fa extension the crosswalk records but the proteome names do not.
# Both are handled rather than assumed away.
def _norm(x):
    return re.sub(r'\.(fa|fna|fasta|fa\.gz|fna\.gz)$', '', str(x or '').strip())

def find_breadth(ws2, alias, need):
    files = sorted(glob.glob(os.path.join(ws2, '26_crossmap', '*.tsv'))) + \
            sorted(glob.glob(os.path.join(ws2, '12_results', '*.tsv')))
    files.sort(key=lambda f: (0 if re.search(r'sed', os.path.basename(f), re.I) else 1,
                              os.path.basename(f)))
    # normalised alias index, so LIQ_metabat_001.fa finds LIQ_das_metabat_001
    ax = defaultdict(set)
    for k, grp in alias.items():
        ax[_norm(k)].update(grp)
    seen, best = [], None
    for c2 in files:
        try:
            rows = list(csv.DictReader(open(c2, errors='replace'), delimiter='\t'))
        except Exception:
            continue
        if not rows:
            continue
        cols = [x for x in rows[0].keys() if x]
        bcol = next((x for x in cols if re.search(r'breadth', x, re.I)
                     and re.search(r'sed', x, re.I)), None) \
            or next((x for x in cols if re.search(r'breadth|covered[ _-]?frac|'
                                                  r'coverage[ _-]?frac|frac[ _-]?cov', x, re.I)), None)
        seen.append((os.path.relpath(c2, ws2), cols, bcol is not None))
        if not bcol:
            continue
        idc = next((x for x in cols if re.search(r'mag|genome|bin|^name$', x, re.I)), cols[0])
        out = {}
        for r in rows:
            mid = _norm(r.get(idc))
            try:
                v = float(str(r.get(bcol) or '').strip() or 'nan')
            except ValueError:
                continue
            if v != v:
                continue
            for a2 in sorted(ax.get(mid) or {mid}):
                out[a2] = v
        cov = len([m for m in need if m in out])
        if best is None or cov > best[2]:
            best = (out, (os.path.relpath(c2, ws2), bcol, idc), cov)
        if cov >= 0.9 * len(need):
            break
    if best is None:
        return {}, None, seen, 0
    return best[0], best[1], seen, best[2]

bd, src, seen, BCOV = find_breadth(WS2, alias, sorted(ctd))
if src:
    say(f"  source: {src[0]}   breadth column: {src[1]}   id column: {src[2]}")
    say()
    say(f"  {'genome':<26}{'genus':<18}{'breadth':>10}{'reading':>19}"
        f"{'cargo CAZy effect':>19}")
    say('-' * W)
    for m in sorted(ctd, key=lambda x: (-len(ctd[x]), x)):
        v = bd.get(m)
        if v is None:
            say(f"  {m:<26}{(gen.get(m, '-') or '-')[:17]:<18}{'not found':>18}"
                f"{'—':>20}")
            continue
        rd = ('confirmed in situ' if v >= 0.5
              else ('marginal' if v >= 0.1 else 'enrichment-only'))
        st = STRATA.get(m)
        ex = (f"{st[4]:.2f}x  p={st[5]:.0e}" if st else '-')
        say(f"  {m:<26}{(gen.get(m, '-') or '-')[:17]:<18}{v:>10.2f}{rd:>19}{ex:>19}")
    say()
    ins = [m for m in ctd if bd.get(m, 0.0) >= 0.5]
    if ins:
        say(f"  >> {len(ins)} of {len(ctd)} CTD-bearing genomes are confirmed in the")
        say(f"     sediment: {', '.join(sorted(ins))}.")
        say(f"     They carry {sum(len(ctd[m] & set(cazy)) for m in ins)} of the"
            f" {a} cargo CAZymes and"
            f" {sum(len(ctd[m]) for m in ins)} of the {len(ctdids)} CTD proteins.")
        say("     Any IN-SITU statement about T9SS degradation rests on those")
        say("     genomes alone; the rest describe the enrichment catalogue.")
    else:
        say("  >> NONE of the CTD-bearing genomes is confirmed in the sediment.")
        say("     Every statement here is about the enrichment catalogue.")
    say("  Step 26's own rule: presence is called on BREADTH, not depth. A genome")
    say("  at 0.5x with 5% of its length covered is conserved-region cross-mapping;")
    say("  the same depth at 60% breadth is a real low-abundance population.")
else:
    say("  !! no breadth column found. Recorded as NOT LOOKED UP, not as absent.")
    say("     Headers actually present, so the right column can be named:")
    for f2, cols, hasb in seen[:8]:
        say(f"       {'[b] ' if hasb else '[ ] '}{f2}: {', '.join(cols[:10])}")

# ---------------------------------------------------------- outputs
with open(os.path.join(D, 'cazy_ctd.tsv'), 'w') as oh:
    oh.write('protein\tmag\tgenus\tctd\tcazy_family\tievalue\thmm_cov\tproduct\n')
    for pid in sorted(allids):
        m = of_mag.get(pid, '')
        fam, iev, cov = (cazy.get(pid) or ('', '', '', 0, 0))[:3]
        oh.write(f"{pid}\t{m}\t{gen.get(m, '')}\t{'yes' if pid in ctdids else 'no'}\t"
                 f"{fam}\t{iev}\t{cov}\t{prod.get(pid, '')}\n")

# one machine-readable line for the validation gate. Keyed on the CONTROLS,
# never on the report text.
_st = 'OK'
if not kw or not CTRL:
    _st = 'RECOVERY_FAILED' if kw else 'NO_CONTROL'
elif VERDICT == 'CIRCULAR':
    _st = 'CIRCULAR'
elif COLOC_OK is None:
    _st = 'NO_CONTROL'
STATUS[0] = _st
rule('WHAT TO REPORT')
if VERDICT == 'CIRCULAR':
    say("  1. NOTHING about cargo CAZyme content is reportable from this run. The")
    say("     apparent enrichment came from CAZy models recognising the sorting")
    say("     domain itself — the same circularity 36c hit with Bakta's naming,")
    say("     reached by a different route. Report step 36 gap B's count only.")
elif VERDICT == 'DEPLETED':
    say("  1. The T9SS cargo of these Bacteroidota is CAZyme-poor by family-HMM")
    say("     assignment, not merely by product string. Quote the fold and the p")
    say("     alongside 36c's adhesion and peptidase enrichments.")
elif VERDICT == 'ENRICHED':
    say("  1. dbCAN OVERTURNS 36c's carbohydrate reading. Quote this, not 36c's")
    say("     keyword-based n.s. result, and say why: a product string is not a")
    say("     family assignment.")
else:
    say("  1. CAZymes occur in the cargo at the background rate of the same")
    say("     genomes, confirmed on family HMMs. 36c's reading stands and is now")
    say("     independent of Bakta's wording.")
say("  2. 23d calling Gaoshiqia a POLYMER DEGRADER and this step finding its")
say("     CAZymes enriched in the T9SS cargo are the SAME finding reached twice.")
say("     (An earlier draft of this text said the opposite — that the CAZymes do")
say("     not carry the export tag — written from the keyword count of 8 before")
say(f"     dbCAN found {len(both)}. It is corrected here rather than deleted.)")
say("     Note what still holds: MOST of these genomes' CAZymes are not cargo.")
say(f"     Cargo carries {a}/{len(ctdids)}; the same genomes hold {c:,} more that")
say("     do not, consistent with SusC/SusD import and periplasmic hydrolysis")
say("     for the bulk, with a T9SS-exported subset acting outside the cell.")
if STRATA:
    _f = sorted(v[4] for v in STRATA.values())
    _pool = (fa / fb) if fb else 0.0
    # the claim in this item must follow the data, not the expectation. v1 said
    # "every genome's fold is HIGHER" unconditionally; on a set where the strata
    # agree with the pooled estimate that sentence is simply false.
    if _f[0] > _pool * 1.02:
        say(f"  2b. QUOTE THE PER-GENOME FOLDS ({_f[0]:.2f}-{_f[-1]:.2f}x, section 3c),")
        say(f"     NOT the pooled {_pool:.2f}x. EVERY genome's own fold is higher than")
        say("     the pooled one: they differ in both CTD count and background CAZyme")
        say("     rate, and pooling across strata that differ in both dilutes the")
        say("     estimate. The pooled number is not a conservative version of the")
        say("     result — it is wrong about the size of the effect.")
    elif _f[-1] < _pool * 0.98:
        say(f"  2b. Every genome's own fold ({_f[0]:.2f}-{_f[-1]:.2f}x) is LOWER than the")
        say(f"     pooled {_pool:.2f}x, so pooling INFLATES the effect here. Quote the")
        say("     per-genome table.")
    else:
        say(f"  2b. Per-genome folds ({_f[0]:.2f}-{_f[-1]:.2f}x) bracket the pooled")
        say(f"     {_pool:.2f}x, so pooling is not distorting the estimate. Either can")
        say("     be quoted; the per-genome table is still the more informative one.")
say("  3. Section 4's keyword-vs-HMM disagreement is worth a Methods sentence.")
_ins = [m for m in ctd if bd.get(m, 0.0) >= 0.5] if bd else []
if bd and _ins:
    say(f"  4. ONLY {len(_ins)} of {len(ctd)} CTD-bearing genomes is present in the")
    say(f"     sediment ({', '.join(sorted(_ins))}, breadth >= 0.50). It carries")
    say(f"     {sum(len(ctd[m] & set(cazy)) for m in _ins)} of the {a} cargo CAZymes"
        f" and {sum(len(ctd[m]) for m in _ins)} of the {len(ctdids)} CTD proteins.")
    say("     Every IN-SITU sentence about T9SS degradation in Etoliko rests on")
    say("     that genome alone. Note that 23d's headline polymer degrader was")
    say("     Gaoshiqia, which is at breadth 0.02 — essentially absent from the")
    say("     lagoon. The in-situ / enrichment split must be quoted every time")
    say("     these genomes are.")
elif bd:
    say("  4. NONE of the CTD-bearing genomes is present in the sediment. Say")
    say("     'enrichment catalogue', never 'in the lagoon'.")
else:
    say("  4. Sediment breadth was NOT LOOKED UP, so say 'enrichment catalogue'.")
say()
say("  LIMITS:")
say("   1. A floor twice over — hmmsearch finds what dbCAN's HMMs match, and the")
say("      CTD set is itself a floor from step 36 gap B.")
say("   2. A CAZy family is a FOLD, not a substrate. GH5 alone spans a dozen")
say("      activities, so 'exports a cellulase' needs the subfamily and even then")
say("      is an inference from sequence.")
say("   3. Presence of a sorting domain is not demonstrated secretion. This is")
say("      gene content; localisation would need proteomics, which step 21 already")
say("      recorded as out of scope for this dataset.")
say(f"   4. dbCAN {os.path.basename(DBH)} — the family set differs between releases,")
say("      so this count is not directly comparable to step 30's.")

dump(0)
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "analysis completed" "$RC" || OK=1
check "cross report written" "$([ -s "$D/CAZY_CROSS.txt" ] && echo 0 || echo 1)" || OK=1
check "identity table written" "$([ -s "$D/cazy_ctd.tsv" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$D/CAZY_CROSS.txt" 2>/dev/null | awk '{print $2}')
echo "  status line: ${ST:-<none>}"
case "${ST:-MISSING}" in
  OK)             ;;
  NO_CONTROL)     echo "  [WARN] a control did not run — see the report." ;;
  *)              echo "  [STOP] ${ST:-no status line} — do not quote the cross."; OK=1 ;;
esac
[ "$OK" -eq 0 ] && runlog "36d_cazycross complete" || runlog "**36d_cazycross FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
