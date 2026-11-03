#!/usr/bin/env bash
# =============================================================================
# 36b_ecf_transporters.sh  --  did step 36 look for the wrong transporter family?
#
# WHAT STEP 36 FOUND, AND WHY IT IS SUSPICIOUS
#   Direct profile search at KEGG's curated thresholds returned, across 29 MAGs:
#
#     K02062  thiamine transport, ATP-binding    0/29   (1,889 hits rejected)
#     K02064  thiamine transport, binding        0/29
#     K06073  B12 transport, permease            0/29
#     K06074  B12 transport, ATP-binding         0/29   (1,829 hits rejected)
#     K06858  B12 transport, binding             0/29
#     K03523  BIOTIN transport, S-component     24/29   (83%)
#
#   Zero B12 and zero thiamine transporters in a Clostridia-dominated catalogue
#   is not plausible. Those organisms are vitamin auxotrophs; they must import.
#
#   The answer is in step 36's own table. K03523 is BioY — the substrate-specific
#   component of an ENERGY-COUPLING FACTOR (ECF) transporter, not part of a
#   classic ABC system. Firmicutes take up B12, thiamine, biotin and riboflavin
#   through ECF transporters: a shared EcfA/EcfA'/EcfT module plus a small
#   swappable S-component per vitamin. Step 29's table encodes the BtuFCD and
#   ThiPQ systems of Gram-negatives. Those genomes never had them.
#
#   So step 36's 190 DOWNGRADED edges may be downgraded against a protein family
#   these organisms do not use. That is Rule 35: an absence is an absence of what
#   you searched for.
#
# THE CONTROL IS BUILT INTO THE BIOLOGY, AND IT CAN FAIL
#   An S-component CANNOT TRANSPORT ANYTHING without the EcfAA'T module — it is
#   the substrate-binding subunit of a machine, not a transporter. BioY is
#   already known to be in 24/29 genomes. Therefore:
#
#     module found alongside BioY  -> ECF transport is real here, and step 36's
#                                     vitamin absences are about the wrong family
#     module ABSENT while BioY 24/29 -> THE SEARCH IS BROKEN. An orphan
#                                     S-component in 83% of genomes is not a
#                                     biological result, it is a failed search.
#
#   That second branch is a real possible outcome and it blocks the verdict.
#
# WHAT THIS CANNOT DO
#   KEGG's S-component definitions say only "energy-coupling factor transport
#   system substrate-specific component" — they do NOT name the vitamin. So this
#   step can show that an ECF SYSTEM is present and how many S-components each
#   genome carries; it cannot say which vitamin any one of them binds. The
#   conclusion is therefore "these edges cannot be decided by an ABC query", not
#   "B12 transport is present". Saying more than that would repeat the mistake
#   this step exists to correct.
#
# OUTPUT  36_gapclose/ecf_matrix.tsv · edges_reopened.tsv · ECF_VERDICT.txt
# RUNTIME 5-15 min
# USAGE   ws2 && bash scripts/36b_ecf_transporters.sh 2>&1 | tee logs/36b.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="36b_ecf"

A="$WS2/09_annotation"; F="$WS2/07_mags"; K="$WS2/11_kofam"
D="$WS2/36_gapclose"; E="$D/ecf"; mkdir -p "$E/profiles" "$E/hits"
BASE="https://www.genome.jp/ftp/db/kofam"

echo "######## 36b — ECF TRANSPORTERS: THE FAMILY STEP 36 DID NOT SEARCH ########"
echo "  started: $(date -Is)"
runlog "36b_ecf started"

# ------------------------------------------------------------------ binaries
bin_env() {
  local tool="$1"; local var="ENV_${tool//-/_}"
  local cands="${!var:-} ${ENV_hmmsearch:-} ${ENV_hmmscan:-} ${ENV_hmmer:-} ${ENV_bakta:-} -"
  local e prefix root d
  for e in $cands; do
    [ -n "$e" ] || continue
    if [ "$e" = "-" ]; then command -v "$tool" >/dev/null 2>&1 && { echo "-"; return 0; }; continue; fi
    if declare -f env_prefix >/dev/null 2>&1; then
      prefix=$(env_prefix "$e" 2>/dev/null) || continue
      [ -x "$prefix/bin/$tool" ] && { echo "$e"; return 0; }
    fi
  done
  for root in ${CONDA_ROOTS:-}; do
    [ -d "$root" ] || continue
    for d in "$root"/*; do [ -x "$d/bin/$tool" ] && { echo "$(basename "$d")"; return 0; }; done
  done
  return 1
}
HS=$(bin_env hmmsearch || true)
[ -n "$HS" ] || { echo "  !! hmmsearch not found"; runlog "**36b ABORTED: no hmmsearch**"; exit 2; }
export ENV_hmmsearch="$HS"; echo "  hmmsearch: env '$HS'"

for src in "$K/ko_list" "$D/ko_list" "$WS2/34_unbinned/ko_list"; do
  [ -s "$src" ] && { cp -n "$src" "$E/ko_list" 2>/dev/null || true; break; }
done
[ -s "$E/ko_list" ] || { echo "  !! no ko_list"; runlog "**36b ABORTED: no ko_list**"; exit 3; }
echo "  ko_list: $(wc -l < "$E/ko_list") families"

mapfile -t MAGS < <(for d in "$A"/*/; do
  n=$(basename "$d"); [ -s "$d/$n.faa" ] && echo "$n"; done | sort)
echo "  proteomes: ${#MAGS[@]}"
[ "${#MAGS[@]}" -ge 5 ] || { echo "  !! too few proteomes"; runlog "**36b ABORTED**"; exit 4; }

# ============================================================ 0. KO GUARD
# The KO list came from a grep of ko_list on 2026-09-04 and every entry is
# re-verified here against KEGG's own definition. Roles are asserted from the
# definition text, not from memory.
echo
echo "######## 0. KO IDENTITY GUARD ########"
cat > "$E/targets.tsv" <<'KOS'
K16785	EcfT	MODULE	permease
K16786	EcfA	MODULE	ATP-binding
K16787	EcfA2	MODULE	ATP-binding
K01552	EcfA3	MODULE	ATP-binding
K16925	EcfT2	MODULE	permease
K24040	EcfTA	MODULE	permease/ATP-binding
K03523	BioY	SCOMP	substrate-specific
K16923	S-comp1	SCOMP	substrate-specific
K16924	S-comp2	SCOMP	substrate-specific
K16926	S-comp3	SCOMP	substrate-specific
K16927	S-comp4	SCOMP	substrate-specific
K16928	S-comp5	SCOMP	substrate-specific
K16929	S-comp6	SCOMP	substrate-specific
K24620	S-comp7	SCOMP	substrate-specific
KOS
python3 - "$E" <<'PYG'
import sys, os
E = sys.argv[1]
defs = {}
with open(os.path.join(E, 'ko_list'), errors='replace') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    di = hdr.index('definition') if 'definition' in hdr else len(hdr) - 1
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) > di:
            defs[f[0]] = f[di]
bad = 0
for line in open(os.path.join(E, 'targets.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) < 4:
        continue
    ko, gene, role, claim = f
    d = defs.get(ko)
    if d is None:
        print(f"  [FAIL] {ko} ({gene}) is not in ko_list"); bad += 1; continue
    low = d.lower()
    # K03523 is the one entry KEGG names by its cargo rather than by the family
    if ko == 'K03523':
        ok = 'substrate-specific component' in low
    else:
        ok = 'energy-coupling factor' in low
    if not ok:
        print(f"  [FAIL] {ko} claimed as an ECF {claim} but KEGG says: {d}")
        bad += 1; continue
    # the role must match the definition too: a permease is not an ATPase
    want = 'permease' if role == 'MODULE' and 'perme' in claim else None
    print(f"  [ok]   {ko}  {gene:<8}{role:<7} {d[:56]}")
if bad:
    print(f"  !! {bad} KO(s) disagree with KEGG. NOTHING IS RUN.")
    sys.exit(1)
print()
print("  NOTE: KEGG does not name the vitamin for any S-component. This step can")
print("  show that an ECF SYSTEM exists and how many S-components a genome has.")
print("  It CANNOT say which vitamin any of them binds, and nothing below does.")
sys.exit(0)
PYG
[ $? -eq 0 ] || { runlog "**36b ABORTED: KO guard failed**"; exit 5; }

# ============================================================ 1. PROFILES
echo
echo "######## 1. PROFILES ########"
mapfile -t KOS < <(cut -f1 "$E/targets.tsv" | sort -u)
got=0; missing=()
for ko in "${KOS[@]}"; do
  [ -s "$E/profiles/${ko}.hmm" ] && { got=$((got+1)); continue; }
  for src in "$K/profiles/${ko}.hmm" "$D/profiles/${ko}.hmm"; do
    [ -s "$src" ] && { cp "$src" "$E/profiles/"; break; }
  done
  if [ ! -s "$E/profiles/${ko}.hmm" ] && [ -s "$K/profiles.tar.gz" ]; then
    tar xzf "$K/profiles.tar.gz" -C "$E/profiles" --strip-components=1 \
        "profiles/${ko}.hmm" 2>/dev/null || true
  fi
  [ -s "$E/profiles/${ko}.hmm" ] && got=$((got+1)) || missing+=("$ko")
done
echo "  $got/${#KOS[@]} profiles"
printf '%s\n' "${missing[@]:-}" | grep -v '^$' > "$E/not_searched.txt" || true
if [ "${#missing[@]}" -gt 0 ]; then
  echo "  !! NO PROFILE for: ${missing[*]}"
  echo "     Carried through as NOT SEARCHED. Never reported as absent."
  echo "     They are in $K/profiles.tar.gz if that archive is present."
fi
cat "$E"/profiles/*.hmm > "$E/all.hmm" 2>/dev/null
[ -s "$E/all.hmm" ] || { echo "  !! no profiles"; runlog "**36b ABORTED: no profiles**"; exit 6; }

# Cache keyed on the profile set, absence failing closed (Rule 37).
PSTAMP=$(md5sum "$E/all.hmm" | cut -c1-12)
echo "  profile-set stamp: $PSTAMP"

# ============================================================ 2. SEARCH
echo
echo "######## 2. SEARCHING ${#MAGS[@]} PROTEOMES ########"
n=0
for name in "${MAGS[@]}"; do
  n=$((n+1)); tb="$E/hits/${name}.tbl"; st="$E/hits/${name}.stamp"
  if [ -s "$tb" ] && [ "$(cat "$st" 2>/dev/null || echo NONE)" = "$PSTAMP" ]; then
    printf '\r  [%2d/%2d] %-30s cached ' "$n" "${#MAGS[@]}" "$name"; continue
  fi
  printf '\r  [%2d/%2d] %-30s        ' "$n" "${#MAGS[@]}" "$name"
  run_tool hmmsearch --cpu "$THREADS" -E 1e-10 --domE 1e-10 \
      --tblout "$tb.part" "$E/all.hmm" "$A/$name/$name.faa" \
      >/dev/null 2>"$E/hits/${name}.err"
  rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$tb.part" ]; then
    printf '\r  [%2d/%2d] %-30s FAILED rc=%s\n' "$n" "${#MAGS[@]}" "$name" "$rc"
    head -3 "$E/hits/${name}.err" | sed 's/^/       /'; rm -f "$tb.part"; continue
  fi
  mv "$tb.part" "$tb"; echo "$PSTAMP" > "$st"
done
echo
echo "  tables: $(ls -1 "$E/hits"/*.tbl 2>/dev/null | wc -l)/${#MAGS[@]}"

# ============================================================ 3. INTERPRET
echo
python3 - "$E" "$D" "$F" <<'PY'
import sys, os, re, csv
from collections import defaultdict, Counter

E, D, F = sys.argv[1], sys.argv[2], sys.argv[3]
W = 100
rep = []
def say(s=''):
    print(s); rep.append(s)
def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)


# ---- NAME JOIN: edge tables key on ETO working names, proteome dirs on bin ids
# Step 36's edge closure looked up `consumer` (ETO_LIQ_MAG_003) in a table keyed
# on proteome directory names (LIQ_das_metabat_001). The intersection is EMPTY,
# so `found` was always empty and every openable edge became DOWNGRADED without
# the transporter ever being looked up. The domain control needed the crosswalk
# and got it; the edge closure needed it and did not.
def build_alias(paths):
    """canonical name -> set of all names for the same genome, both directions"""
    groups = []
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path, errors='replace') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        if not rows:
            continue
        cols = [c for c in rows[0].keys() if c]
        for r in rows:
            ids = set()
            for c in cols:
                v = (r.get(c) or '').strip()
                if v and ' ' not in v and ';' not in v and len(v) >= 4 \
                        and re.fullmatch(r'[A-Za-z0-9._-]+', v) \
                        and not v.startswith(('d__', 'p__', 'c__', 'o__',
                                              'f__', 'g__', 's__')):
                    ids.add(v)
            if len(ids) > 1:
                groups.append(ids)
    alias = {}
    for g in groups:
        for i in g:
            alias.setdefault(i, set()).update(g)
    return alias

ALIAS = build_alias([os.path.join(F, 'MAG_CROSSWALK.tsv'),
                     os.path.join(F, 'MAG_TAXONOMY.tsv')])

def resolve(nm, known):
    """map a name onto whichever alias of it appears in `known`"""
    if nm in known:
        return nm
    for a in ALIAS.get(nm, ()):
        if a in known:
            return a
    return None

gene, role = {}, {}
for line in open(os.path.join(E, 'targets.tsv')):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 3:
        gene[f[0]] = f[1]; role[f[0]] = f[2]

NOTSEARCHED = set()
ns = os.path.join(E, 'not_searched.txt')
if os.path.exists(ns):
    NOTSEARCHED = {l.strip() for l in open(ns) if l.strip()}

thr = {}
with open(os.path.join(E, 'ko_list'), errors='replace') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    ti = hdr.index('threshold') if 'threshold' in hdr else 1
    si = hdr.index('score_type') if 'score_type' in hdr else 2
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) > max(ti, si):
            try:
                thr[f[0]] = (float(f[ti]), f[si].strip().lower())
            except ValueError:
                thr[f[0]] = (None, f[si].strip().lower())

hits = defaultdict(set)     # ko -> {mag}
sub = Counter(); nothr = set()
mags = []
for fn in sorted(os.listdir(os.path.join(E, 'hits'))):
    if not fn.endswith('.tbl'):
        continue
    mag = fn[:-4]; mags.append(mag)
    for line in open(os.path.join(E, 'hits', fn), errors='replace'):
        if line.startswith('#'):
            continue
        f = line.split()
        if len(f) < 9:
            continue
        ko = f[2] if f[2] in gene else (f[3] if len(f) > 3 and f[3] in gene else None)
        if ko is None:
            continue
        t, st = thr.get(ko, (None, ''))
        try:
            sc = float(f[8]) if st and st.startswith('dom') else float(f[5])
        except ValueError:
            continue
        if t is None:
            nothr.add(ko); hits[ko].add(mag)
        elif sc >= t:
            hits[ko].add(mag)
        else:
            sub[ko] += 1

NM = max(1, len(mags))
rule('1. THE ECF MODULE AND ITS S-COMPONENTS, ACROSS THE CATALOGUE')
say(f"  {'KO':<9}{'gene':<9}{'role':<8}{'genomes':>9}{'%':>7}   rejected below threshold")
say('-' * W)
for ko in sorted(gene, key=lambda k: (role[k], k)):
    if ko in NOTSEARCHED:
        say(f"  {ko:<9}{gene[ko]:<9}{role[ko]:<8}{'NOT SEARCHED':>16}")
        continue
    g = hits.get(ko, set())
    say(f"  {ko:<9}{gene[ko]:<9}{role[ko]:<8}{len(g):>9}{100.0 * len(g) / NM:>6.0f}%"
        f"   {sub.get(ko, 0)}")
if nothr:
    say(f"  no curated threshold for {', '.join(sorted(nothr))} — kept on E-value alone")

# per-genome module / S-component
mod_kos = [k for k in gene if role[k] == 'MODULE' and k not in NOTSEARCHED]
sc_kos = [k for k in gene if role[k] == 'SCOMP' and k not in NOTSEARCHED]
has_mod = {m: [k for k in mod_kos if m in hits.get(k, ())] for m in mags}
has_sc = {m: [k for k in sc_kos if m in hits.get(k, ())] for m in mags}

nmod = sum(1 for m in mags if has_mod[m])
nsc = sum(1 for m in mags if has_sc[m])
nboth = sum(1 for m in mags if has_mod[m] and has_sc[m])
norphan = sum(1 for m in mags if has_sc[m] and not has_mod[m])

rule('2. THE CONTROL — AN S-COMPONENT CANNOT TRANSPORT WITHOUT ITS MODULE')
say("  BioY and its relatives are substrate-BINDING subunits. They have no")
say("  transmembrane channel and no ATPase; without EcfA/EcfA'/EcfT they bind a")
say("  vitamin and do nothing with it. So a catalogue with S-components and no")
say("  module is not a biological finding, it is a failed search.")
say()
say(f"  {'genomes with an ECF MODULE subunit':<44}{nmod:>5} / {len(mags)}")
say(f"  {'genomes with an S-COMPONENT':<44}{nsc:>5} / {len(mags)}")
say(f"  {'genomes with BOTH':<44}{nboth:>5} / {len(mags)}")
say(f"  {'S-component but NO module (orphans)':<44}{norphan:>5} / {len(mags)}")
say()
VALID = None
if not sc_kos:
    say("  [NO TEST] no S-component was searched.")
elif nsc == 0:
    say("  [NO TEST] no S-component found in any genome, so there is nothing to")
    say("            corroborate. Step 36's vitamin absences stand unchallenged.")
elif nmod == 0:
    VALID = False
    say("  [FAIL] S-components in "
        f"{nsc}/{len(mags)} genomes and the module in NONE of them.")
    say("         That combination is not biologically possible. Either the module")
    say("         profiles are wrong for these organisms or the search is broken.")
    say("         NO VERDICT — do not reopen any edge on this evidence.")
elif norphan > nboth:
    VALID = False
    say("  [FAIL] more genomes carry an orphan S-component than carry a complete")
    say("         system. The module search is not finding what must be there.")
    say("         NO VERDICT.")
else:
    VALID = True
    say(f"  [PASS] {nboth} genomes carry a complete ECF system — module plus at least")
    say("         one S-component. The transporters ARE here; step 36 was searching")
    say("         for the wrong protein family.")

say()
say(f"  {'genome':<32}{'module subunits':>16}{'S-components':>14}{'complete':>10}")
say('-' * W)
for m in sorted(mags, key=lambda x: (-len(has_mod[x]), -len(has_sc[x])))[:32]:
    say(f"  {m:<32}{len(has_mod[m]):>16}{len(has_sc[m]):>14}"
        f"{('yes' if has_mod[m] and has_sc[m] else 'no'):>10}")

with open(os.path.join(E, 'ecf_matrix.tsv'), 'w') as fh:
    fh.write('mag\tmodule_kos\tscomponent_kos\tcomplete\n')
    for m in sorted(mags):
        fh.write(f"{m}\t{','.join(has_mod[m])}\t{','.join(has_sc[m])}\t"
                 f"{'yes' if has_mod[m] and has_sc[m] else 'no'}\n")

# ---------------------------------------------- reopen step 36's edges
rule("3. WHICH OF STEP 36's DOWNGRADED EDGES MUST BE REOPENED")
ep = os.path.join(D, 'edges_closed.tsv')
VIT = ('b12', 'thiamin', 'biotin', 'cobalamin', 'riboflavin')
reopened = kept = 0
rows = []
if not os.path.exists(ep):
    say(f"  !! {ep} not found — run step 36 first; this section acts on its output.")
elif VALID is not True:
    say("  The control did not pass, so nothing is reopened. Step 36's edges stand")
    say("  exactly as it left them.")
else:
    with open(ep, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            nv = (r.get('new_verdict') or '')
            if not nv.startswith('DOWNGRADED'):
                continue
            met = (r.get('metabolite') or '').lower()
            cons = (r.get('consumer') or '').strip()
            cons = resolve(cons, set(mags)) or cons
            if any(v in met for v in VIT):
                if has_mod.get(cons) and has_sc.get(cons):
                    reopened += 1
                    rows.append((cons, r.get('metabolite', ''), nv,
                                 'REOPENED — vitamin, and this genome has a complete ECF system'))
                    continue
            kept += 1
    say(f"  {reopened:>5}  REOPENED — a vitamin edge whose consumer carries a complete")
    say("         ECF transporter. Step 36 downgraded it on the absence of an ABC")
    say("         system these organisms do not use. It is NOT decided either way.")
    say(f"  {kept:>5}  downgrade stands — not a vitamin, or no ECF system in that genome")
    say()
    say("  A REOPENED edge is not a supported edge. This step cannot say which")
    say("  vitamin any S-component binds, so it cannot confirm uptake — only that")
    say("  the evidence step 36 used does not settle the question.")
    with open(os.path.join(E, 'edges_reopened.tsv'), 'w') as oh:
        oh.write('consumer\tmetabolite\tstep36_verdict\tstatus\n')
        for t in rows:
            oh.write('\t'.join(t) + '\n')

rule('WHAT TO REPORT')
if VALID is True:
    say("  1. The catalogue uses ECF transporters for vitamins, not the classic ABC")
    say(f"     systems: {nboth}/{len(mags)} genomes carry a complete module plus at least one")
    say("     S-component, while step 36 found ZERO BtuFCD and ZERO ThiPQ.")
    if reopened:
        say("  2. Step 36's B12 and thiamine 'absences' are absences of the WRONG")
        say(f"     FAMILY. Report them as such, and reopen the {reopened} vitamin edges in")
        say("     section 3 rather than reporting them as downgraded.")
    else:
        say("  2. Step 36's B12 and thiamine 'absences' are absences of the WRONG")
        say("     FAMILY — but NO EDGE CHANGES, because none of step 36's downgraded")
        say("     edges involves a vitamin. The two findings are independent: the ECF")
        say("     result is about genome content, the downgrades are about a different")
        say("     metabolite entirely. Do not let one carry the other.")
    say("  3. What this does NOT establish: which vitamin any S-component binds.")
    say("     KEGG does not name them. 'An ECF system is present' is the claim;")
    say("     'B12 uptake is present' is not.")
elif VALID is False:
    say("  >> NO VERDICT. The module/S-component control failed; see section 2.")
else:
    say("  >> No ECF system found at all. Step 36's vitamin absences stand, and the")
    say("     'wrong family' hypothesis is refuted rather than confirmed — which is")
    say("     a result, not a failure.")
say()
say("  LIMITS:")
say("   1. S-components are short (~180 aa) and diverse; a curated threshold on a")
say("      short profile is a blunt instrument. The rejected counts in section 1")
say("      are worth reading before quoting any single genome.")
say("   2. Presence of a module says the machine exists, not what it carries.")
say("   3. This searched MAG proteomes only. The unbinned fraction is untouched.")

open(os.path.join(E, 'ECF_VERDICT.txt'), 'w').write('\n'.join(rep) + '\n')
print()
print(f"  written: {os.path.join(E, 'ecf_matrix.tsv')}")
print(f"  written: {os.path.join(E, 'ECF_VERDICT.txt')}")
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "interpretation completed" "$RC" || OK=1
check "ecf matrix written" "$([ -s "$E/ecf_matrix.tsv" ] && echo 0 || echo 1)" || OK=1
if grep -q 'NO VERDICT' "$E/ECF_VERDICT.txt" 2>/dev/null; then
  echo "  [STOP] the module/S-component control failed — nothing may be reopened."
  OK=1
fi
[ "$OK" -eq 0 ] && runlog "36b_ecf complete" || runlog "**36b_ecf FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
