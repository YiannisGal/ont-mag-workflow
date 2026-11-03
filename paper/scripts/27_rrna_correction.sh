#!/usr/bin/env bash
# =============================================================================
# 27_rrna_correction.sh  --  the correction step 23 asked for and nobody did
#
# WHAT STEP 23 SAID, VERBATIM FROM ITS OWN OUTPUT
#   "USE THIS TO CORRECT THE rRNA PROFILES. A genus with 5 operons contributes
#    5x the SSU reads per cell of a genus with 1."
#   It was never done. This does it.
#
# WHY IT MATTERS HERE MORE THAN USUAL
#   The spread in this catalogue is as wide as it gets: Lysinibacillus 10
#   operons, Chakrabartyella 9, Gaoshiqia 3, Methanofollis 0-1. An SSU read
#   count is a count of OPERONS, not of cells. Dividing by copy number converts
#   one to the other.
#
#   The direction is not neutral. THE METHANOGENS HAVE THE LOWEST COPY NUMBERS
#   IN THE CATALOGUE, so every rRNA-based number in this chapter understates
#   them — and methanogenesis is the chapter's subject.
#
# THREE LIMITS, ALL OF WHICH CONSTRAIN HOW THIS MAY BE REPORTED
#
#   1. COPY NUMBERS ARE FLOORS. rRNA operons are near-identical repeats that
#      collapse during assembly. A MAG's count is a minimum, so the correction
#      is a MINIMUM correction and the true effect is larger, never smaller.
#      That at least makes the bias one-directional and safe to state.
#
#   2. ONLY LINEAGES WITH A MAG CAN BE CORRECTED. This is the awkward one.
#      Correcting part of a profile and renormalising against an uncorrected
#      remainder produces a number that is neither raw nor corrected. So the
#      correction is computed and reported WITHIN the correctable subset, and
#      the uncorrected remainder is reported separately as what it is. A
#      whole-community correction needs rrnDB, an external download, and is not
#      attempted here on my own initiative.
#
#   3. TWO MAGs WITH ZERO DETECTED OPERONS cannot be corrected at all —
#      dividing by zero is not a small problem. They are excluded and named.
#
# OUTPUT   27_rrna_corrected/<SAMPLE>_corrected_profile.tsv
#          27_rrna_corrected/<SAMPLE>_rank_changes.tsv   (one pair per profile)
# RUNTIME  seconds.
# USAGE    ws2 && bash scripts/27_rrna_correction.sh 2>&1 | tee logs/27.log
#          SILVA_DIR=... to point at a profile other than the default
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="27_rrna_correction"
R="$WS2/27_rrna_corrected"; E="$WS2/23_extras"; F="$WS2/07_mags"
SILVA_DIR="${SILVA_DIR:-$WS2/13_profile/silva}"
mkdir -p "$R"

echo "######## 27 — rRNA COPY-NUMBER CORRECTION ########"
echo "  started: $(date -Is)"
runlog "27_rrna_correction started"
[ -s "$E/rrna_copy_number.tsv" ] || { echo "  !! run 23_genome_extras.sh first"; exit 1; }
[ -d "$SILVA_DIR" ] || { echo "  !! no SILVA profile at $SILVA_DIR (set SILVA_DIR)"; exit 1; }

# THE FIRST RUN CORRECTED ONLY ONE SAMPLE AND DID NOT SAY SO. It globbed the
# SILVA directory, took the first match, and reported LIQ_ssu_profile.tsv — the
# cultures. The sediment profile, which is the one the chapter is actually
# about, was never touched and nothing in the output revealed that. A script
# that silently picks one of several inputs is a script that will eventually be
# quoted as though it used all of them.
PROFILES=$(find "$SILVA_DIR" -maxdepth 2 -name '*ssu_profile.tsv' -size +0c 2>/dev/null | sort)
[ -n "$PROFILES" ] || { echo "  !! no *ssu_profile.tsv found under $SILVA_DIR"; exit 1; }
echo "  profiles to correct:"
echo "$PROFILES" | sed "s|$WS2/||; s|^|    |"

for P in $PROFILES; do
TAG=$(basename "$P" | sed 's/_ssu_profile\.tsv$//; s/\.tsv$//')
echo
echo "################ $TAG ################"
python3 - "$P" "$TAG" "$E" "$F" "$R" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict

prof_file, TAG, E, F, R = sys.argv[1:6]

# ---- taxonomy and copy number ---------------------------------------------
tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return ((tax.get(m, {}) or {}).get('genus','') or '')
def compl(m):
    try: return float((tax.get(m, {}) or {}).get('completeness','') or 0)
    except ValueError: return 0.0

# genus -> best copy number. Where several MAGs share a genus, take the MAXIMUM:
# every count is a floor, so the largest observed floor is the least wrong.
copies, src = {}, {}
for r in csv.DictReader(open(os.path.join(E, 'rrna_copy_number.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id') or ''
    m = wn if wn in tax else next((k for k in tax if mid(k) == wn), wn)
    g = (gen(m) or '').strip()
    if not g: continue
    # THE DIVISOR MUST BE THE 16S COUNT, NOT THE LARGEST OF THE THREE.
    # The first run took max(16S, 23S, 5S) on the reasoning that every count is
    # a floor so the biggest floor is least wrong. That is right for "how many
    # operons does this genome have" and WRONG here. This profile is built from
    # SSU reads, and what generates an SSU read is a 16S gene. Lysinibacillus
    # has 13 16S, 10 23S and 11 5S; for correcting 16S read counts the answer is
    # 13, and it was only right by luck. For a genome where 23S exceeded 16S the
    # max would have over-divided and understated that lineage.
    n = None
    for c in ('n_16S', '16S', 'ssu', 'n_ssu'):
        if r.get(c):
            try: n = int(float(r[c])); break
            except ValueError: pass
    if n is None:
        for c in ('floor', 'operons', 'copies'):
            if r.get(c):
                try: n = int(float(r[c])); break
                except ValueError: pass
    if n is None: continue
    if compl(m) < 90:      # an incomplete genome loses operons outright
        continue
    if g not in copies or n > copies[g]:
        copies[g] = n; src[g] = mid(m)

zero = sorted(g for g, n in copies.items() if n == 0)
usable = {g: n for g, n in copies.items() if n > 0}
print(f"  copy numbers from MAGs >=90% complete: {len(copies)} genera "
      f"({len(usable)} usable, {len(zero)} with zero detected operons)")
if zero:
    print(f"  EXCLUDED, cannot divide by zero: {', '.join(zero)}")
    print("   -> these are the genomes whose operon collapsed entirely during")
    print("      assembly. Their true copy number is >=1, so excluding them")
    print("      removes lineages that would have been corrected DOWNWARD.")

# ---- SILVA genus profile ---------------------------------------------------
rows = list(csv.DictReader(open(prof_file), delimiter='\t'))
if not rows: print("  !! profile is empty"); sys.exit(1)
hdr = list(rows[0])
gcol = next((h for h in hdr if 'genus' in h.lower()), None) or \
       next((h for h in hdr if 'taxon' in h.lower() or 'name' in h.lower()), hdr[0])
ncol = None
for h in hdr:
    hl = h.lower()
    if any(x in hl for x in ('read','count','n_','abund','pct','percent','fraction')):
        try:
            float(rows[0][h] or 0); ncol = h; break
        except (ValueError, TypeError): pass
if ncol is None:
    for h in hdr[1:]:
        try: float(rows[0][h] or 0); ncol = h; break
        except (ValueError, TypeError): pass
if ncol is None: print("  !! no numeric column in the profile"); sys.exit(1)
print(f"  columns: genus='{gcol}'  value='{ncol}'")

raw = {}
for r in rows:
    g = re.sub(r'^[a-z]__', '', (r.get(gcol) or '').strip())
    if not g: continue
    try: v = float(r.get(ncol) or 0)
    except ValueError: continue
    raw[g] = raw.get(g, 0.0) + v

# ---- correct within the correctable subset --------------------------------
matched = {g: v for g, v in raw.items() if g in usable}
if not matched:
    print("  !! no SILVA genus matched a MAG genus — nothing to correct"); sys.exit(1)
tot_raw_all = sum(raw.values()) or 1.0
tot_raw_sub = sum(matched.values()) or 1.0
corr = {g: v / usable[g] for g, v in matched.items()}
tot_cor_sub = sum(corr.values()) or 1.0

print()
print("=" * 104)
print("  CORRECTION WITHIN THE SUBSET THAT HAS A MAG")
print("=" * 104)
print(f"  {len(matched)} of {len(raw)} genera in the profile have a MAG with a")
print(f"  usable operon count. They are {100*tot_raw_sub/tot_raw_all:.1f}% of the "
      f"profile's reads.")
print(f"  The other {100*(1-tot_raw_sub/tot_raw_all):.1f}% CANNOT be corrected and are excluded")
print("  from the percentages below. These are shares WITHIN the correctable")
print("  subset — they are not community percentages and must never be quoted")
print("  as if they were.")
print()
print(f"  {'genus':<26}{'operons':>8}{'raw %':>9}{'corrected %':>13}{'change':>9}  source MAG")
out = []
for g in sorted(matched, key=lambda z: -corr[z]):
    rp = 100 * matched[g] / tot_raw_sub
    cp = 100 * corr[g] / tot_cor_sub
    print(f"  {g[:25]:<26}{usable[g]:>8}{rp:>9.2f}{cp:>13.2f}{cp-rp:>+9.2f}  {src.get(g,'')}")
    out.append(dict(genus=g, operons=usable[g], raw_reads=f"{matched[g]:.4f}",
                    raw_pct_of_subset=f"{rp:.3f}", corrected_pct_of_subset=f"{cp:.3f}",
                    change_pct_points=f"{cp-rp:+.3f}", source_mag=src.get(g, '')))
with open(os.path.join(R, f'{TAG}_corrected_profile.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(out[0]))
    w.writeheader(); w.writerows(out)

# ---- rank changes ----------------------------------------------------------
rr = sorted(matched, key=lambda z: -matched[z])
rc = sorted(matched, key=lambda z: -corr[z])
moves = [(g, rr.index(g) + 1, rc.index(g) + 1) for g in matched]
moves.sort(key=lambda t: -(t[1] - t[2]))
print()
print("=" * 104)
print("  WHAT THE CORRECTION CHANGES")
print("=" * 104)
big = [t for t in moves if abs(t[1] - t[2]) >= 2]
if big:
    print(f"  {'genus':<26}{'raw rank':>10}{'corrected':>11}{'move':>7}")
    for g, a, b in big:
        print(f"  {g[:25]:<26}{a:>10}{b:>11}{a-b:>+7}")
else:
    print("  No genus moves more than one rank. The correction is real but does")
    print("  not reorder this profile, which is itself worth stating: copy-number")
    print("  bias exists here and does not change the headline.")
with open(os.path.join(R, f'{TAG}_rank_changes.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['genus','raw_rank','corrected_rank','move'])
    for g, a, b in moves: w.writerow([g, a, b, a - b])

meth = [g for g in matched if any(x in g.lower() for x in ('methano','methan'))]
if meth:
    rp = sum(100 * matched[g] / tot_raw_sub for g in meth)
    cp = sum(100 * corr[g] / tot_cor_sub for g in meth)
    print()
    print(f"  METHANOGENS: {rp:.2f}% of the correctable subset raw, {cp:.2f}% corrected "
          f"({cp-rp:+.2f} points).")
    if cp > rp:
        print("  They rise, as predicted from their operon counts being the lowest in")
        print("  the catalogue. An rRNA profile systematically understates them, and")
        print("  that is a methodological point worth making in its own right —")
        print("  it applies to every amplicon survey of a methanogenic system.")

print()
print("  HOW TO REPORT THIS")
print("  Quote raw percentages for the community and corrected percentages only")
print("  within the subset, saying which. The correction uses FLOOR copy numbers,")
print("  so it is a MINIMUM correction: the real effect is in the same direction")
print("  and larger. That one-sidedness is what makes it safe to state at all.")
PY
done

echo
echo "######## VALIDATION GATE ########"
NC=$(ls "$R"/*_corrected_profile.tsv 2>/dev/null | wc -l)
NR=$(ls "$R"/*_rank_changes.tsv 2>/dev/null | wc -l)
NP=$(echo "$PROFILES" | wc -l)
echo "  profiles in: $NP   corrected tables out: $NC   rank tables out: $NR"
check "one corrected table per profile" $([ "$NC" -eq "$NP" ] && echo 0 || echo 1)
check "one rank table per profile"      $([ "$NR" -eq "$NP" ] && echo 0 || echo 1)
runlog "27_rrna_correction complete"
echo
echo "######## 27 DONE ########"
