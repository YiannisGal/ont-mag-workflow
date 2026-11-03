#!/usr/bin/env bash
# =============================================================================
# 13c_silva_report.sh  --  re-read the SILVA profile from disk and report it
#
# WHY THIS EXISTS
#   The step-13b output was printed to a terminal and to a log, and the strict
#   re-run overwrote the permissive one. Neither of those is a result you can
#   cite. This script reads the profile TSVs that 13b already wrote and produces
#   the tables that actually go in the thesis, so the numbers can be regenerated
#   at any time without re-running a 1-3 h mapping job.
#
# THE IDENTITY CUTOFF, AND WHY THE FIRST RUN WAS WRONG
#   13b ran first at MINID=0.80. That was chosen for ONT error rate and it is
#   wrong for rRNA, because the SSU gene is CONSERVED: at 80% identity a read
#   matches the conserved core of almost any bacterium, and the "best hit" is
#   then decided by database composition rather than by the organism. The
#   symptom was unmistakable — *Streptococcus* at 33.49% of an anoxic, sulfidic
#   lagoon sediment. Streptococcus is a host-associated genus with enormous
#   database representation and no business dominating this sample.
#
#   Re-run at MINID=0.92 / MINLEN=1000 it dropped out of the top hits. Record
#   both runs: the failure is itself a methods point worth one sentence, because
#   the same trap catches anyone mapping long reads to rRNA databases.
#
# WHAT THIS PRINTS
#   1. domain composition per sample
#   2. top genera, as relative rRNA abundance (NOT cell abundance)
#   3. THE KEY TABLE — abundant lineages with no MAG, i.e. the named gaps in the
#      catalogue
#   4. STRUCTURAL screens: single-genus dominance, and facultative-aerobe load
#
# USAGE   ws2 && bash scripts/13c_silva_report.sh 2>&1 | tee logs/13c_silva_report.log
#         SILVA_DIR=/path/to/other/run bash scripts/13c_silva_report.sh
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="13c_silva_report"
S="${SILVA_DIR:-$WS2/13_profile/silva}"; F="$WS2/07_mags"

echo "######## 13c — SILVA SSU REPORT ########"
echo "  reading: $S"
[ -d "$S" ] || { echo "  !! no SILVA output at $S — run 13b_silva_ssu.sh first"; exit 1; }
ls -1 "$S" | sed 's/^/    /'

python3 - "$S" "$F" <<'PY'
import sys, os, csv
from collections import defaultdict

S, F = sys.argv[1], sys.argv[2]

# STRUCTURAL implausibility screen.
# The first version matched genus names against a hand list by EXACT string. It
# missed "Escherichia-Shigella" (list had "escherichia"), and Erwinia,
# Citrobacter and Enterobacter were never on the list at all — so it printed
# "consistent with a genuine environmental profile" over a sediment profile that
# is >90% Enterobacterales. A screen built from remembered names will always
# have that hole. These two screens do not depend on my recall:
#   1. DOMINANCE  — any single genus above 20% of an environmental sample
#   2. AEROBE LOAD — orders that cannot respire in permanently anoxic sediment
AEROBIC_ORDERS = ('enterobacterales','pseudomonadales','burkholderiales',
                  'xanthomonadales','rhizobiales','sphingomonadales',
                  'caulobacterales','micrococcales','corynebacteriales',
                  'lactobacillales','staphylococcales','neisseriales',
                  'rickettsiales','acetobacterales','flavobacteriales')
AEROBIC_FAMILIES = ('enterobacteriaceae','erwiniaceae','yersiniaceae',
                    'morganellaceae','pectobacteriaceae','pseudomonadaceae',
                    'moraxellaceae','streptococcaceae','staphylococcaceae',
                    'lactobacillaceae','micrococcaceae','acetobacteraceae')
DOMINANCE_LIMIT = 20.0

mag_taxa = set()
p = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter='\t'):
        for k in ('phylum','class','order','family','genus'):
            v = (r.get(k) or '').strip()
            if v: mag_taxa.add(v.lower())
print(f"\n  taxa represented in the 29-MAG catalogue: {len(mag_taxa)}")

found_any = False
for code in ('SED','LIQ'):
    fp = os.path.join(S, f'{code}_ssu_profile.tsv')
    if not os.path.exists(fp):
        print(f"\n  {code}: no profile at {fp}")
        continue
    found_any = True
    rows = list(csv.DictReader(open(fp), delimiter='\t'))
    by_rank = defaultdict(list)
    for r in rows:
        try: by_rank[r['rank']].append((r['taxon'], int(r['reads']), float(r['pct_of_assigned'])))
        except (ValueError, KeyError): pass

    print("\n" + "=" * 92)
    print(f"  {code}  —  SSU rRNA PROFILE")
    print("=" * 92)

    dom = by_rank.get('domain', [])
    dtot = sum(n for _, n, _ in dom) or 1
    print(f"  reads assigned at domain level: {dtot:,}")
    print("  domain composition: " + ', '.join(f"{t} {100*n/dtot:.1f}%"
          for t, n, _ in sorted(dom, key=lambda x: -x[1])[:5]))

    gen = sorted(by_rank.get('genus', []), key=lambda x: -x[1])
    gtot = sum(n for _, n, _ in gen) or 1
    print(f"\n  TOP 25 GENERA (% of genus-assigned SSU reads)")
    print(f"  {'%':>7}  {'reads':>9}  {'MAG?':<6} genus")
    for t, n, _ in gen[:25]:
        inmag = 'yes' if t.lower() in mag_taxa else '—'
        flag = ''
        tl = t.lower()
        if 100*n/gtot > DOMINANCE_LIMIT: flag = '  <-- DOMINANT, verify'
        elif any(a[:-3] in tl for a in ('erwiniaceae','enterobacteriaceae')): flag = '  <-- aerobe'
        print(f"  {100*n/gtot:>6.2f}%  {n:>9,}  {inmag:<6} {t}{flag}")

    covered = sum(n for t, n, _ in gen if t.lower() in mag_taxa)
    print(f"\n  {100*covered/gtot:.1f}% of genus-assigned SSU reads belong to genera that")
    print(f"  the MAG catalogue contains. {100-100*covered/gtot:.1f}% do not.")

    missing = [(t, n) for t, n, _ in gen if t.lower() not in mag_taxa][:20]
    print(f"\n  ABUNDANT LINEAGES WITH NO MAG — the named gaps in the catalogue")
    print(f"  {'%':>7}  {'reads':>9}  genus")
    for t, n in missing:
        print(f"  {100*n/gtot:>6.2f}%  {n:>9,}  {t}")

    # ---- screen 1: dominance
    print()
    top_g, top_n, _ = gen[0] if gen else ('', 0, 0)
    top_pct = 100*top_n/gtot
    print(f"  SCREEN 1 — DOMINANCE: top genus is {top_g} at {top_pct:.2f}%")
    # An enrichment culture is SUPPOSED to be dominated — that is what
    # enrichment means. Applying an environmental-diversity expectation to LIQ
    # produced a false alarm on Anaerotignum at 21%. Only SED is environmental.
    is_env = (code == 'SED')
    dom_fail = is_env and top_pct > DOMINANCE_LIMIT
    if not is_env and top_pct > DOMINANCE_LIMIT:
        print("  >> This is an ENRICHMENT CULTURE. Dominance is the expected outcome")
        print("     of enrichment, not a warning sign. Screen not applied.")
    if dom_fail:
        print(f"  >> ABOVE {DOMINANCE_LIMIT:.0f}%. A single genus dominating an environmental")
        print("     sample this heavily is rare and usually indicates contamination, a")
        print("     reference sink, or a mapping artefact. VERIFY before reporting.")
    else:
        print("  >> No single genus dominates — normal for an environmental sample.")

    # ---- screen 2: aerobe load, by order/family, substring not exact match
    aer = 0.0
    for t, n, _ in gen:
        tl = t.lower()
        if any(a[:-1] in tl for a in AEROBIC_FAMILIES) or tl in AEROBIC_FAMILIES:
            aer += 100*n/gtot
    # genus names alone rarely carry the family, so also fold in the family rank
    fam = sorted(by_rank.get('family', []), key=lambda x: -x[1])
    ftot = sum(n for _, n, _ in fam) or 1
    aer_f = sum(100*n/ftot for t, n, _ in fam if t.lower() in AEROBIC_FAMILIES)
    ord_ = sorted(by_rank.get('order', []), key=lambda x: -x[1])
    otot = sum(n for _, n, _ in ord_) or 1
    aer_o = sum(100*n/otot for t, n, _ in ord_ if t.lower() in AEROBIC_ORDERS)
    print()
    print(f"  SCREEN 2 — AEROBE LOAD (organisms that cannot respire in this habitat)")
    print(f"    at family rank: {aer_f:.1f}%    at order rank: {aer_o:.1f}%")
    if fam:
        print(f"    top families: " + ', '.join(f"{t} {100*n/ftot:.1f}%"
              for t, n, _ in fam[:5]))
    aer_fail = max(aer_f, aer_o) > 25
    if aer_fail:
        print("  >> Etoliko is permanently anoxic and sulfidic. Facultative aerobes at")
        print("     this level did not grow there. Treat this profile as UNVERIFIED and")
        print("     run 13d_silva_verify.sh before using any of it.")
    else:
        print("  >> Consistent with a permanently anoxic habitat.")

    if dom_fail or aer_fail:
        print()
        print(f"  ***  {code} PROFILE FAILS SCREENING — DO NOT REPORT IT YET  ***")
        print("  ***  run: bash scripts/13d_silva_verify.sh                   ***")

if not found_any:
    print("\n  No profile tables found. Run 13b_silva_ssu.sh first:")
    print("    MINID=0.92 MINLEN=1000 bash scripts/13b_silva_ssu.sh")
    raise SystemExit(1)

print()
print("=" * 92)
print("  HOW TO REPORT THIS")
print("=" * 92)
print("  'SSU rRNA reads were recovered from the unassembled metagenome by mapping")
print("   to SILVA SSURef NR99 (minimap2 map-ont, >=92% identity over >=1000 bp).")
print("   Relative rRNA abundance is reported; because rRNA operon copy number")
print("   varies roughly 1-15 between taxa, these are not cell abundances. The same")
print("   bias applies to the 16S amplicon data, which is why the two datasets are")
print("   directly comparable to each other.'")
print()
print("  AND THE ONE SENTENCE ABOUT THE CUTOFF — worth including, because it is a")
print("  trap others will fall into:")
print("  'An initial identity threshold of 80%, chosen for Nanopore error rate,")
print("   produced implausible host-associated taxa at high abundance; the SSU gene")
print("   is sufficiently conserved that a permissive threshold assigns reads by")
print("   database representation rather than by organism. The threshold was raised")
print("   to 92% over a 1000 bp minimum alignment.'")
PY

echo
echo "######## 13c DONE ########"
