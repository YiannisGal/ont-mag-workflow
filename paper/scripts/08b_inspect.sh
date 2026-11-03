#!/usr/bin/env bash
# =============================================================================
# 08b_inspect.sh  --  read MAG_QC.tsv properly before making two decisions
#
# DECISION 1: did polishing fix the gene-length problem?
#   The 244 aa baseline is workspace_1's median over its 38 RETAINED MAGs — a
#   quality-filtered set. Comparing it to a median over all 47 union bins (which
#   include ~16 low-quality ones) is not like-for-like: low-completeness bins are
#   fragmented, fragmentation truncates ORFs, and those bins pull the median
#   down. This script reports the median over ALL bins, over HQ+MQ only, and per
#   sample, so the comparison is made on matched sets.
#
# DECISION 2: which DAS_Tool score threshold?
#   The step-07b rule prefers the stricter threshold when the looser one adds
#   chimeric bins. But a chimera already classed 'low' would be discarded anyway
#   and costs nothing. What matters is the count of NON-CHIMERIC HQ+MQ bins —
#   the actual publishable MAG set. That is recomputed here.
#
# READ-ONLY. Runs in seconds. Changes nothing.
# USAGE   ws2 && bash scripts/08b_inspect.sh 2>&1 | tee logs/08b_inspect.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
Q="$WS2/06_binqc"
[ -s "$Q/MAG_QC.tsv" ] || { echo "!! $Q/MAG_QC.tsv not found — run 08 first"; exit 1; }

python3 - "$Q/MAG_QC.tsv" <<'PY'
import sys, csv, statistics as st

rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
def fl(r, k):
    try: return float(r[k])
    except: return None

SCORES = ['0.5', '0.4', '0.3']
BASE = 244.0

def med(vals):
    vals = [v for v in vals if v is not None]
    return st.median(vals) if vals else float('nan')

# ====================================================== 1. GENE LENGTH
print("=" * 72)
print("  AVERAGE GENE LENGTH — matched-set comparison")
print("=" * 72)
print(f"  workspace_1 baseline: {BASE:.0f} aa (median of 38 RETAINED MAGs)")
print(f"  expected prokaryotic range: 290-320 aa\n")
print(f"  {'set':<34} {'n':>4} {'median':>8} {'min':>6} {'max':>6} {'>=290':>7}")

def line(label, sub):
    v = [fl(r, 'avg_gene_length') for r in sub]
    v = [x for x in v if x is not None]
    if not v:
        print(f"  {label:<34} {0:>4}"); return
    n290 = sum(1 for x in v if x >= 290)
    print(f"  {label:<34} {len(v):>4} {st.median(v):>8.1f} {min(v):>6.0f} "
          f"{max(v):>6.0f} {n290:>4}/{len(v)}")

hqmq = [r for r in rows if r['mimag'] in ('high', 'medium')]
clean = [r for r in hqmq if r['gunc_chimeric'] != 'yes']
line("ALL union bins", rows)
line("  MIMAG high+medium", hqmq)
line("    ...and non-chimeric", clean)
line("  MIMAG low only", [r for r in rows if r['mimag'] == 'low'])
print()
line("SED — all", [r for r in rows if r['sample'] == 'SED'])
line("SED — high+medium", [r for r in hqmq if r['sample'] == 'SED'])
line("LIQ — all", [r for r in rows if r['sample'] == 'LIQ'])
line("LIQ — high+medium", [r for r in hqmq if r['sample'] == 'LIQ'])

m_all, m_hq = med([fl(r,'avg_gene_length') for r in rows]), \
              med([fl(r,'avg_gene_length') for r in hqmq])
print()
print(f"  LIKE-FOR-LIKE VERDICT (HQ+MQ set vs workspace_1's retained set):")
print(f"    {BASE:.0f} aa  ->  {m_hq:.1f} aa    ({m_hq-BASE:+.1f} aa, "
      f"{100*(m_hq-BASE)/BASE:+.1f}%)")
if   m_hq >= 290: print("    >> CONFIRMED — within the expected prokaryotic range.")
elif m_hq >= 270: print("    >> PARTIAL — large improvement, still short of 290 aa.")
else:             print("    >> NOT CONFIRMED on the matched set.")

# is residual shortness explained by fragmentation / completeness?
print()
print("  Is the remaining shortfall explained by bin quality?")
print(f"  {'completeness band':<22} {'n':>4} {'med AGL':>9} {'med coding density':>19}")
bands = [(90,101,'>=90%'),(70,90,'70-90%'),(50,70,'50-70%'),(0,50,'<50%')]
for lo, hi, lab in bands:
    sub = [r for r in rows if (fl(r,'completeness') or -1) >= lo
                           and (fl(r,'completeness') or -1) < hi]
    if not sub: continue
    print(f"  {lab:<22} {len(sub):>4} {med([fl(r,'avg_gene_length') for r in sub]):>9.1f} "
          f"{med([fl(r,'coding_density') for r in sub]):>19.3f}")
print("  (If AGL rises steadily with completeness, the residual gap is")
print("   fragmentation in poor bins, not uncorrected sequencing error.)")

# ====================================================== 2. CHIMERAS
print()
print("=" * 72)
print("  GUNC CHIMERIC BINS — do they cost us anything?")
print("=" * 72)
chim = [r for r in rows if r['gunc_chimeric'] == 'yes']
if not chim:
    print("  none flagged")
else:
    print(f"  {'bin':<26} {'smp':<4} {'MIMAG':<7} {'comp':>6} {'cont':>6} "
          f"{'CSS':>6}  {'s0.5':>5} {'s0.4':>5} {'s0.3':>5}")
    for r in sorted(chim, key=lambda r: (r['sample'], r['bin'])):
        print(f"  {r['bin']:<26} {r['sample']:<4} {r['mimag']:<7} "
              f"{r['completeness'][:6]:>6} {r['contamination'][:6]:>6} "
              f"{r['gunc_css'][:6]:>6}  "
              f"{r.get('in_s0.5',''):>5} {r.get('in_s0.4',''):>5} {r.get('in_s0.3',''):>5}")
    print()
    print("  A chimera already classed 'low' is discarded by the MIMAG filter")
    print("  anyway and costs nothing. One classed high/medium is the real")
    print("  problem: it would enter the MAG set looking clean to CheckM2.")

# ====================================================== 3. THRESHOLD
print()
print("=" * 72)
print("  THRESHOLD DECISION — on the PUBLISHABLE set (HQ+MQ, non-chimeric)")
print("=" * 72)
print(f"  {'sample':<7} {'score':>5} {'high':>5} {'medium':>7} {'publishable':>12} "
      f"{'chimeric HQ/MQ':>15} {'discarded low':>14}")
for sample in ('SED', 'LIQ'):
    for s in SCORES:
        sub = [r for r in rows if r['sample'] == sample and r.get(f'in_s{s}') == 'yes']
        ok  = [r for r in sub if r['mimag'] in ('high','medium')
                              and r['gunc_chimeric'] != 'yes']
        hi  = sum(1 for r in ok if r['mimag'] == 'high')
        me  = sum(1 for r in ok if r['mimag'] == 'medium')
        cx  = sum(1 for r in sub if r['mimag'] in ('high','medium')
                                 and r['gunc_chimeric'] == 'yes')
        lo  = sum(1 for r in sub if r['mimag'] == 'low')
        print(f"  {sample:<7} {s:>5} {hi:>5} {me:>7} {len(ok):>12} {cx:>15} {lo:>14}")
    print()
print("  Read the 'publishable' column only. 'discarded low' bins cost nothing")
print("  but disk. 'chimeric HQ/MQ' is the column that should veto a threshold.")

# ====================================================== 4. THE MAG SET
print()
print("=" * 72)
print("  CANDIDATE MAG SET (HQ+MQ, non-chimeric, any threshold)")
print("=" * 72)
print(f"  {'bin':<26} {'smp':<4} {'binner':<8} {'MIMAG':<7} {'comp':>6} {'cont':>6} "
      f"{'AGL':>6} {'CD':>6} {'tRNA':>5} {'rRNA':>6}")
for r in sorted(clean, key=lambda r: (r['sample'], -(fl(r,'completeness') or 0))):
    rr = ''.join(x for x, k in (('5', 'has_5S'), ('6', 'has_16S'), ('3', 'has_23S'))
                 if r[k] == 'yes') or '-'
    print(f"  {r['bin']:<26} {r['sample']:<4} {r['binner']:<8} {r['mimag']:<7} "
          f"{(r['completeness'] or '')[:6]:>6} {(r['contamination'] or '')[:6]:>6} "
          f"{(r['avg_gene_length'] or '')[:6]:>6} {(r['coding_density'] or '')[:6]:>6} "
          f"{r['n_trna_types']:>5} {rr:>6}")
print()
print(f"  TOTAL candidate MAGs: {len(clean)}   "
      f"(workspace_1 published: 22 in the tree, 38 in the catalogue)")
PY
