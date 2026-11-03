#!/usr/bin/env bash
# =============================================================================
# 08_binqc.sh  --  CheckM2 + GUNC + MIMAG completeness, on the union bin set
#
# PURPOSE — two jobs in one pass.
#
#   JOB 1: THE HEADLINE NUMBER.
#     CheckM2's Average_Gene_Length column. workspace_1's 38 MAGs had a MEDIAN
#     of 244 aa against an expected 290-320 aa for prokaryotes, and a single
#     closed 4.90 Mb contig showed 281 aa — proving fragmentation was not the
#     cause. The diagnosis was frameshifts from Flye's --nano-raw error model
#     applied to Q20+ sup reads, compounded by no Medaka polishing. workspace_2
#     fixed both. This column is the test. It is reported whether or not it
#     confirms the hypothesis.
#
#   JOB 2: THE THRESHOLD DECISION.
#     Step 07b showed SED 0.4 == 0.5 exactly (6 bins) while 0.3 adds 4, and LIQ
#     rises 25 -> 33 -> 37 but reintroduces a 12.23 Mb bin at 0.4 and below.
#     Bin counts cannot settle this. What settles it is how many ADDITIONAL bins
#     at each threshold meet MIMAG medium/high quality — a decision rule fixed
#     BEFORE seeing the result (recorded in the lab notebook at step 07b).
#
#     Rather than run CheckM2 three times, the three threshold sets are merged
#     into a non-redundant union (deduplicated on exact contig membership),
#     assessed once, and each threshold's yield is derived by lookup.
#
# QC APPLIED
#   CheckM2      completeness / contamination / Average_Gene_Length / coding density
#   GUNC         chimerism — detects lineage mixing that CheckM2 is blind to,
#                because a chimera assembled from non-marker regions duplicates
#                no single-copy genes and so reads as low contamination
#   barrnap      5S / 16S / 23S rRNA presence (MIMAG high-quality requirement)
#   tRNAscan-SE  distinct tRNA types (MIMAG high-quality requires >=18)
#
# MIMAG (Bowers et al. 2017)
#   High     >=90% complete, <5% contaminated, 5S+16S+23S present, >=18 tRNAs
#   Medium   >=50% complete, <10% contaminated
#   Low      everything else
#
# RUNTIME  2-5 h. GUNC against proGenomes 2.1 (12 Gbp) dominates. Use tmux.
#
# USAGE
#   ws2 && bash scripts/08_binqc.sh 2>&1 | tee logs/08_binqc_$(date +%Y%m%d_%H%M).log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="08_binqc"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
C="$WS2/05_consensus"; Q="$WS2/06_binqc"
U="$Q/union_bins"
mkdir -p "$Q" "$U"

SCORES="0.5 0.4 0.3"     # 0.5 FIRST: the primary set's bin names are canonical

echo "######## 08 — BIN QC (CheckM2 / GUNC / MIMAG) ########"
echo "  started: $(date -Is)  threads: $THREADS"
runlog "08_binqc started"

# ------------------------------------------------------------- DB PREFLIGHT
# 00_config.sh names the CheckM2 database CHECKM2DB (no underscore); accept both.
export CHECKM2_DB="${CHECKM2_DB:-${CHECKM2DB:-}}"

fail=0
for v in CHECKM2_DB GUNC_DB; do
  p="${!v:-}"
  if [ -z "$p" ]; then echo "  !! \$$v is not set in 00_config.sh"; fail=1
  elif [ ! -e "$p" ]; then echo "  !! \$$v points at a missing path: $p"; fail=1
  else echo "  [db] $v = $p"; fi
done
[ "$fail" -eq 0 ] || { echo "  Fix the database paths in 00_config.sh and re-run."; exit 1; }

record_version "checkm2"     "${ENV_checkm2:-checkm2}"       checkm2 --version
record_version "gunc"        "${ENV_gunc:-gunc_env}"         gunc --version
record_version "barrnap"     "${ENV_barrnap:-NB_Genome}"     barrnap --version
record_version "tRNAscan-SE" "${ENV_tRNAscan_SE:-NB_Genome}" tRNAscan-SE -h

# ==========================================================================
# 1. BUILD THE NON-REDUNDANT UNION OF THE THREE THRESHOLD SETS
# ==========================================================================
echo
echo "######## BUILDING UNION BIN SET ########"
if [ -s "$Q/union_membership.tsv" ] && ls "$U"/*.fa >/dev/null 2>&1; then
  echo "  union already built — delete $Q/union_membership.tsv to force"
else
  rm -f "$U"/*.fa
  python3 - "$C" "$U" "$Q" $SCORES <<'PY'
import sys, os, glob, hashlib, csv, re

consensus, udir, qdir = sys.argv[1], sys.argv[2], sys.argv[3]
scores = sys.argv[4:]

def contig_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                ids.append(line[1:].split()[0])
    return ids

def bin_dirs(score):
    """0.5 primary set lives in 05_consensus/bins; sweeps in 05_consensus/sweep."""
    out = {}
    for s in ("SED", "LIQ"):
        d = os.path.join(consensus, "sweep", f"{s}_s{score}_DASTool_bins")
        if not os.path.isdir(d):
            d = os.path.join(consensus, "bins")       # fallback: primary run
        out[s] = d
    return out

union = {}          # key -> dict(name, sample, binner, path, contigs)
membership = {}     # key -> set(scores)
used_names = set()

for score in scores:                       # 0.5 first => its names win
    dirs = bin_dirs(score)
    for sample, d in dirs.items():
        for f in sorted(glob.glob(os.path.join(d, "*.fa"))):
            base = os.path.basename(f)[:-3]
            if not base.startswith(sample + "_"):
                continue                    # primary dir holds both samples
            ids = contig_ids(f)
            if not ids:
                continue
            key = hashlib.md5(("\n".join(sorted(ids))).encode()).hexdigest()
            membership.setdefault(key, set()).add(score)
            if key in union:
                continue
            m = re.match(r'^(SED|LIQ)_(?:das_)?(metabat|vamb|semibin)_(\d+)$', base)
            binner = m.group(2) if m else "unknown"
            num = m.group(3) if m else "000"
            name = f"{sample}_das_{binner}_{num}"
            suffix = 0
            while name in used_names:       # same name, different content
                suffix += 1
                name = f"{sample}_das_{binner}_{num}{chr(96+suffix)}"
            used_names.add(name)
            union[key] = dict(name=name, sample=sample, binner=binner,
                              path=f, contigs=ids)

for key, b in union.items():
    dst = os.path.join(udir, b["name"] + ".fa")
    with open(b["path"]) as src, open(dst, "w") as out:
        out.write(src.read())

with open(os.path.join(qdir, "union_membership.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["bin", "sample", "binner", "n_contigs", "bp"]
               + [f"in_s{s}" for s in scores] + ["source_file"])
    for key, b in sorted(union.items(), key=lambda kv: kv[1]["name"]):
        bp = sum(len(l.strip()) for l in open(b["path"]) if not l.startswith('>'))
        w.writerow([b["name"], b["sample"], b["binner"], len(b["contigs"]), bp]
                   + [("yes" if s in membership[key] else "no") for s in scores]
                   + [b["path"]])

print(f"  union: {len(union)} distinct bins across thresholds {', '.join(scores)}")
for s in scores:
    n = sum(1 for k in membership if s in membership[k])
    print(f"    score {s}: {n} bins")
PY
fi
NBINS=$(ls -1 "$U"/*.fa 2>/dev/null | wc -l)
echo "  union bin FASTAs written: $NBINS"
[ "$NBINS" -gt 0 ] || { echo "  !! no union bins — check 05_consensus/"; exit 1; }

# ==========================================================================
# 2. CHECKM2
# ==========================================================================
echo
echo "######## CheckM2 ($(date +%H:%M)) ########"
if [ -s "$Q/checkm2/quality_report.tsv" ]; then
  echo "  already done"
else
  run_tool checkm2 predict --threads "$THREADS" --input "$U" -x fa \
      --output-directory "$Q/checkm2" --database_path "$CHECKM2_DB" --force \
    || echo "  !! CheckM2 failed — see log"
fi

# ==========================================================================
# 3. GUNC  (chimerism — CheckM2 cannot see this class of error)
# ==========================================================================
echo
echo "######## GUNC ($(date +%H:%M)) — the long one ########"
mkdir -p "$Q/gunc"
if ls "$Q/gunc"/GUNC*maxCSS*.tsv >/dev/null 2>&1; then
  echo "  already done"
else
  run_tool gunc run --input_dir "$U" --file_suffix .fa --db_file "$GUNC_DB" \
      --out_dir "$Q/gunc" --threads "$THREADS" --temp_dir "$TMPDIR" \
    || echo "  !! GUNC failed — see log; MIMAG will proceed without chimerism"
fi

# ==========================================================================
# 4. rRNA + tRNA  (MIMAG high-quality requirements)
# ==========================================================================
echo
echo "######## barrnap + tRNAscan-SE ($(date +%H:%M)) ########"
mkdir -p "$Q/rrna" "$Q/trna"
printf 'bin\tn_trna_types\thas_5S\thas_16S\thas_23S\n' > "$Q/rna_summary.tsv"
i=0
for f in "$U"/*.fa; do
  b=$(basename "$f" .fa); i=$((i+1))
  printf '\r  [%3d/%3d] %s                    ' "$i" "$NBINS" "$b"

  if [ ! -s "$Q/rrna/${b}.gff" ]; then
    { run_tool barrnap --kingdom bac --threads 4 --quiet "$f"
      run_tool barrnap --kingdom arc --threads 4 --quiet "$f"; } > "$Q/rrna/${b}.gff" 2>/dev/null
  fi
  h5=$(grep -c '5S_rRNA'  "$Q/rrna/${b}.gff" 2>/dev/null || true)
  h16=$(grep -c '16S_rRNA' "$Q/rrna/${b}.gff" 2>/dev/null || true)
  h23=$(grep -c '23S_rRNA' "$Q/rrna/${b}.gff" 2>/dev/null || true)

  # -G = general/domain-agnostic: domain is not known until GTDB-Tk at step 09
  if [ ! -s "$Q/trna/${b}.txt" ]; then
    run_tool tRNAscan-SE -G -q -o "$Q/trna/${b}.txt" "$f" >/dev/null 2>&1
  fi
  nt=$(awk 'NR>3 && $5!="Undet" && $5!="Sup" {print $5}' "$Q/trna/${b}.txt" 2>/dev/null \
       | sort -u | wc -l)

  printf '%s\t%s\t%s\t%s\t%s\n' "$b" "${nt:-0}" \
    "$([ "${h5:-0}"  -gt 0 ] && echo yes || echo no)" \
    "$([ "${h16:-0}" -gt 0 ] && echo yes || echo no)" \
    "$([ "${h23:-0}" -gt 0 ] && echo yes || echo no)" >> "$Q/rna_summary.tsv"
done
echo

# ==========================================================================
# 5. MERGE  ->  MAG_QC.tsv  + MIMAG classification
# ==========================================================================
echo
echo "######## MERGING QC TABLES ########"
python3 - "$Q" $SCORES <<'PY'
import sys, os, csv, glob, statistics

Q = sys.argv[1]; scores = sys.argv[2:]

def read_tsv(path, key):
    out = {}
    if not path or not os.path.exists(path): return out
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            k = row.get(key) or row.get(key.lower())
            if k: out[k.replace('.fa','')] = row
    return out

memb = read_tsv(os.path.join(Q, 'union_membership.tsv'), 'bin')
cm2  = read_tsv(os.path.join(Q, 'checkm2', 'quality_report.tsv'), 'Name')
rna  = read_tsv(os.path.join(Q, 'rna_summary.tsv'), 'bin')

g = glob.glob(os.path.join(Q, 'gunc', 'GUNC*maxCSS*.tsv'))
gunc = read_tsv(g[0], 'genome') if g else {}
if not gunc:
    print("  !! no GUNC output found — chimerism column will be blank")

def f(d, k, default=''):
    v = (d or {}).get(k, default)
    return v if v not in (None, '') else default

rows = []
for b, m in sorted(memb.items()):
    c = cm2.get(b, {}); r = rna.get(b, {}); u = gunc.get(b, {})
    comp = float(f(c, 'Completeness', 'nan') or 'nan')
    cont = float(f(c, 'Contamination', 'nan') or 'nan')
    agl  = float(f(c, 'Average_Gene_Length', 'nan') or 'nan')
    ntr  = int(f(r, 'n_trna_types', 0) or 0)
    rrna_ok = all(f(r, k) == 'yes' for k in ('has_5S', 'has_16S', 'has_23S'))
    css  = f(u, 'clade_separation_score', '')
    chim = (css != '' and float(css) > 0.45)

    if comp == comp and cont == cont:
        if comp >= 90 and cont < 5 and rrna_ok and ntr >= 18: mimag = 'high'
        elif comp >= 50 and cont < 10:                        mimag = 'medium'
        else:                                                 mimag = 'low'
    else:
        mimag = 'unknown'

    rows.append(dict(
        bin=b, sample=m['sample'], binner=m['binner'],
        n_contigs=m['n_contigs'], bp=m['bp'],
        completeness=f(c,'Completeness'), contamination=f(c,'Contamination'),
        avg_gene_length=f(c,'Average_Gene_Length'),
        coding_density=f(c,'Coding_Density'),
        total_coding_seqs=f(c,'Total_Coding_Sequences'),
        n_trna_types=ntr, has_5S=f(r,'has_5S'), has_16S=f(r,'has_16S'),
        has_23S=f(r,'has_23S'),
        gunc_css=css, gunc_chimeric=('yes' if chim else ('no' if css!='' else '')),
        mimag=mimag,
        **{f'in_s{s}': m.get(f'in_s{s}','') for s in scores}))

cols = list(rows[0].keys()) if rows else []
with open(os.path.join(Q, 'MAG_QC.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
    w.writerows(rows)

# ---------------------------------------------------- THE HEADLINE NUMBER
agls = [float(r['avg_gene_length']) for r in rows if r['avg_gene_length']]
print()
print("#" * 62)
print("#  AVERAGE GENE LENGTH — the reason workspace_2 exists")
print("#" * 62)
if agls:
    med = statistics.median(agls)
    print(f"  workspace_1 baseline (38 MAGs) : 244 aa median")
    print(f"  workspace_2 union    ({len(agls):>2} bins) : {med:.1f} aa median"
          f"   (range {min(agls):.0f}-{max(agls):.0f})")
    print(f"  expected prokaryotic range     : 290-320 aa")
    n_ok = sum(1 for a in agls if a >= 290)
    print(f"  bins at or above 290 aa        : {n_ok}/{len(agls)}"
          f"   (workspace_1: 7/38)")
    print()
    if med >= 290:
        print("  >> CONFIRMED. The --nano-raw + no-Medaka diagnosis was correct;")
        print("     gene calls are now in the expected range and every downstream")
        print("     functional result from workspace_1 was undercounting genes.")
    elif med >= 270:
        print("  >> PARTIAL. Substantially improved but short of 290 aa. Report the")
        print("     improvement honestly and investigate residual error: check")
        print("     per-bin coding density and whether low-coverage bins dominate.")
    else:
        print("  >> NOT CONFIRMED. The hypothesis does not survive. Do NOT present")
        print("     workspace_2 as a correction of workspace_1 on this basis.")
        print("     Alternative causes to test: coverage-dependent consensus")
        print("     failure, or a genuine biological signal (many small-ORF taxa).")
else:
    print("  !! no CheckM2 results — cannot evaluate")

# ------------------------------------------------- MIMAG PER THRESHOLD
print()
print("#" * 62)
print("#  MIMAG YIELD PER DAS_Tool SCORE THRESHOLD")
print("#" * 62)
print(f"  {'sample':<6} {'score':>5} {'high':>5} {'medium':>7} {'HQ+MQ':>6} "
      f"{'low':>5} {'chimeric':>9}")
best = {}
for sample in ('SED', 'LIQ'):
    for s in scores:
        sub = [r for r in rows if r['sample'] == sample and r.get(f'in_s{s}') == 'yes']
        hi = sum(1 for r in sub if r['mimag'] == 'high')
        me = sum(1 for r in sub if r['mimag'] == 'medium')
        lo = sum(1 for r in sub if r['mimag'] == 'low')
        ch = sum(1 for r in sub if r['gunc_chimeric'] == 'yes')
        print(f"  {sample:<6} {s:>5} {hi:>5} {me:>7} {hi+me:>6} {lo:>5} {ch:>9}")
        best.setdefault(sample, []).append((s, hi + me, ch))
    print()

print("  DECISION RULE (fixed at step 07b, before these numbers were seen):")
print("  adopt the threshold maximising HQ+MQ, and prefer the STRICTER threshold")
print("  when the gain is <=1 bin or when it adds chimeric bins.")
for sample, opts in best.items():
    opts_sorted = sorted(opts, key=lambda t: (-t[1], -float(t[0])))
    top = opts_sorted[0]
    strictest_within_1 = sorted([o for o in opts if o[1] >= top[1] - 1],
                                key=lambda t: -float(t[0]))[0]
    print(f"    {sample}: suggested score_threshold = {strictest_within_1[0]} "
          f"({strictest_within_1[1]} HQ+MQ, {strictest_within_1[2]} chimeric)")
print()
print("  This is a SUGGESTION from the rule, not an automatic adoption.")
print("  Inspect the chimeric and low rows before committing in the notebook.")
PY

# ==========================================================================
# 6. VALIDATION
# ==========================================================================
echo
echo "######## VALIDATION GATE ########"
check "union bin set built"        $([ -s "$Q/union_membership.tsv" ] && echo 0 || echo 1)
check "CheckM2 quality report"     $([ -s "$Q/checkm2/quality_report.tsv" ] && echo 0 || echo 1)
check "GUNC output present"        $(ls "$Q/gunc"/GUNC*maxCSS*.tsv >/dev/null 2>&1 && echo 0 || echo 1)
check "rRNA/tRNA summary"          $([ -s "$Q/rna_summary.tsv" ] && echo 0 || echo 1)
check "merged MAG_QC.tsv"          $([ -s "$Q/MAG_QC.tsv" ] && echo 0 || echo 1)

if [ -s "$Q/checkm2/quality_report.tsv" ]; then
  nq=$(( $(wc -l < "$Q/checkm2/quality_report.tsv") - 1 ))
  check "CheckM2 assessed every union bin ($nq/$NBINS)" \
    $([ "$nq" -eq "$NBINS" ] && echo 0 || echo 1)
fi
if [ -s "$Q/MAG_QC.tsv" ]; then
  nagl=$(awk -F'\t' 'NR>1 && $8!=""' "$Q/MAG_QC.tsv" | wc -l)
  check "Average_Gene_Length populated ($nagl bins)" \
    $([ "$nagl" -gt 0 ] && echo 0 || echo 1)
fi

runlog "08_binqc complete"
echo
echo "######## 08 DONE — next: 09_taxonomy.sh (GTDB-Tk R232) ########"
echo "Do not adopt a threshold until we have looked at the chimeric rows."
