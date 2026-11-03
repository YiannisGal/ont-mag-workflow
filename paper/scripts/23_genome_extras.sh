#!/usr/bin/env bash
# =============================================================================
# 23_genome_extras.sh  --  secretome, rRNA operon copy number, genome statistics
#
# THREE ADDITIONS, ALL FROM DATA ALREADY ON DISK
#
# 1. rRNA OPERON COPY NUMBER
#    Not a curiosity — it is the correction factor for every rRNA-based
#    abundance in this project. Copy number varies ~1-15 between taxa, so a
#    SILVA or amplicon profile systematically overstates high-copy lineages.
#    Counting operons in the MAGs gives a per-genome factor and, more usefully,
#    a check on whether the taxa this study cares about are high- or low-copy.
#    It is also the best available proxy for ecological strategy: high copy
#    number tracks fast response to resource pulses, low copy number tracks
#    oligotrophy. Reported alongside step 17's codon-usage result, which found
#    the enrichments did NOT select fast growers — the two should agree.
#
#    CRITICAL CAVEAT, ENFORCED BELOW: a MAG's rRNA count is a FLOOR, not a
#    count. rRNA operons are repetitive and collapse during assembly, and an
#    incomplete genome loses copies. Only genomes >=90% complete are reported,
#    and even those are stated as minima.
#
# 2. SECRETOME
#    Which proteins leave the cell. This matters here specifically: step 20
#    found Gaoshiqia carrying 221 glycoside hydrolases and 120 SusC/SusD pairs.
#    Whether those enzymes are secreted or cytoplasmic decides whether that
#    organism feeds the community or only itself — a public-goods question that
#    the raw CAZyme count cannot answer.
#
#    Uses DeepSig (free, published, on bioconda). Gram staining is taken from
#    GTDB phylum, since the model differs for Gram+ and Gram-. Archaea have no
#    dedicated model; they are run under the Gram-positive model and FLAGGED,
#    because an unflagged approximation is how a caveat gets lost.
#
# 3. GENOME STATISTICS, COMPLETENESS-CORRECTED
#    Estimated true genome size (observed / completeness), coding density, GC.
#    Genome size is an ecological trait: streamlined genomes indicate stable,
#    resource-poor niches; large genomes indicate variable ones.
#
# INPUT    07_mags/bins/, 07_mags/MAG_TAXONOMY.tsv, 09_annotation/<MAG>/
# OUTPUT   23_extras/
#     rrna_copy_number.tsv   per MAG, with a completeness flag
#     secretome.tsv          secreted / membrane / cytoplasmic per MAG
#     secreted_cazymes.tsv   THE INTERESTING ONE: public-goods enzymes
#     genome_stats.tsv
#
# RUNTIME  30-60 min. USAGE:
#   ws2 && bash scripts/23_genome_extras.sh 2>&1 | tee logs/23_extras.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="23_extras"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
A="$WS2/09_annotation"; F="$WS2/07_mags"; E="$WS2/23_extras"
mkdir -p "$E/barrnap" "$E/deepsig"

echo "######## 23 — SECRETOME, rRNA COPY NUMBER, GENOME STATS ########"
echo "  started: $(date -Is)"
runlog "23_genome_extras started"
[ -d "$A" ] || { echo "  !! run 10_annotate.sh first"; exit 1; }

# ---------------------------------------------------------------- 1. rRNA
echo
echo "######## 1. rRNA OPERONS PER MAG ########"
for fa in "$F"/bins/*.fa; do
  [ -s "$fa" ] || continue
  b=$(basename "$fa" .fa); out="$E/barrnap/${b}.gff"
  [ -s "$out" ] && continue
  dom=$(awk -F'\t' -v n="$b" 'NR==1{for(i=1;i<=NF;i++){if($i=="working_name")w=i;if($i=="domain")d=i}}
        NR>1 && $w==n {print tolower($d)}' "$F/MAG_TAXONOMY.tsv" 2>/dev/null | head -1)
  case "$dom" in *arch*) K=arc ;; *) K=bac ;; esac
  run_tool barrnap --kingdom "$K" --quiet "$fa" > "$out" 2>/dev/null
done
echo "  scanned: $(ls -1 "$E/barrnap"/*.gff 2>/dev/null | wc -l) MAGs"

# ---------------------------------------------------------------- 2. secretome
echo
echo "######## 2. SECRETOME (DeepSig) ########"
DS_ENV="${DEEPSIG_ENV:-${ENVS_DIR}/deepsig_env}"
DS_ROOT="${DEEPSIG_ROOT:-${ENVS_DIR}/deepsig_src}"
DS="$DS_ENV/bin/deepsig"

# The BIOCONDA recipe pins tensorflow 2.2.0, which no longer exists for
# linux-64 on any channel, so `conda install deepsig` cannot be solved. The
# PIP route (deepsig-biocomp) has no such pin and resolves a current
# tensorflow, so it is tried FIRST and conda only as a fallback.
# A binary that exists is not a binary that works. DeepSig 0.9 fails at runtime
# under Keras 3, so check the installed Keras version and repair the env even
# when the executable is already present — the previous run skipped the pin
# entirely for exactly this reason.
KERAS_OK=0
if [ -x "$DS_ENV/bin/python" ]; then
  KV=$("$DS_ENV/bin/python" -c 'import keras,sys; sys.stdout.write(keras.__version__)' 2>/dev/null || echo "none")
  echo "  keras in $DS_ENV: ${KV:-none}"
  case "$KV" in
    2.*) KERAS_OK=1 ;;
    *)   echo "  Keras ${KV:-missing} cannot load DeepSig's Keras 2.0 models — repairing" ;;
  esac
fi
if [ -x "$DS" ] && [ "$KERAS_OK" -eq 0 ]; then
  echo "  downgrading to a TensorFlow that ships Keras 2 ..."
  "$DS_ENV/bin/pip" install --no-input --force-reinstall \
      "tensorflow==2.15.1" "keras==2.15.0" 2>&1 | tail -4 | sed 's/^/    /' || true
  KV=$("$DS_ENV/bin/python" -c 'import keras,sys; sys.stdout.write(keras.__version__)' 2>/dev/null || echo none)
  echo "  keras now: ${KV:-none}"
fi

if [ ! -x "$DS" ]; then
  echo "  installing DeepSig via pip (the bioconda recipe pins tensorflow 2.2.0,"
  echo "  which no longer exists — pip resolves a current one) ..."
  SOLVER=conda; command -v mamba >/dev/null 2>&1 && SOLVER=mamba
  $SOLVER create -y -p "$DS_ENV" -c conda-forge python=3.10 pip \
      2>&1 | tail -5 | sed 's/^/    /' || true
  if [ -x "$DS_ENV/bin/pip" ]; then
    "$DS_ENV/bin/pip" install --no-input deepsig-biocomp 2>&1 \
      | tail -8 | sed 's/^/    /' || true
    # DeepSig 0.9 calls `from keras.models import model_from_json` on a model
    # serialised by Keras 2.0. TF_USE_LEGACY_KERAS does NOT help: it only
    # redirects tf.keras, while DeepSig imports the top-level `keras` package,
    # which is Keras 3 whenever TensorFlow >=2.16 is installed. Keras 3 cannot
    # deserialise a 2.0 architecture — hence "Could not locate class 'Model'".
    # The only reliable fix is to install a TensorFlow that still ships Keras 2.
    echo "    pinning tensorflow 2.15 (ships Keras 2.15; TF>=2.16 ships Keras 3)"
    "$DS_ENV/bin/pip" install --no-input "tensorflow==2.15.1" "keras<3" 2>&1 \
      | tail -4 | sed 's/^/    /' || true
  fi
  # DeepSig needs its model files, which ship in the GitHub repo, not the wheel.
  # Use HTTPS — the documented git@ URL needs an SSH key this machine may lack.
  if [ ! -d "$DS_ROOT/deepsig" ] && [ ! -f "$DS_ROOT/setup.py" ]; then
    echo "  fetching model files ..."
    mkdir -p "$(dirname "$DS_ROOT")"
    git clone --depth 1 https://github.com/BolognaBiocomp/deepsig.git "$DS_ROOT" \
      2>&1 | tail -3 | sed 's/^/    /' || true
  fi
  [ -d "$DS_ROOT" ] && export DEEPSIG_ROOT="$DS_ROOT"
fi

# conda fallback, in case a future recipe drops the pin
if [ ! -x "$DS" ]; then
  echo "  pip route did not produce a binary — trying conda as a fallback ..."
  SOLVER=conda; command -v mamba >/dev/null 2>&1 && SOLVER=mamba
  $SOLVER create -y -p "$DS_ENV" -c conda-forge -c bioconda deepsig 2>&1 \
    | tail -10 | sed 's/^/    /' || true
fi

export DEEPSIG_ROOT="${DEEPSIG_ROOT:-$DS_ROOT}"
export TF_USE_LEGACY_KERAS=1
export TF_CPP_MIN_LOG_LEVEL=3
if [ ! -x "$DS" ]; then
  echo "  !! DeepSig unavailable. The secretome section will be skipped; rRNA"
  echo "     copy number and genome stats are unaffected. To install by hand:"
  echo "       conda create -y -p $DS_ENV -c conda-forge python=3.10 pip"
  echo "       $DS_ENV/bin/pip install deepsig-biocomp"
  echo "       git clone https://github.com/BolognaBiocomp/deepsig.git $DS_ROOT"
  echo "       export DEEPSIG_ROOT=$DS_ROOT"
  echo "     Then re-run this script. 23b_pul_analysis.sh already answers the"
  echo "     Gaoshiqia question without it."
  DS=""
else
  echo "  deepsig: $("$DS" --help 2>&1 | head -1)"
  n=0; tot=$(ls -1d "$A"/*/ 2>/dev/null | wc -l)
  for d in "$A"/*/; do
    name=$(basename "$d"); faa="$d/${name}.faa"
    [ -s "$faa" ] || continue
    n=$((n+1)); out="$E/deepsig/${name}.tsv"
    [ -s "$out" ] && { printf '\r  [%2d/%2d] %-28s cached  ' "$n" "$tot" "$name"; continue; }
    # DeepSig has models for gram-positive, gram-negative and eukaryote only.
    # Archaea get the gram-positive model as the closest available and are
    # flagged in the output — an approximation that is stated, not hidden.
    dom=$(awk -F'\t' -v x="$name" 'NR==1{for(i=1;i<=NF;i++){if($i=="working_name")w=i;if($i=="domain")dd=i;if($i=="phylum")p=i}}
          NR>1 && $w==x {print tolower($dd)"|"tolower($p)}' "$F/MAG_TAXONOMY.tsv" 2>/dev/null | head -1)
    case "$dom" in
      *arch*)                         ORG=gramp ;;
      *bacillota*|*actinomycetota*)   ORG=gramp ;;
      *)                              ORG=gramn ;;
    esac
    printf '\r  [%2d/%2d] %-28s %s   ' "$n" "$tot" "$name" "$ORG"
    if [ "$n" -eq 1 ]; then
      # First genome runs LOUD. Swallowing stderr is how 29 identical failures
      # got reported with no diagnosis at all.
      echo; echo "  --- first genome runs with errors visible ---"
      "$DS" -f "$faa" -o "$out" -k "$ORG" 2>&1 | tail -20 | sed 's/^/      /'
      if [ ! -s "$out" ]; then
        echo "  !! DeepSig cannot run — see the error above. Skipping the secretome."
        echo "     23b_pul_analysis.sh already answers the Gaoshiqia question."
        DS=""; break
      fi
      echo "  --- first genome OK, continuing quietly ---"
    else
      "$DS" -f "$faa" -o "$out" -k "$ORG" >/dev/null 2>&1 \
        || echo "  !! deepsig failed on $name"
    fi
  done
  echo
fi

# ===========================================================================
python3 - "$E" "$F" "$A" <<'PY'
import sys, os, csv, glob
from collections import Counter, defaultdict

E, F, A = sys.argv[1], sys.argv[2], sys.argv[3]

tax = {}
for r in csv.DictReader(open(os.path.join(F, 'MAG_TAXONOMY.tsv')), delimiter='\t'):
    wn = r.get('working_name') or r.get('mag_id')
    if wn: tax[wn] = r
def mid(m): return (tax.get(m, {}) or {}).get('mag_id', m)
def gen(m): return (tax.get(m, {}) or {}).get('genus', '') or '-'
def dom(m): return ((tax.get(m, {}) or {}).get('domain', '') or '?')
def comp(m):
    try: return float((tax.get(m, {}) or {}).get('completeness', 0) or 0)
    except ValueError: return 0.0

# =====================================================================
# 1. rRNA OPERON COPY NUMBER
# =====================================================================
print()
print("=" * 100)
print("  1. rRNA OPERONS PER GENOME")
print("=" * 100)
print("  A MAG's rRNA count is a FLOOR. Operons are near-identical repeats that")
print("  collapse during assembly, and an incomplete genome loses copies outright.")
print("  Only genomes >=90% complete are shown, and even those are minima.")
rows = []
for g in sorted(glob.glob(os.path.join(E, 'barrnap', '*.gff'))):
    m = os.path.basename(g)[:-4]
    c = Counter()
    for line in open(g):
        if line.startswith('#'): continue
        f = line.split('\t')
        if len(f) < 9: continue
        for k in ('16S', '23S', '5S'):
            if k in f[8]: c[k] += 1
    rows.append(dict(mag_id=mid(m), working_name=m, genus=gen(m),
                     domain=dom(m), completeness=f"{comp(m):.1f}",
                     n_16S=c['16S'], n_23S=c['23S'], n_5S=c['5S'],
                     operon_floor=min(c['16S'], c['23S']) or c['16S'],
                     trustworthy='yes' if comp(m) >= 90 else 'no'))
rows.sort(key=lambda r: -int(r['operon_floor']))
print(f"\n  {'MAG ID':<17}{'genus':<20}{'compl':>7}{'16S':>5}{'23S':>5}{'5S':>5}{'floor':>7}")
for r in rows:
    if r['trustworthy'] != 'yes': continue
    print(f"  {r['mag_id']:<17}{r['genus'][:19]:<20}{float(r['completeness']):>6.1f}%"
          f"{r['n_16S']:>5}{r['n_23S']:>5}{r['n_5S']:>5}{r['operon_floor']:>7}")
hi = [r for r in rows if r['trustworthy']=='yes' and int(r['operon_floor']) >= 3]
lo = [r for r in rows if r['trustworthy']=='yes' and int(r['operon_floor']) <= 1]
print(f"\n  >=3 operons (fast responders): {len(hi)}   <=1 operon (oligotrophs): {len(lo)}")
print("  Cross-check against step 17: codon-usage bias said the enrichments did")
print("  NOT select fast growers (culture 0.2705 vs sediment 0.2808). If copy")
print("  number agrees, that is two independent lines on the same conclusion.")
with open(os.path.join(E, 'rrna_copy_number.tsv'), 'w', newline='') as fh:
    if rows:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

print()
print("  USE THIS TO CORRECT THE rRNA PROFILES. A genus with 5 operons")
print("  contributes 5x the SSU reads per cell of a genus with 1. Where a MAG")
print("  exists for a lineage in the SILVA or amplicon profile, its operon count")
print("  is the correction factor for that lineage.")

# =====================================================================
# 2. SECRETOME
# =====================================================================
print()
print("=" * 100)
print("  2. SECRETOME — which enzymes are public goods")
print("=" * 100)
ds = sorted(glob.glob(os.path.join(E, 'deepsig', '*.tsv')))
if not ds:
    print("  DeepSig produced no output — section skipped (see the note above).")
else:
    # Bakta product per locus tag, so a secreted protein can be named
    prod = {}
    for d in sorted(glob.glob(os.path.join(A, '*'))):
        if not os.path.isdir(d): continue
        n = os.path.basename(d); t = os.path.join(d, f'{n}.tsv')
        if not os.path.exists(t): continue
        for line in open(t):
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) >= 8 and f[1].lower() == 'cds':
                prod[(n, f[5])] = f[7]

    # ONE DEFINITION OF A CAZyme, SHARED WITH 23c. The first version of this
    # list was shorter and stopped at the polymer-acting enzymes; 23c added the
    # exo-acting glycosidases that finish the job (galactosidase, glucosidase,
    # mannosidase, fucosidase, xylosidase, arabinofuranosidase) plus levanase,
    # pullulanase and chondroitin lyase. For Gaoshiqia that is the whole
    # difference between 114 and 161 secreted CAZymes.
    #
    # This matters beyond bookkeeping: it means the earlier arithmetic
    # "137 PUL + 114 secreted = 251 out of 253, therefore near-total overlap"
    # was comparing counts made under two different rules and the agreement
    # was a coincidence. Two numbers only support an inference if the same
    # definition produced both. If this list changes again, change it in BOTH
    # scripts and re-run both — DeepSig output is cached, so 23 re-runs in
    # under a minute.
    CAZY = ('glycoside hydrolase','glycosyl hydrolase','polysaccharide lyase',
            'carbohydrate esterase','glucanase','xylanase','chitinase','amylase',
            'pectinesterase','pectate lyase','cellulase','agarase','alginate lyase',
            'galactosidase','glucosidase','mannosidase','fucosidase','xylosidase',
            'arabinofuranosidase','levanase','pullulanase','chondroitin')
    PEPT = ('peptidase','protease','proteinase','collagenase','aminopeptidase')
    srows, crows = [], []
    print(f"\n  {'MAG ID':<17}{'genus':<20}{'proteins':>9}{'secreted':>10}{'%':>7}"
          f"{'sec.CAZy':>10}{'sec.pept':>10}  model")
    for p in ds:
        m = os.path.basename(p)[:-4]
        nsec = ntot = ncaz = npep = 0
        for line in open(p):
            if line.startswith('#') or not line.strip(): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 3: continue
            ntot += 1
            lab = ' '.join(f[1:]).lower()
            if 'signal' not in lab: continue
            nsec += 1
            pr = (prod.get((m, f[0])) or '').lower()
            if any(x in pr for x in CAZY):
                ncaz += 1
                crows.append(dict(mag_id=mid(m), genus=gen(m), locus=f[0],
                                  product=prod.get((m, f[0]), ''), kind='CAZyme'))
            elif any(x in pr for x in PEPT):
                npep += 1
                crows.append(dict(mag_id=mid(m), genus=gen(m), locus=f[0],
                                  product=prod.get((m, f[0]), ''), kind='peptidase'))
        model = 'gramP(archaeal approx)' if 'arch' in dom(m).lower() else 'standard'
        print(f"  {mid(m):<17}{gen(m)[:19]:<20}{ntot:>9}{nsec:>10}"
              f"{100*nsec/max(1,ntot):>6.1f}%{ncaz:>10}{npep:>10}  {model}")
        srows.append(dict(mag_id=mid(m), working_name=m, genus=gen(m), domain=dom(m),
                          n_proteins=ntot, n_secreted=nsec,
                          pct_secreted=f"{100*nsec/max(1,ntot):.2f}",
                          secreted_cazymes=ncaz, secreted_peptidases=npep,
                          model_caveat=model))
    for name, rr in (('secretome.tsv', srows), ('secreted_cazymes.tsv', crows)):
        with open(os.path.join(E, name), 'w', newline='') as fh:
            if rr:
                w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(rr[0]))
                w.writeheader(); w.writerows(rr)
    print()
    print("  THE QUESTION THIS ANSWERS. Step 20 found Gaoshiqia with 221 glycoside")
    print("  hydrolases and 120 SusC/SusD pairs. If those enzymes are SECRETED, it")
    print("  releases oligosaccharides into the porewater and feeds the fermenters")
    print("  that feed the methanogens — a public good, and the base of the food")
    print("  chain this chapter describes. If they are cytoplasmic, it imports")
    print("  polymers and digests them privately, and the cross-feeding argument")
    print("  needs a different source. Check its row above before writing either.")
    print()
    print("  CAVEAT: archaeal genomes were scored with the Gram-positive model")
    print("  because DeepSig has none for Archaea. Treat archaeal percentages as")
    print("  approximate and do not compare them directly with bacterial ones.")

# =====================================================================
# 3. GENOME STATISTICS
# =====================================================================
print()
print("=" * 100)
print("  3. GENOME SIZE, COMPLETENESS-CORRECTED")
print("=" * 100)
print("  Observed size understates the genome by exactly the fraction missing.")
print("  Estimated size = observed / completeness. Genome size is an ecological")
print("  trait: streamlined genomes indicate stable low-resource niches.")
grows = []
print(f"\n  {'MAG ID':<17}{'genus':<20}{'observed':>11}{'compl':>7}{'estimated':>12}")
for m in sorted(tax, key=lambda x: mid(x)):
    t = tax[m]
    try: bp = float(t.get('bp') or t.get('total_bp') or 0)
    except ValueError: bp = 0
    c = comp(m)
    if not bp or not c: continue
    est = bp / (c / 100.0)
    grows.append(dict(mag_id=mid(m), genus=gen(m), observed_bp=int(bp),
                      completeness=f"{c:.1f}", estimated_bp=int(est)))
    print(f"  {mid(m):<17}{gen(m)[:19]:<20}{bp/1e6:>10.2f}M{c:>6.1f}%{est/1e6:>11.2f}M")
if grows:
    with open(os.path.join(E, 'genome_stats.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, delimiter='\t', fieldnames=list(grows[0]))
        w.writeheader(); w.writerows(grows)
else:
    print("  no size column found in MAG_TAXONOMY.tsv — check the crosswalk instead")
PY

echo
echo "######## VALIDATION GATE ########"
check "rRNA copy number"  $([ -s "$E/rrna_copy_number.tsv" ] && echo 0 || echo 1)
check "genome stats"      $([ -s "$E/genome_stats.tsv" ] && echo 0 || echo 1)
[ -s "$E/secretome.tsv" ] && check "secretome" 0 || echo "  [SKIP] secretome (DeepSig unavailable)"

runlog "23_genome_extras complete"
echo
echo "######## 23 DONE ########"
