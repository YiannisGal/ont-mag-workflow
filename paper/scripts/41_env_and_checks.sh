#!/usr/bin/env bash
# 41_env_and_checks.sh -- namespace probe for the SMETANA environment file, plus four
# short checks. Nothing here runs long; nothing here writes outside $WS2/41_env_checks.
#
# A. SBML namespace probe   -- what ID space do gapseq models actually use?
# B. in-situ 10 genome list -- derived from cross_abundance.tsv, not hand-typed
# C. cross_abundance re-read at full precision + the methanogen/Methanosarcina census
# D. Sedimentibacter and Chakrabartyella SED x LIQ ANI
# E. anaerobic aromatic degradation SCREEN (not a result -- see the banner it prints)
set -uo pipefail

WS2="${WS2:-${WS2}}"
OUT="$WS2/41_env_checks"
export TMPDIR="${TMPDIR41:-${SCRATCH_ROOT}/tmp41}"
mkdir -p "$OUT" "$TMPDIR" || { echo "FATAL: cannot create $OUT"; exit 1; }
CULT="${CULT:-$HOME/Documents/UniqueELT_Solid/Manuscript_Analysis}"

say(){ printf '\n=== %s ===\n' "$*"; }
ok(){  printf '  [OK]   %s\n' "$*"; }
no(){  printf '  [FAIL] %s\n' "$*"; }
note(){ printf '  ...    %s\n' "$*"; }

echo "41_env_and_checks  $(date -Is)"
echo "WS2  = $WS2"
echo "OUT  = $OUT"

# ---------------------------------------------------------------- A. namespace
say "A. SBML namespace probe"
# Find a COMPLETED gapseq model: a *.xml that is not a *-draft.xml.
mapfile -t MODELS < <(find "$CULT" -name '*.xml' ! -name '*-draft.xml' -size +20k 2>/dev/null | sort)
if [ "${#MODELS[@]}" -eq 0 ]; then
  no "no completed gapseq model found under $CULT"
  no "cannot determine the namespace -- the environment file MUST NOT be written yet"
else
  M="${MODELS[0]}"
  ok "reading $M  ($(du -h "$M" | cut -f1))"
  python3 - "$M" "$OUT" <<'PY'
import re, sys, collections
path, out = sys.argv[1], sys.argv[2]
txt = open(path, encoding='utf-8', errors='replace').read()

# exchange reactions, however they are spelled
rx = sorted(set(re.findall(r'id="(R_)?(EX_[A-Za-z0-9_\-\[\]]+)"', txt)))
rx = sorted({r[1] for r in rx})
print(f"  ...    exchange reactions found: {len(rx)}")
for r in rx[:8]:
    print(f"           {r}")

pat = collections.Counter()
for r in rx:
    body = r[3:]
    if   re.match(r'^cpd\d{5}', body):        pat['ModelSEED  cpd#####'] += 1
    elif re.match(r'^C\d{5}',   body):        pat['KEGG       C#####']   += 1
    elif re.match(r'^M_',       body):        pat['SBML       M_ prefix']+= 1
    else:                                     pat['other/BiGG-like']     += 1
print("  ...    namespace tally:")
for k, v in pat.most_common():
    print(f"           {v:>6}  {k}")

if not rx:
    print("  [FAIL] no exchange reactions -- wrong file or unexpected SBML flavour")
    sys.exit(0)
top, n = pat.most_common(1)[0]
frac = n / len(rx)
if frac < 0.90:
    print(f"  [FAIL] namespace is mixed ({frac:.0%} {top}) -- do not write the environment file")
else:
    print(f"  [OK]   namespace is {top} at {frac:.0%}")

# metabolite id -> name, so compounds can be looked up BY NAME rather than by a
# hardcoded ID table. This is the whole point of the probe.
names = dict(re.findall(r'<species[^>]*\bid="([^"]+)"[^>]*\bname="([^"]*)"', txt))
if not names:
    names = dict(re.findall(r'<species[^>]*\bname="([^"]*)"[^>]*\bid="([^"]+)"', txt))
    names = {v: k for k, v in names.items()}
print(f"  ...    species with names: {len(names)}")

with open(f"{out}/sbml_exchange_metabolites.tsv", "w") as fh:
    fh.write("exchange_rxn\tmetabolite_id\tname\n")
    for r in rx:
        body = r[3:]
        cands = [k for k in names if body.split('_e')[0] in k] or [body]
        mid = cands[0]
        fh.write(f"{r}\t{mid}\t{names.get(mid,'')}\n")
print(f"  [OK]   wrote {out}/sbml_exchange_metabolites.tsv")

WANT = {
 "ammonium":["ammonium","nh4","ammonia"], "phosphate":["phosphate","orthophosphate"],
 "sulfate":["sulfate"], "sulfide":["h2s","hydrogen sulfide","sulfide"],
 "bicarbonate":["bicarbonate","hco3","co2","carbon dioxide"],
 "water":["h2o","water"], "proton":["h+","proton"],
 "potassium":["k+","potassium"], "sodium":["na+","sodium"],
 "magnesium":["mg2+","magnesium"], "calcium":["ca2+","calcium"],
 "chloride":["cl-","chloride"], "iron2":["fe2+","ferrous","fe2"],
 "cellulose/glucan":["cellobiose","cellulose","glucan"],
 "xylan/xylose":["xylose","xylan","d-xylose"],
 "pectin/galacturonate":["d-galacturonate","galacturonate","pectin"],
 "acetate":["acetate"], "h2":["h2","hydrogen"],
 "EXCLUDED-methanol":["methanol"], "EXCLUDED-TMA":["trimethylamine"],
 "EXCLUDED-choline":["choline"], "EXCLUDED-B12":["cobalamin","b12"],
}
low = {k: (v or "").strip().lower() for k, v in names.items()}
print("  ...    can the intended environment be expressed in this namespace?")
print("         EXACT = name equals the key. SUBSTR = fell back to containment,")
print("         which is ambiguous ('h2' is inside 'h2s' and 'h2o') and needs a human.")
miss, ambig = [], []
for label, keys in WANT.items():
    exact  = sorted({k for k, v in low.items() if v in keys})
    substr = sorted({k for k, v in low.items() if any(kk in v for kk in keys)}) if not exact else []
    hit, how = (exact, "EXACT") if exact else (substr, "SUBSTR")
    req = not label.startswith("EXCLUDED")
    if not hit:
        if req: miss.append(label)
        print(f"           [--    ] {label:<24}   0 candidate(s)")
        continue
    if how == "SUBSTR" and len(hit) > 1 and req:
        ambig.append(f"{label} ({len(hit)}: {', '.join(hit[:4])})")
    flag = "OK    " if how == "EXACT" and len(hit) == 1 else \
           ("OK*   " if how == "EXACT" else "AMBIG ")
    print(f"           [{flag}] {label:<24} {len(hit):>3} {how:<6} e.g. {hit[0]}"
          + (f"  <-- {low[hit[0]]!r}" if how == "SUBSTR" else ""))
bad = False
if miss:
    bad = True
    print(f"  [FAIL] {len(miss)} required compound(s) not found by name: {', '.join(miss)}")
if ambig:
    bad = True
    print(f"  [FAIL] {len(ambig)} required compound(s) matched only ambiguously:")
    for a in ambig: print(f"           {a}")
    print("         A wrong ID here does not error -- SMETANA runs and returns no")
    print("         cross-feeding. Pick these by hand from sbml_exchange_metabolites.tsv.")
if bad:
    print("  [FAIL] DO NOT write the environment file until the above are resolved.")
else:
    print("  [OK]   every required compound resolves to exactly one metabolite")
PY
fi

# ------------------------------------------------------- B. the in-situ 10 list
say "B. in-situ genome set, derived not typed"
CA=""
for c in "$WS2"/26_crossmap/cross_abundance.tsv "$WS2"/*/cross_abundance.tsv; do
  [ -f "$c" ] && { CA="$c"; break; }
done
if [ -z "$CA" ]; then
  no "cross_abundance.tsv not found under $WS2"
else
  ok "using $CA"
  python3 - "$CA" "$OUT" <<'PY'
import sys, csv
path, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(path), delimiter='\t'))
print(f"  ...    {len(rows)} rows, columns: {', '.join(rows[0].keys())}")
key = next((k for k in rows[0] if k.strip().lower() == 'sed_breadth'), None)
if not key:
    print("  [FAIL] no sed_breadth column"); sys.exit(0)
def f(r):
    try: return float(r[key])
    except: return -1.0
rows.sort(key=f, reverse=True)

# find the largest gap in the top half -- report it, do not assume 0.85
vals = [f(r) for r in rows]
gaps = sorted(((vals[i]-vals[i+1], i) for i in range(len(vals)-1)), reverse=True)
print("  ...    three largest gaps in the sorted breadth distribution:")
for g, i in gaps[:3]:
    print(f"           {g:.4f}  between rank {i+1} ({vals[i]:.4f}) and rank {i+2} ({vals[i+1]:.4f})")
g, i = gaps[0]
n = i + 1
print(f"  ...    largest gap puts {n} genomes above the cut")
if n != 10:
    print(f"  [FAIL] expected 10 in-situ genomes, the data's own largest gap gives {n}")
    print("         Do NOT proceed on the assumption of 10 -- reconcile with step 26 first.")
else:
    print("  [OK]   the largest gap in the data gives exactly 10 -- matches step 26")

idc = next((k for k in rows[0] if k.strip().lower() in ('mag_id','working_name','genome')), list(rows[0])[0])
# The filename carries the ACTUAL count. Only n==10 gets the canonical name that
# downstream steps look for -- a file called in_situ_10.txt must never hold 11.
fname = "in_situ_10.txt" if n == 10 else f"in_situ_{n}_UNRECONCILED.txt"
with open(f"{out}/{fname}","w") as fh:
    for r in rows[:n]:
        fh.write(r[idc] + "\n")
if n == 10:
    print(f"  [OK]   wrote {out}/{fname}")
else:
    print(f"  [FAIL] wrote {out}/{fname} -- NOT the canonical name.")
    print("         Nothing downstream will pick this up, which is deliberate.")
with open(f"{out}/cross_abundance_full_precision.tsv","w") as fh:
    fh.write("\t".join(rows[0].keys())+"\n")
    for r in rows: fh.write("\t".join(str(r[k]) for k in rows[0].keys())+"\n")
print(f"  [OK]   wrote {out}/cross_abundance_full_precision.tsv  (sorted, no column -t damage)")
PY
fi

# ---------------------------------------------- C. methanogen / genus census
say "C. methanogen census and the Methanosarcina count"
TAX=""
for t in "$WS2"/MAG_TAXONOMY.tsv "$WS2"/*/MAG_TAXONOMY.tsv; do [ -f "$t" ] && { TAX="$t"; break; }; done
if [ -z "$TAX" ]; then no "MAG_TAXONOMY.tsv not found"; else
  ok "using $TAX"
  python3 - "$TAX" "$OUT" <<'PY'
import sys, csv, collections
path, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(path), delimiter='\t'))
gk = next((k for k in rows[0] if k.strip().lower()=='genus'), None)
nk = next((k for k in rows[0] if k.strip().lower()=='working_name'), list(rows[0])[0])
dk = next((k for k in rows[0] if k.strip().lower() in ('domain','kingdom')), None)
if not gk: print("  [FAIL] no genus column"); sys.exit(0)
meth = [r for r in rows if 'methano' in (r[gk] or '').lower()]
print(f"  ...    {len(rows)} genomes; {len(meth)} with a Methano* genus")
for r in sorted(meth, key=lambda r: r[gk]):
    print(f"           {r[nk]:<28} {r[gk]}")
c = collections.Counter(r[gk] for r in meth)
print("  ...    by genus:")
for k,v in sorted(c.items()): print(f"           {v}  {k}")
ms = c.get('Methanosarcina', 0)
if ms == 2:   print("  [OK]   Methanosarcina = 2 -- cross_abundance was right, the 'three' was a drafting error")
elif ms == 3: print("  [OK]   Methanosarcina = 3 -- MEDIUM_DEFINITION was right, re-check cross_abundance labels")
else:         print(f"  [FAIL] Methanosarcina = {ms} -- neither document matches; resolve by hand")
if dk:
    arch = [r for r in rows if 'archae' in (r[dk] or '').lower()]
    print(f"  ...    archaea by domain column: {len(arch)}")
    for r in arch: print(f"           {r[nk]:<28} {r[gk]}")
PY
fi

# ------------------------------------------------------------------- D. ANI
say "D. Sedimentibacter / Chakrabartyella  SED x LIQ"
if ! command -v skani >/dev/null 2>&1; then
  no "skani not on PATH -- skipping"
else
  mapfile -t FNA < <(find "$WS2"/07_mags -maxdepth 2 -name '*.fna' \
      ! -name 'all_mags_contigs.fna' ! -name '*all_mags*' 2>/dev/null | sort)
  if [ "${#FNA[@]}" -lt 2 ]; then no "fewer than 2 per-MAG FASTAs found"; else
    ok "${#FNA[@]} per-MAG FASTAs (aggregate files excluded by name)"
    printf '%s\n' "${FNA[@]}" > "$TMPDIR/all.txt"
    timeout 900 skani dist --ql "$TMPDIR/all.txt" --rl "$TMPDIR/all.txt" \
      -t 4 -o "$OUT/mag_x_mag_ani.tsv" 2>"$OUT/skani.err"
    rc=$?
    if [ $rc -ne 0 ]; then no "skani exit $rc -- see $OUT/skani.err"; else
      n=$(( $(wc -l < "$OUT/mag_x_mag_ani.tsv") - 1 ))
      ok "$n pairs above skani's screening threshold -> $OUT/mag_x_mag_ani.tsv"
      echo "  ...    pairs at >=95% ANI, excluding self-hits:"
      awk -F'\t' 'NR>1 && $3>=95 && $1!=$2 {printf "           %6.2f%%  %s  vs  %s\n",$3,$1,$2}' \
        "$OUT/mag_x_mag_ani.tsv" | sed 's#/[^ ]*/##g' | sort -rn | head -20
      echo "  ...    (the Chakrabartyella pair should appear at ~99.75% -- that is the control)"
    fi
  fi
fi

# --------------------------------------------------- E. aromatic degradation
say "E. anaerobic aromatic degradation -- SCREEN ONLY"
cat <<'BANNER'
  !! This is a product-string screen, which is exactly what step 16 was rejected for.
  !! A positive here justifies a proper KOfam pass at curated thresholds.
  !! A negative here is WEAK evidence and closes nothing.
BANNER
ANN=$(find "$WS2"/09_annotation -maxdepth 2 -name '*.tsv' 2>/dev/null | head -1)
if [ -z "$ANN" ]; then no "no Bakta .tsv under 09_annotation -- cannot screen"; else
  ok "screening $(find "$WS2"/09_annotation -maxdepth 2 -name '*.tsv' | wc -l) annotation tables"
  grep -Hi -E 'benzoyl-CoA|benzoate|4-hydroxybenzo|phenylphosphate|phenol|catechol|resorcinol|cyclohex-1,5-diene' \
    "$WS2"/09_annotation/*/*.tsv 2>/dev/null \
    | awk -F'\t' '{split($1,a,"/"); print a[length(a)-1]}' | sort | uniq -c | sort -rn \
    > "$OUT/aromatic_screen_by_genome.txt"
  if [ -s "$OUT/aromatic_screen_by_genome.txt" ]; then
    tot=$(awk '{s+=$1} END{print s+0}' "$OUT/aromatic_screen_by_genome.txt")
    ng=$(wc -l < "$OUT/aromatic_screen_by_genome.txt")
    ok "$tot hits across $ng genomes -> $OUT/aromatic_screen_by_genome.txt"
    head -12 "$OUT/aromatic_screen_by_genome.txt" | sed 's/^/           /'
    if [ -f "$OUT/in_situ_10.txt" ]; then
      echo "  ...    restricted to the in-situ set:"
      grep -Ff "$OUT/in_situ_10.txt" "$OUT/aromatic_screen_by_genome.txt" | sed 's/^/           /' \
        || echo "           none of the in-situ 10 carry a hit"
    fi
    echo "  ...    VERDICT: hits present -> a curated KOfam pass on bcr/bam is worth running."
  else
    note "no product-string hits -- weak evidence only, does not close the question"
  fi
fi

say "DONE"
echo "outputs in $OUT"
echo "STEP_STATUS: COMPLETE  $(date -Is)"
