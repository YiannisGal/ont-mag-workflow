#!/usr/bin/env bash
# =============================================================================
# 40_integration.sh  --  four cheap jobs that need the idle 23 cores
#
# gapseq `find` is single-threaded and using ~3% of this machine. Everything
# here runs alongside it without competing: no heavy I/O, no shared temp, and
# nothing written outside $WS2.
#
#   A. ANI: 29 MAGs x the culturomics isolates      (closes step 12)
#   B. the 30 free-virus contigs sitting inside MAGs (from 39b's topology split)
#   C. SMETANA into its own venv, verified by SOLVING (step 38's second half)
#   D. the Phred-string filenames, listed before anything is moved
#
# ---------------------------------------------------------------------------
# A. ANI — THE LINK BETWEEN THE TWO CHAPTERS
#
#   Step 12 flagged "ANI: liquid MAGs vs culturomics isolates" as the join
#   between the metagenome and culturomics arms and it was never run. The
#   isolate genomes are in the culturomics manuscript directory:
#   Final_Genomes/{YGCA01, YGCA01_Clean, YGCA01_Contaminant, YGCA07, YGCA10,
#   YGCC01}.fasta plus four Clostridia references.
#
#   Those files are COPIED into $WS2 before use. The culturomics directory is
#   read and never written — the lesson from 38b, where using another project's
#   directory in place cost 2,435 files.
#
#   POSITIVE CONTROL, and it is a good one: this dataset contains a known
#   99.8% ANI pair — ETO_SED_MAG_002 (sediment) and ETO_LIQ_MAG_012 (culture),
#   the *Chakrabartyella* strain recovered from both. If the ANI run does not
#   return ~99.8% for that pair, the tool or the join is wrong and no other
#   number in section A may be believed. A negative control comes free: an
#   archaeal MAG against a bacterial isolate must produce no ANI value at all.
#
#   The audit's specific question: two isolates were cultured and never binned
#   (*Paraclostridium benzoelyticum*, *Terrisporobacter petrolearius*), and the
#   MAG set has Peptostreptococcales only in Anaerovoracaceae. If ANI confirms
#   no MAG matches those isolates, "cultivation recovered organisms the
#   metagenome missed" becomes a measured statement.
#
# ---------------------------------------------------------------------------
# B. THE 30 FREE VIRUSES INSIDE MAGs
#
#   39b's topology breakdown of the MAG-level geNomad run:
#       38  Provirus
#       29  No terminal repeats
#        1  DTR
#
#   The 38 proviruses were re-assigned and reported. The other 30 are calls of
#   FREE viruses on contigs that the binner placed INSIDE a MAG. Either
#   geNomad is over-calling (step 28 already caught it flagging chromosomal
#   regions in *Gaoshiqia*), or those genomes carry binned viral contigs that
#   inflate their size and gene counts.
#
#   The test is distributional, not anecdotal: concentrated in one or two
#   genomes means a binning problem in those genomes; spread thinly across many
#   means detector noise. The contig->MAG map from 39b is exact (2,549/2,549,
#   zero collisions), so this is a lookup, not an inference.
#
# ---------------------------------------------------------------------------
# C. SMETANA — INSTALLED, THEN MADE TO SOLVE
#
#   Into a NEW venv under $WS2/envs. Never into a shared env; never into the
#   culturomics pixi env. Then the same discipline as the pre-flight: an import
#   proves nothing, so it is asked to solve an LP and a MILP with known optima
#   and a wrong answer is reported WRONG rather than present.
#
# ---------------------------------------------------------------------------
# D. THE PHRED-STRING FILENAMES
#
#   Volume III: "stray files in workspace_2/ — Phred quality strings used as
#   filenames; quarantine script written, source untraced." A filename made of
#   Phred characters is a redirection that lost its target.
#
#   THIS MOVES NOTHING BY DEFAULT. It lists candidates with size and mtime so
#   the set can be inspected first. Set QUARANTINE=1 to move them into
#   $WS2/_quarantine/ — moved, never deleted, because a file whose origin is
#   unknown is a file whose importance is unknown.
#
# OUTPUT  40_integration/{INTEGRATION.txt, ani_matrix.tsv, viral_in_mags.tsv,
#                         phred_candidates.tsv}
# RUNTIME 5-20 min, dominated by the ANI all-vs-all
# USAGE   ws2 && bash scripts/40_integration.sh 2>&1 | tee logs/40.log
#         # QUARANTINE=1 to actually move the Phred files
#         # SKIP_SMETANA=1 if you would rather not pip install right now
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="40_integration"
F="$WS2/07_mags"; V="$WS2/18_viruses"
D="$WS2/40_integration"; mkdir -p "$D"
ISO="$D/isolates"; mkdir -p "$ISO"

echo "######## 40 — INTEGRATION JOBS ########"
echo "  started: $(date -Is)"
rm -f "$D/INTEGRATION.txt" "$D/ani_matrix.tsv" "$D/viral_in_mags.tsv" \
      "$D/phred_candidates.tsv" "$D/smetana_probe.txt"

# ---------------------------------------------------------------- locate things
CULT="$(find "$HOME/Documents" -maxdepth 4 -type d -name 'Final_Genomes' \
        2>/dev/null | sort | head -1)"
echo "  culturomics genomes : ${CULT:-<not found>}   [read-only]"

ANITOOL=""; ANIKIND=""
for t in skani fastANI; do
  p="$(command -v "$t" 2>/dev/null)"
  if [ -z "$p" ]; then
    p="$(find "$HOME" -maxdepth 8 -type f -name "$t" -perm -u+x 2>/dev/null \
         | sort | head -1)"
  fi
  if [ -n "$p" ]; then ANITOOL="$p"; ANIKIND="$t"; break; fi
done
echo "  ANI tool            : ${ANITOOL:-<none found>} ${ANIKIND:+($ANIKIND)}"
echo

# =============================================================== A. ANI
echo "######## A — ANI: MAGs vs CULTUROMICS ISOLATES ########"
NISO=0
if [ -n "$CULT" ]; then
  # copy, never use in place
  for f in "$CULT"/*.fasta "$CULT"/*.fa "$CULT"/*.fna; do
    [ -e "$f" ] || continue
    cp -n "$f" "$ISO/" 2>/dev/null
  done
  # the four Clostridia reference genomes, if present beside them
  REFD="$(dirname "$CULT")/Clostridia_References"
  if [ -d "$REFD" ]; then
    for f in "$REFD"/*.fna "$REFD"/*.fasta; do
      [ -e "$f" ] || continue
      cp -n "$f" "$ISO/" 2>/dev/null
    done
  fi
  NISO=$(find "$ISO" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fna' \
         -o -name '*.fasta' \) | wc -l)
fi
echo "  isolate/reference genomes copied in: $NISO"

BINDIR=""
for sub in bins final MAGs; do
  if compgen -G "$F/$sub/*.fa" >/dev/null; then BINDIR="$F/$sub"; break; fi
done
NMAG=$(ls -1 "$BINDIR"/*.fa 2>/dev/null | wc -l)
echo "  MAGs: $NMAG from $BINDIR"

if [ -n "$ANITOOL" ] && [ "$NISO" -gt 0 ] && [ "$NMAG" -gt 0 ]; then
  ls -1 "$BINDIR"/*.fa > "$D/mag_list.txt"
  find "$ISO" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fna' \
       -o -name '*.fasta' \) | sort > "$D/iso_list.txt"
  # both lists together: MAG-vs-MAG gives the 99.8% positive control for free
  cat "$D/mag_list.txt" "$D/iso_list.txt" > "$D/all_list.txt"
  echo "  running $ANIKIND on $(wc -l < "$D/all_list.txt") genomes ..."
  if [ "$ANIKIND" = "skani" ]; then
    # `triangle` emits a lower-triangular MATRIX (a count line, then one row per
    # genome). v1 asked for that and then parsed it as a pair table, so row
    # labels and numbers were read as genome names — "39 pairs" from 39 genomes
    # was the giveaway. `dist --ql --rl` emits the header'd pair table:
    #   Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query
    "$ANITOOL" dist --ql "$D/all_list.txt" --rl "$D/all_list.txt" \
      -t "${THREADS:-4}" -o "$D/skani_raw.tsv" > "$D/ani_tool.log" 2>&1
    RCA=$?
  else
    "$ANITOOL" --ql "$D/all_list.txt" --rl "$D/all_list.txt" \
      -t "${THREADS:-4}" -o "$D/fastani_raw.tsv" > "$D/ani_tool.log" 2>&1
    RCA=$?
  fi
  echo "  $ANIKIND exit $RCA"
else
  echo "  [SKIP] need an ANI tool, isolates and MAGs; one is missing."
fi
echo

# =============================================================== C. SMETANA
echo "######## C — SMETANA INTO ITS OWN VENV ########"
SMENV="$WS2/envs/smetana"
if [ "${SKIP_SMETANA:-0}" = "1" ]; then
  echo "  skipped (SKIP_SMETANA=1)"
elif [ -x "$SMENV/bin/python" ]; then
  echo "  venv already exists at $SMENV"
else
  python3 -m venv "$SMENV" > "$D/venv.log" 2>&1 \
    && echo "  venv created: $SMENV" || echo "  !! venv creation failed"
fi
if [ -x "$SMENV/bin/pip" ] && [ "${SKIP_SMETANA:-0}" != "1" ]; then
  "$SMENV/bin/pip" install --quiet --upgrade pip >> "$D/venv.log" 2>&1
  echo "  installing smetana + swiglpk (this is the only network step) ..."
  "$SMENV/bin/pip" install --quiet smetana swiglpk >> "$D/venv.log" 2>&1
  echo "  pip exit $?"
fi
# an import proves nothing; make it solve
if [ -x "$SMENV/bin/python" ]; then
  "$SMENV/bin/python" - > "$D/smetana_probe.txt" 2>&1 <<'PYX'
res = {}
for m in ("smetana", "reframed", "swiglpk", "optlang", "cobra"):
    try:
        __import__(m); res[m] = "present"
    except Exception as e:
        res[m] = f"MISSING ({type(e).__name__})"
for k, v in sorted(res.items()):
    print(f"MOD\t{k}\t{v}")
try:
    __import__("zzz_no_such_module_40"); print("MOD\tPROBE\tBROKEN")
except Exception:
    print("MOD\tPROBE\tok")
try:
    from optlang import glpk_interface as gi
    x = gi.Variable("x", lb=0); y = gi.Variable("y", lb=0)
    m1 = gi.Model(name="lp")
    m1.add([gi.Constraint(x + y, ub=4), gi.Constraint(x + 3 * y, ub=6)])
    m1.objective = gi.Objective(3 * x + 2 * y, direction="max")
    m1.optimize(); print(f"LP\t{m1.objective.value:.4f}")
    xi = gi.Variable("xi", lb=0, ub=2.5, type="integer")
    yi = gi.Variable("yi", lb=0, type="integer")
    m2 = gi.Model(name="milp")
    m2.add([gi.Constraint(xi + yi, ub=4), gi.Constraint(xi + 3 * yi, ub=6)])
    m2.objective = gi.Objective(3 * xi + 2 * yi, direction="max")
    m2.optimize(); print(f"MILP\t{m2.objective.value:.4f}")
except Exception as e:
    print(f"SOLVE\tFAILED\t{type(e).__name__}: {e}")
PYX
fi
echo

# =============================================================== the report
python3 - "$D" "$F" "$V" "$WS2" "$ANIKIND" "$BINDIR" "${QUARANTINE:-0}" <<'PY'
import sys, os, re, csv, glob, time
from collections import Counter, defaultdict

D, F, V, WS2, ANIKIND, BINDIR, QUAR = sys.argv[1:8]
W = 100
rep = []
STATUS = ['INCOMPLETE']


def say(s=''):
    print(s)
    rep.append(s)


def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t); say('=' * W)


def dump(st):
    STATUS[0] = st
    with open(os.path.join(D, 'INTEGRATION.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'INTEGRATION.txt')}")
    # v1 omitted this. A failed control printed "Stopping." and then kept going,
    # so the report said [FAIL] in the middle and STEP_STATUS: OK at the end.
    # A stop that does not stop is worse than no stop: it produces a document
    # that contradicts itself and a gate that passes.
    sys.exit(0)


def base(p):
    return re.sub(r'\.(fa|fna|fasta)$', '', os.path.basename(p))


# ---- crosswalk: working_name <-> mag_id, read from the file, not assumed
w2m, m2w, dom = {}, {}, {}
tp = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(tp):
    with open(tp, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            wn = (r.get('working_name') or '').strip()
            mi = (r.get('mag_id') or '').strip()
            dd = (r.get('domain') or '').strip()
            if wn:
                if mi:
                    w2m[wn] = mi; m2w[mi] = wn
                if dd:
                    dom[wn] = dd

# =====================================================================  A
rule('A. ANI — MAGs vs CULTUROMICS ISOLATES')
raw = None
for cand in ('skani_raw.tsv', 'fastani_raw.tsv'):
    p = os.path.join(D, cand)
    if os.path.exists(p) and os.path.getsize(p) > 0:
        raw = p; break
if not raw:
    say("  no ANI output — the tool did not run or produced nothing.")
    say("  Section A is NOT TESTED, which is different from 'no matches'.")
else:
    pairs = {}
    with open(raw, errors='replace') as fh:
        head = fh.readline()
        fh.seek(0)
        if ANIKIND == 'skani' and 'ANI' in head:
            for r in csv.DictReader(fh, delimiter='\t'):
                a = base(r.get('Ref_file') or r.get('Ref_name') or '')
                b = base(r.get('Query_file') or r.get('Query_name') or '')
                try:
                    v = float(r.get('ANI') or 0)
                except ValueError:
                    continue
                if a and b and a != b:
                    pairs[(a, b)] = v
        else:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 3:
                    continue
                a, b = base(f[0]), base(f[1])
                try:
                    v = float(f[2])
                except ValueError:
                    continue
                if a != b:
                    pairs[(a, b)] = v
    say(f"  {len(pairs):,} genome pairs with an ANI value ({ANIKIND})")
    say()

    # ---- POSITIVE CONTROL: the known 99.8% Chakrabartyella pair
    p1, p2 = m2w.get('ETO_SED_MAG_002'), m2w.get('ETO_LIQ_MAG_012')
    say("  CONTROL — the known 99.8% ANI pair must come back at ~99.8%")
    say(f"    ETO_SED_MAG_002 = {p1 or '?'}")
    say(f"    ETO_LIQ_MAG_012 = {p2 or '?'}")
    got = None
    if p1 and p2:
        got = pairs.get((p1, p2), pairs.get((p2, p1)))
    if got is None:
        say("  [FAIL] that pair is absent from the output. Either the join is")
        say("  wrong or the tool dropped it. No ANI number below is trustworthy.")
        dump('ANI_CONTROL_MISSING')
    say(f"    measured: {got:.2f}%")
    if abs(got - 99.8) > 0.5:
        say(f"  [FAIL] expected ~99.8%, got {got:.2f}%. Something is wrong with")
        say("  the tool, its settings, or the genome files. Stopping.")
        dump('ANI_CONTROL_FAILED')
    say("  [PASS] the control pair is recovered, so the run is measuring what")
    say("  it should.")
    say()

    mags = {base(p) for p in glob.glob(os.path.join(BINDIR, '*.fa'))}
    isos = {base(p) for p in glob.glob(os.path.join(D, 'isolates', '*'))}
    cross = sorted(((v, a, b) for (a, b), v in pairs.items()
                    if (a in mags and b in isos) or (a in isos and b in mags)),
                   key=lambda t: (-t[0], t[1], t[2]))
    say(f"  MAG x isolate pairs with any ANI value: {len(cross):,}")
    say()
    if cross:
        say(f"  {'ANI %':>8}  {'MAG':<28}{'ETO name':<20}{'isolate'}")
        say('-' * W)
        seen = set()
        for v, a, b in cross[:20]:
            mg, iso = (a, b) if a in mags else (b, a)
            if (mg, iso) in seen:
                continue
            seen.add((mg, iso))
            say(f"  {v:>8.2f}  {mg:<28}{w2m.get(mg, ''):<20}{iso}")
        say()
    say("  NOTE ON ABSENCE. skani only emits a pair when it clears its internal")
    say("  screening threshold (~80% ANI by default). So a MAG-isolate pair that")
    say("  does not appear was COMPUTED AND FELL BELOW THAT LINE — it is not")
    say("  'untested'. That distinction is the difference between a controlled")
    say("  negative and a gap, and this chapter has been bitten by it before.")
    say()
    hits95 = [t for t in cross if t[0] >= 95.0]
    say(f"  pairs at >=95% ANI (species level): {len(hits95)}")
    if not hits95:
        say()
        say("  >> NO MAG matches any isolate at species level.")
        say("     That is the audit's F-finding made quantitative: cultivation")
        say("     recovered organisms the metagenome did not. *Paraclostridium")
        say("     benzoelyticum* and *Terrisporobacter petrolearius* have no")
        say("     counterpart in the 29-genome catalogue, and the MAG set carries")
        say("     Peptostreptococcales only in Anaerovoracaceae.")
        say("     LIMIT: ANI below 95% means not the same species; it does not")
        say("     mean the organism is absent from the lagoon.")
    else:
        say()
        say("  >> At least one isolate has a species-level match in the MAG set.")
        say("     Those pairs are the direct bridge between the two chapters and")
        say("     should be reported with the ANI value and both genome names.")

# =====================================================================  B
rule('B. FREE VIRUSES ON CONTIGS THAT WERE BINNED INTO MAGs')
say("  39b split the MAG-level geNomad calls: 38 Provirus, 29 'No terminal")
say("  repeats', 1 DTR. The 38 were re-assigned. The other 30 are FREE virus")
say("  calls on contigs the binner put inside a MAG — either geNomad")
say("  over-calling (step 28 caught it flagging chromosomal regions in")
say("  Gaoshiqia) or genuine viral contigs inflating those genomes.")
say()
cmap = {}
cp = os.path.join(WS2, '39_loose_ends', 'provirus', 'contig_map.tsv')
if os.path.exists(cp):
    with open(cp, errors='replace') as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2:
                cmap[f[0]] = f[1]
say(f"  contig -> MAG map: {len(cmap):,} contigs (from 39b, 100% round trip)")

vs = sorted(glob.glob(os.path.join(V, '**', 'all_mags_contigs*virus_summary.tsv'),
                      recursive=True))
vs = [x for x in vs if 'provirus_genes' not in x]
if not cmap or not vs:
    say("  [SKIP] need 39b's contig_map.tsv and the geNomad summary. Run 39b first.")
else:
    rows = []
    with open(vs[0], errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            rows.append(r)
    cols = {k.lower(): k for k in (rows[0] if rows else {}) if k}
    tk = cols.get('topology'); sk = cols.get('seq_name') or (list(rows[0])[0] if rows else '')
    free = [r for r in rows if tk and 'provirus' not in str(r.get(tk, '')).lower()]
    say(f"  free-virus rows (topology != Provirus): {len(free)}")
    per = Counter(); unres = []
    with open(os.path.join(D, 'viral_in_mags.tsv'), 'w') as out:
        out.write('contig\tbin_id\tmag_id\ttopology\n')
        for r in free:
            sid = str(r.get(sk, '')).split('|')[0]
            b = sid.rpartition('~')[2] if '~' in sid else sid
            binid = cmap.get(b)
            if binid:
                per[binid] += 1
                out.write(f"{b}\t{binid}\t{w2m.get(binid,'')}\t{r.get(tk,'')}\n")
            else:
                unres.append(sid)
    say(f"  resolved to a MAG: {sum(per.values())}   unresolved: {len(unres)}")
    say()
    if per:
        say(f"  {'bin_id':<30}{'ETO name':<20}{'free virus contigs':>20}")
        say('-' * W)
        for b, n in sorted(per.items(), key=lambda kv: (-kv[1], kv[0])):
            say(f"  {b:<30}{w2m.get(b, ''):<20}{n:>20}")
        say()
        tot = sum(per.values())
        top = max(per.values())
        say(f"  genomes affected: {len(per)}   most in one genome: {top}"
            f"   ({100.0 * top / tot:.0f}% of the total)")
        say()
        # NO VERDICT FROM AN INVENTED THRESHOLD. 23c reported 46.6% against a
        # 40% line I had made up, and the lesson (Rule 8) was: report the
        # number, not the side of the line it landed on. So this prints the
        # distribution and names the reading each shape would support.
        nmag = len(glob.glob(os.path.join(BINDIR, '*.fa'))) or len(per)
        say(f"  affected genomes    : {len(per)} of {nmag}"
            f"  ({100.0 * len(per) / max(1, nmag):.0f}% of the catalogue)")
        say(f"  contigs             : {tot}")
        say(f"  most in one genome  : {top}  ({100.0 * top / tot:.0f}% of them)")
        counts = sorted(per.values(), reverse=True)
        say(f"  per-genome counts   : {', '.join(str(c) for c in counts[:12])}"
            + (' ...' if len(counts) > 12 else ''))
        say()
        say("  HOW TO READ IT — the threshold is not mine to set, so both")
        say("  readings are stated and the numbers above decide between them:")
        say("    * a few genomes holding most of the contigs => a BINNING problem")
        say("      in those genomes. Their size, gene count and CAZyme totals may")
        say("      be inflated by viral sequence and should be re-checked before")
        say("      being quoted.")
        say("    * a thin spread across many genomes => geNomad OVER-CALLING,")
        say("      consistent with step 28 where it flagged chromosomal regions")
        say("      in Gaoshiqia as MGEs. Then it is a caveat on the free-virus")
        say("      calls, not contamination of the catalogue.")
        say()
        say("  Either way these 30 are NOT proviruses and must not be added to")
        say("  the 38. `provirus_host_links.tsv` remains authoritative for")
        say("  n_proviruses; step 25 settled that and this does not reopen it.")
        say()

# =====================================================================  C
rule('C. SMETANA — INSTALLED AND MADE TO SOLVE')
sp = os.path.join(D, 'smetana_probe.txt')
if not os.path.exists(sp):
    say("  not probed (SKIP_SMETANA=1, or the venv could not be created).")
else:
    lp = milp = None
    broken = False
    for line in open(sp, errors='replace'):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'MOD' and len(f) > 2:
            if f[1] == 'PROBE':
                broken = (f[2] == 'BROKEN')
            else:
                say(f"    {f[1]:<12}{f[2]}")
        elif f[0] == 'LP' and len(f) > 1:
            try: lp = float(f[1])
            except ValueError: pass
        elif f[0] == 'MILP' and len(f) > 1:
            try: milp = float(f[1])
            except ValueError: pass
        elif f[0] == 'SOLVE':
            say(f"    solve failed: {' '.join(f[1:])[:80]}")
    say()
    if broken:
        say("  [FAIL] the negative control reported a nonexistent module as")
        say("  present. Nothing above may be believed.")
    else:
        say(f"  LP   = {lp}   (expected 12.0)")
        say(f"  MILP = {milp}   (expected 9.0)")
        if lp is not None and abs(lp - 12.0) < 1e-3 and \
           milp is not None and abs(milp - 9.0) < 1e-3:
            say("  [PASS] SMETANA's environment solves a MILP correctly. Step 38's")
            say("  simulation half is no longer blocked on tooling.")
        else:
            say("  [FAIL] the solver is absent or returns the wrong optimum.")
            say("  smetana will import and then die at the first optimisation.")

# =====================================================================  D
rule('D. THE PHRED-STRING FILENAMES')
say("  A filename made of Phred quality characters is a redirection that lost")
say("  its target. Volume III has carried this as 'source untraced'.")
say()
PH = re.compile(r'^[!-J]{20,}$')
cands = []
for dirpath, dirnames, filenames in os.walk(WS2):
    dirnames[:] = [d for d in dirnames
                   if d not in ('_quarantine', 'envs', '.git', 'tmp')]
    for fn in filenames:
        if PH.match(fn):
            p = os.path.join(dirpath, fn)
            try:
                st = os.stat(p)
                cands.append((os.path.relpath(p, WS2), st.st_size, st.st_mtime))
            except OSError:
                pass
cands.sort()
with open(os.path.join(D, 'phred_candidates.tsv'), 'w') as fh:
    fh.write('path\tbytes\tmtime\n')
    for p, s, m in cands:
        fh.write(f"{p}\t{s}\t{time.strftime('%Y-%m-%d %H:%M', time.localtime(m))}\n")
say(f"  candidates: {len(cands)}")
for p, s, m in cands[:15]:
    say(f"    {s:>12,}  {time.strftime('%Y-%m-%d %H:%M', time.localtime(m))}  {p[:60]}")
if len(cands) > 15:
    say(f"    ... and {len(cands) - 15} more (full list in phred_candidates.tsv)")
say()
if cands:
    dirs = Counter(os.path.dirname(p) or '.' for p, _, _ in cands)
    say("  where they are:")
    for dd, n in sorted(dirs.items(), key=lambda kv: (-kv[1], kv[0]))[:8]:
        say(f"    {n:>4}  {dd}")
    say()
    ts = sorted(m for _, _, m in cands)
    say(f"  created between {time.strftime('%Y-%m-%d %H:%M', time.localtime(ts[0]))}"
        f" and {time.strftime('%Y-%m-%d %H:%M', time.localtime(ts[-1]))}")
    say("  That window is the lead on the source: whichever script ran then")
    say("  wrote a redirection whose target expanded to a quality string.")
    say()
if QUAR == '1':
    q = os.path.join(WS2, '_quarantine')
    os.makedirs(q, exist_ok=True)
    moved = 0
    for p, _, _ in cands:
        src = os.path.join(WS2, p)
        dst = os.path.join(q, p.replace(os.sep, '__'))
        try:
            os.replace(src, dst); moved += 1
        except OSError as e:
            say(f"    could not move {p}: {e}")
    say(f"  MOVED {moved} file(s) to _quarantine/ — moved, never deleted.")
else:
    say("  NOTHING WAS MOVED. Re-run with QUARANTINE=1 once the list looks right.")
    say("  A file whose origin is unknown is a file whose importance is unknown.")

dump('OK')
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "report written" "$([ -s "$D/INTEGRATION.txt" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$D/INTEGRATION.txt" 2>/dev/null | cut -d' ' -f2-)
echo "  status: ${ST:-<none>}"
[ "${ST:-MISSING}" = "OK" ] || { echo "  [STOP] ${ST:-no status line}"; OK=1; }
[ "$OK" -eq 0 ] && runlog "40_integration complete" || runlog "**40_integration FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
