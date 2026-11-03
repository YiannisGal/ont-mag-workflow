#!/usr/bin/env bash
# =============================================================================
# 38_preflight_gems.sh  --  CAN step 38 RUN ON THIS MACHINE AT ALL?
#
# WHY A PRE-FLIGHT AND NOT THE STEP
#   Step 38 is genome-scale metabolic models (CarveMe or gapseq) plus SMETANA
#   cross-feeding simulation across 29 MAGs. It is the heaviest remaining item
#   and the only one whose answer is a FLUX FEASIBILITY rather than a marker
#   presence — which is exactly why it is worth doing and exactly why it is
#   expensive. Every hour spent on it is wasted if this machine has no MILP
#   solver, because both CarveMe's gap-filling and SMETANA's SMETANA score are
#   mixed-integer problems. cobrapy will import happily with no solver
#   installed and fail at the first `optimize()`, hours in.
#
#   So: answer the go/no-go question in five minutes, with a solve, not an
#   import. Nothing here writes a model, downloads a database, or touches a
#   shared conda environment.
#
# WHAT IT ESTABLISHES, IN ORDER
#   1. Which python interpreters exist (system + every conda env prefix), and
#      which of them can import a solver binding at all.
#   2. Whether any of them can SOLVE. An LP and then a MILP with a known
#      optimum are built and solved; a solver that returns the wrong answer is
#      reported as WRONG, not as present. This is the control — the whole
#      pre-flight is worthless if "solver found" is inferred from an import.
#   3. Whether the model-building tools exist (carveme, gapseq) and their own
#      hard dependencies (diamond for CarveMe; blast/hmmer/exonerate/R for
#      gapseq), and whether smetana exists.
#   4. Whether a NEW environment could be created if something is missing —
#      without modifying any shared env. This tests the writable location and
#      PyPI reachability, because "just pip install it" is an assumption on a
#      machine where ${CONDA_ENVS_SHARED} is not writable.
#   5. Whether the INPUTS are ready: 29 proteomes, their sizes, and the
#      community-size problem SMETANA faces (its global score is combinatorial;
#      29 members is not a run you start without knowing that).
#
# WHAT IT DOES NOT DO
#   It does not build a model, does not gap-fill, does not simulate, and does
#   not install anything. It ends in one of three verdicts — GO, GO AFTER
#   INSTALL (with the exact missing pieces), or NO-GO (with the reason) — so
#   that the decision to spend the hours is made on evidence.
#
# OUTPUT  38_preflight/PREFLIGHT.txt
# RUNTIME 1-5 min
# USAGE   ws2 && bash scripts/38_preflight_gems.sh 2>&1 | tee logs/38pre.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="38_preflight"
A="$WS2/09_annotation"; F="$WS2/07_mags"
D="$WS2/38_preflight"; mkdir -p "$D"

echo "######## 38 PRE-FLIGHT — GEMs + SMETANA ########"
echo "  started: $(date -Is)"
rm -f "$D/PREFLIGHT.txt" "$D/interpreters.tsv" "$D/tools.tsv"

# ---------------------------------------------------------------- interpreters
# Every python we could legitimately use, WITHOUT modifying any of them.
: "${CONDA_ROOTS:=${CONDA_ENVS_SHARED} /opt/conda/envs $HOME/miniconda3/envs}"
{
  command -v python3 2>/dev/null
  for root in $CONDA_ROOTS; do
    for p in "$root"/*/bin/python3; do [ -x "$p" ] && echo "$p"; done
  done
} 2>/dev/null | sort -u > "$D/interpreters.tsv"
NPY=$(wc -l < "$D/interpreters.tsv")
echo "  python interpreters found: $NPY"

# ---------------------------------------------------------------- CLI tools
for t in carveme carve gapseq smetana diamond blastp hmmsearch exonerate \
         barrnap bedtools Rscript git curl; do
  loc=""
  if command -v "$t" >/dev/null 2>&1; then
    loc="$(command -v "$t")"
  else
    for root in $CONDA_ROOTS; do
      for p in "$root"/*/bin/"$t"; do [ -x "$p" ] && { loc="$p"; break 2; }; done
    done
  fi
  printf '%s\t%s\n' "$t" "${loc:-MISSING}"
done | sort > "$D/tools.tsv"

# ---------------------------------------------------------------- R solvers
# gapseq is an R program and its fill step is a MILP solved through R, not
# through Python. Probing only python interpreters answers a narrower question
# than the one asked -- the 19b rpoB/tuf error: a control applied outside its
# domain of applicability. Ask R the same two problems.
: > "$D/r_solver.txt"
RS="$(command -v Rscript 2>/dev/null || true)"
if [ -n "$RS" ]; then
  "$RS" --vanilla -e '
    libs <- c("Rglpk","glpkAPI","cplexAPI","sybil","cobrar","Rsymphony")
    for (l in libs) cat(sprintf("MOD\t%s\t%s\n", l,
        ifelse(requireNamespace(l, quietly=TRUE), "present", "-")))
    if (requireNamespace("Rglpk", quietly=TRUE)) {
      obj <- c(3,2); mat <- matrix(c(1,1,1,3), nrow=2, byrow=TRUE)
      lp <- Rglpk::Rglpk_solve_LP(obj, mat, c("<=","<="), c(4,6), max=TRUE)
      cat(sprintf("LP\t%.4f\n", lp$optimum))
      mi <- Rglpk::Rglpk_solve_LP(obj, rbind(mat, c(1,0)), c("<=","<=","<="),
                                  c(4,6,2), types=c("I","I"), max=TRUE)
      cat(sprintf("MILP\t%.4f\n", mi$optimum))
    }' 2>/dev/null | grep -E '^(MOD|LP|MILP)' > "$D/r_solver.txt" || true
fi

# ---------------------------------------------------------------- gapseq evidence
# A tool that is not on $PATH may still have RUN here. gapseq's fill step turns
# <name>-draft.xml into <name>.xml; that pair is proof a MILP was solved on this
# machine, whatever the probes above say. Look for the evidence, not the binary.
find "$HOME" -maxdepth 6 -name '*-draft.xml' -not -path '*/.*' 2>/dev/null \
  | sort | head -40 > "$D/gapseq_evidence.txt" || true

echo
python3 - "$D" "$A" "$F" "$WS2" <<'PY'
import sys, os, re, glob, subprocess, json, shutil, tempfile

D, A, F, WS2 = sys.argv[1:5]
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


def dump(verdict):
    STATUS[0] = verdict
    with open(os.path.join(D, 'PREFLIGHT.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'PREFLIGHT.txt')}")


# =====================================================================  1
rule('1. SOLVER — probed by SOLVING, not by importing')
say("  Both CarveMe's gap-filling and SMETANA's score are MIXED-INTEGER")
say("  problems. cobrapy imports fine with no solver and dies at the first")
say("  optimize(). So each interpreter is asked to solve two problems with")
say("  known answers:")
say("      LP   : max 3x + 2y  s.t. x+y<=4, x+3y<=6, x,y>=0      -> 12.0 at (4,0)")
say("      MILP : same, x,y INTEGER, plus x<=2.5                  ->  9.0 at (2,1)")
say("  A solver that answers anything else is reported WRONG, not present.")
say()

PROBE = r'''
import json, sys
res = {"mods": [], "lp": None, "milp": None, "engine": None, "err": ""}
for m in ("cplex", "gurobipy", "swiglpk", "highspy", "pyscipopt", "mip",
          "optlang", "cobra", "reframed", "carveme", "smetana"):
    try:
        __import__(m)
        res["mods"].append(m)
    except Exception:
        pass
# negative control: a module that cannot exist must NOT be reported present
try:
    __import__("zzz_no_such_module_38")
    res["mods"].append("PROBE_BROKEN")
except Exception:
    pass
try:
    import optlang
    for name in ("cplex_interface", "gurobi_interface", "scipy_interface",
                 "glpk_interface", "coinor_cbc_interface"):
        try:
            iface = getattr(__import__("optlang." + name, fromlist=["x"]), "Model", None)
        except Exception:
            continue
        if iface is None:
            continue
        try:
            mod = __import__("optlang." + name, fromlist=["x"])
            Model, Var, Con, Obj = mod.Model, mod.Variable, mod.Constraint, mod.Objective
            # ---- LP
            x = Var("x", lb=0); y = Var("y", lb=0)
            mm = Model(name="lp")
            mm.add([Con(x + y, ub=4), Con(x + 3 * y, ub=6)])
            mm.objective = Obj(3 * x + 2 * y, direction="max")
            mm.optimize()
            res["lp"] = round(float(mm.objective.value), 4)
            # ---- MILP
            xi = Var("xi", lb=0, ub=2.5, type="integer")
            yi = Var("yi", lb=0, type="integer")
            m2 = Model(name="milp")
            m2.add([Con(xi + yi, ub=4), Con(xi + 3 * yi, ub=6)])
            m2.objective = Obj(3 * xi + 2 * yi, direction="max")
            m2.optimize()
            res["milp"] = round(float(m2.objective.value), 4)
            res["engine"] = name
            break
        except Exception as e:
            res["err"] = f"{name}: {type(e).__name__}: {e}"[:160]
except Exception as e:
    res["err"] = (res["err"] + " | optlang: " + str(e))[:200]
print("JSON38" + json.dumps(res))
'''

interps = [l.strip() for l in open(os.path.join(D, 'interpreters.tsv'))
           if l.strip()]
rows = []
for py in interps:
    env = os.path.basename(os.path.dirname(os.path.dirname(py)))
    try:
        p = subprocess.run([py, '-c', PROBE], capture_output=True, text=True,
                           timeout=180)
        line = next((l for l in p.stdout.split('\n') if l.startswith('JSON38')), '')
        r = json.loads(line[6:]) if line else {}
    except Exception as e:
        r = {'err': f'{type(e).__name__}: {e}'[:120]}
    rows.append((env, py, r))

say(f"  {'environment':<26}{'solver modules present':<40}{'LP':>7}{'MILP':>7}{'engine':>16}")
say('-' * W)
GO = []
probe_broken = False
for env, py, r in sorted(rows, key=lambda t: (t[0], t[1])):
    mods = [m for m in r.get('mods', [])
            if m in ('cplex', 'gurobipy', 'swiglpk', 'highspy', 'pyscipopt',
                     'mip', 'optlang')]
    if 'PROBE_BROKEN' in r.get('mods', []):
        probe_broken = True
    lp, milp = r.get('lp'), r.get('milp')
    say(f"  {env[:25]:<26}{(', '.join(mods) or '-')[:39]:<40}"
        f"{(f'{lp:.1f}' if lp is not None else '-'):>7}"
        f"{(f'{milp:.1f}' if milp is not None else '-'):>7}"
        f"{(r.get('engine') or '-')[:15]:>16}")
    if lp is not None and abs(lp - 12.0) < 1e-3 and \
       milp is not None and abs(milp - 9.0) < 1e-3:
        GO.append((env, py, r))
say()

# ---- controls on the probe itself
say("  CONTROL — negative: a deliberately nonexistent module was probed in")
say("  every interpreter. " + ("!! IT WAS REPORTED PRESENT — the probe is broken"
                               " and nothing above may be believed."
                               if probe_broken else
                               "It was reported absent everywhere, so 'missing'"))
if not probe_broken:
    say("  above means missing and not merely unprobed.")
say()
wrong = [(e, r) for e, p, r in rows
         if (r.get('lp') is not None and abs(r['lp'] - 12.0) >= 1e-3)
         or (r.get('milp') is not None and abs(r['milp'] - 9.0) >= 1e-3)]
if wrong:
    say("  !! WRONG ANSWERS — these environments have a solver that does not")
    say("     return the known optimum. Do not use them:")
    for e, r in sorted(wrong):
        say(f"       {e:<26}LP={r.get('lp')}  MILP={r.get('milp')}  "
            f"(expected 12.0 / 9.0)")
    say()

if probe_broken:
    GO = []
    say("  >> SECTION 1 IS VOID. The probe reported a module that cannot exist,")
    say("     so every 'present' above is unreliable and no solver claim may be")
    say("     made from this run. Fix the probe before believing any of it.")
elif GO:
    say(f"  >> {len(GO)} environment(s) SOLVE both problems correctly:")
    for e, p, r in GO:
        say(f"       {e:<26}{r.get('engine')}   {p}")
    say("     A MILP solver exists on this machine. Step 38 is not blocked on")
    say("     the solver question.")
else:
    say("  >> NO environment on this machine solves a MILP. Neither CarveMe's")
    say("     gap-filling nor SMETANA can run as things stand. This is THE")
    say("     blocking fact and it is now established in five minutes rather")
    say("     than after a day of model building.")
say()

# =====================================================================  1b
rule('1b. R — GAPSEQ IS AN R PROGRAM AND PYTHON IS NOT THE ONLY PLACE A '
     'SOLVER LIVES')
say("  Section 1 asked python. gapseq's `fill` step is a MILP solved through R")
say("  (glpkAPI / cplexAPI via sybil or cobrar), so a python-only probe answers")
say("  a narrower question than the one asked — and reporting its answer as")
say("  'no solver on this machine' would be a control applied outside its")
say("  domain of applicability, which is the 19b rpoB/tuf error.")
say()
rmods, rlp, rmilp = [], None, None
rp = os.path.join(D, 'r_solver.txt')
if os.path.exists(rp):
    for line in open(rp, errors='replace'):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'MOD' and len(f) > 2 and f[2] == 'present':
            rmods.append(f[1])
        elif f[0] == 'LP' and len(f) > 1:
            try: rlp = float(f[1])
            except ValueError: pass
        elif f[0] == 'MILP' and len(f) > 1:
            try: rmilp = float(f[1])
            except ValueError: pass
if not os.path.exists(rp) or os.path.getsize(rp) == 0:
    say("  Rscript not found, or R answered nothing. R is NOT ruled out — it was")
    say("  not successfully probed, which is a different statement.")
else:
    say(f"  R solver/model libraries present : {', '.join(rmods) if rmods else 'none of those probed'}")
    say(f"  R LP   : {rlp if rlp is not None else 'not exercised'}   (expected 12.0)")
    say(f"  R MILP : {rmilp if rmilp is not None else 'not exercised'}   (expected 9.0)")
R_SOLVES = (rlp is not None and abs(rlp - 12.0) < 1e-3
            and rmilp is not None and abs(rmilp - 9.0) < 1e-3)
if R_SOLVES:
    say("  >> R SOLVES A MILP correctly. gapseq's gap-filling is not blocked.")
elif rmods:
    say("  >> R carries solver libraries but the MILP was not exercised (Rglpk")
    say("     absent). Presence is not a solve — treat this as UNCONFIRMED.")
say()

# ---- evidence beats a PATH lookup
ev = []
evp = os.path.join(D, 'gapseq_evidence.txt')
if os.path.exists(evp):
    ev = [l.strip() for l in open(evp, errors='replace') if l.strip()]
say("  EVIDENCE CHECK — has a gap-fill ever completed on this machine?")
say("  gapseq turns <name>-draft.xml into <name>.xml, and that transition IS the")
say("  MILP. A tool absent from $PATH may still have run.")
filled = []
for d0 in ev:
    final = d0.replace('-draft.xml', '.xml')
    if os.path.exists(final):
        filled.append(final)
if filled:
    say(f"  >> {len(filled)} gap-FILLED model(s) found on disk:")
    for f0 in filled[:6]:
        say(f"       {f0}")
    say("     A MILP has demonstrably been solved here. If sections 1 and 1b")
    say("     both say 'no solver', they are looking in the wrong place — find")
    say("     the R installation and the gapseq clone that produced these before")
    say("     installing anything.")
elif ev:
    say(f"  >> {len(ev)} draft model(s) found but no gap-filled counterpart —")
    say("     gapseq's find/draft stages ran and `fill` did not. That is the")
    say("     stage that needs the solver.")
    for f0 in ev[:6]:
        say(f"       {f0}")
else:
    say("  >> no gapseq draft models anywhere under $HOME. No evidence either way.")
say()

# =====================================================================  2
rule('2. MODEL-BUILDING AND SIMULATION TOOLS')
tools = {}
with open(os.path.join(D, 'tools.tsv')) as fh:
    for line in fh:
        k, _, v = line.rstrip('\n').partition('\t')
        tools[k] = v
say(f"  {'tool':<12}{'needed for':<50}{'location'}")
say('-' * W)
NEED = [
    ('carveme',   'CarveMe route: build the GEMs'),
    ('carve',     'CarveMe route: the actual CLI name'),
    ('diamond',   'CarveMe: protein->reaction mapping (HARD dep)'),
    ('gapseq',    'gapseq route: build the GEMs'),
    ('blastp',    'gapseq: pathway prediction (HARD dep)'),
    ('exonerate', 'gapseq: sequence handling (HARD dep)'),
    ('bedtools',  'gapseq: interval handling (HARD dep)'),
    ('barrnap',   'gapseq: rRNA (HARD dep)'),
    ('Rscript',   'gapseq: the prediction scripts are R (HARD dep)'),
    ('smetana',   'the actual cross-feeding simulation'),
    ('hmmsearch', 'already in use here; listed for completeness'),
]
for t, why in NEED:
    loc = tools.get(t, 'MISSING')
    say(f"  {t:<12}{why[:49]:<50}{'MISSING' if loc == 'MISSING' else loc[:38]}")
say()
have_carve = tools.get('carve', 'MISSING') != 'MISSING' or \
             tools.get('carveme', 'MISSING') != 'MISSING'
have_gapseq = tools.get('gapseq', 'MISSING') != 'MISSING'
have_smet = tools.get('smetana', 'MISSING') != 'MISSING'
gapseq_deps = [t for t in ('blastp', 'exonerate', 'bedtools', 'barrnap', 'Rscript')
               if tools.get(t, 'MISSING') == 'MISSING']

if have_gapseq and not gapseq_deps:
    say("  >> gapseq is present WITH all its hard dependencies. The audit noted")
    say("     it was already used for isolate auxotrophy/complementarity, so the")
    say("     route with the least new tooling is gapseq for the models.")
elif have_gapseq:
    say(f"  >> gapseq is present but is MISSING: {', '.join(gapseq_deps)}.")
    say("     gapseq fails partway through without these; it does not warn first.")
if have_carve:
    say("  >> CarveMe is present." + ("" if tools.get('diamond', 'MISSING') != 'MISSING'
                                      else " But diamond is MISSING, and CarveMe"
                                           " cannot map proteins without it."))
if not have_carve and not have_gapseq:
    say("  >> NEITHER model builder is installed. Nothing can produce a GEM here")
    say("     without an install.")
if not have_smet:
    say("  >> smetana is MISSING. Even with models in hand there is nothing to")
    say("     run the cross-feeding simulation with.")
say()

# =====================================================================  3
rule('3. COULD THE MISSING PIECES BE INSTALLED — WITHOUT TOUCHING A SHARED ENV?')
say("  The shared conda envs are NOT writable and one has already had its")
say("  metadata broken by a write attempt. So the only acceptable install is a")
say("  NEW environment in a writable location. That is tested, not assumed.")
say()
cand = [os.path.join(WS2, 'envs'), os.path.join(WS2, 'tmp'), '${SCRATCH_ROOT}']
writable = []
for c in cand:
    try:
        os.makedirs(c, exist_ok=True)
        t = tempfile.NamedTemporaryFile(dir=c, delete=True)
        t.close()
        writable.append(c)
    except Exception:
        pass
say(f"  writable locations for a new env : "
    f"{', '.join(writable) if writable else 'NONE FOUND'}")
for c in writable:
    try:
        st = os.statvfs(c)
        say(f"    {c:<40}{st.f_bavail * st.f_frsize / 1e9:>8.1f} GB free")
    except Exception:
        pass
say()
say(f"  python -m venv available          : "
    f"{'yes' if subprocess.run([sys.executable, '-m', 'venv', '--help'], capture_output=True).returncode == 0 else 'NO'}")
net = 'not tested'
try:
    import urllib.request
    urllib.request.urlopen('https://pypi.org/simple/', timeout=15).read(64)
    net = 'reachable'
except Exception as e:
    code = getattr(e, 'code', None)
    net = (f'HTTP {code} — host reachable, index refused the request'
           if code else f'NOT reachable ({type(e).__name__})')
say(f"  PyPI                              : {net}")
mamba = shutil.which('mamba') or shutil.which('conda')
say(f"  conda/mamba on PATH               : {mamba or 'MISSING'}")
say()
say("  NOTE. smetana and carveme are pip-installable; gapseq is not (it is a")
say("  git clone plus R packages). A solver is the harder half: GLPK arrives")
say("  with `pip install swiglpk`, HiGHS with `pip install highspy`, and both")
say("  do MILP. CPLEX and Gurobi need licences. If section 1 found no solver,")
say("  swiglpk or highspy into a NEW venv is the cheapest way to one — and it")
say("  touches nothing shared.")
say()

# =====================================================================  4
rule('4. ARE THE INPUTS READY, AND HOW BIG IS THE JOB?')
# A GLOB IS A CLAIM ABOUT WHICH FILES EXIST. `09_annotation/*/*.faa` matched 58
# files for 29 MAGs on the real workspace, because Bakta writes BOTH
# <name>.faa and <name>.hypotheticals.faa into each directory — so every genome
# was counted twice and its hypotheticals file, being small, was then reported
# as a "fragmentary genome". Take the primary proteome by exact name only.
faas = []
for d in sorted(os.listdir(A)) if os.path.isdir(A) else []:
    f = os.path.join(A, d, d + '.faa')
    if os.path.exists(f):
        faas.append(f)
extra = sorted(set(glob.glob(os.path.join(A, '*', '*.faa'))) - set(faas))
sizes = []
for f in faas:
    n = 0
    try:
        with open(f, errors='replace') as fh:
            for line in fh:
                if line.startswith('>'):
                    n += 1
    except Exception:
        pass
    sizes.append((os.path.basename(os.path.dirname(f)), n))
say(f"  proteomes available for model building : {len(sizes)}"
    f"   (<dir>/<dir>.faa only)")
if extra:
    say(f"  {len(extra)} other .faa file(s) in the same directories were NOT counted")
    say(f"  as genomes — e.g. {os.path.basename(extra[0])}. Bakta writes a")
    say("  .hypotheticals.faa beside each proteome; counting it doubles the")
    say("  genome count and then reports it as a fragmentary genome.")
# cross-check against the frozen catalogue rather than trusting the directory
xw = os.path.join(F, 'MAG_CROSSWALK.tsv')
ncat = 0
if os.path.exists(xw):
    with open(xw, errors='replace') as fh:
        ncat = max(0, sum(1 for _ in fh) - 1)
if ncat:
    say(f"  MAG_CROSSWALK.tsv rows (the frozen catalogue)  : {ncat}")
    if len(sizes) != ncat:
        say(f"  !! MISMATCH: {len(sizes)} proteome directories against {ncat}")
        say("     catalogue rows. Resolve this BEFORE building models — a GEM set")
        say("     that does not correspond to the catalogue cannot be reported")
        say("     against it. Directories with no catalogue row are pre-DAS_Tool")
        say("     bins or an older binning pass.")
    else:
        say("  >> proteome directories and catalogue rows agree.")
if sizes:
    ns = sorted(s[1] for s in sizes)
    say(f"  proteins per genome: min {ns[0]:,}  median {ns[len(ns) // 2]:,}  "
        f"max {ns[-1]:,}  total {sum(ns):,}")
    thin = [g for g, n in sorted(sizes) if n < 800]
    if thin:
        say(f"  {len(thin)} genome(s) under 800 proteins — a GEM built from a")
        say("  fragmentary proteome gap-fills into a fictional organism:")
        say(f"    {', '.join(thin[:8])}"
            + (f"  ... and {len(thin) - 8} more" if len(thin) > 8 else ""))
        say("  Decide explicitly whether these enter the community model.")
say()
n = len(sizes) or 29
say("  SMETANA's scale problem, stated before the run and not after:")
say(f"    - `detailed` mode enumerates species pairs      : {n * (n - 1) // 2:,} pairs")
say(f"    - `global` mode samples subcommunities of the whole {n}-member set;")
say("      its cost grows with community size and it is routinely run on")
say(f"      communities of 3-10, not {n}. Expect to define subcommunities on a")
say("      biological basis (the sulfate-reducer / methanogen / fermenter")
say(f"      partitions this chapter already has) rather than running all {n}.")
say(f"    - model building alone: roughly {n} x 2-10 min = "
    f"{n * 2 / 60:.1f}-{n * 10 / 60:.1f} h for CarveMe, longer for gapseq.")
say()
say("  A MEDIA FILE IS ALSO REQUIRED and none has been specified. SMETANA's")
say("  answer depends on it entirely: a rich medium makes everything")
say("  independent and every cross-feeding score collapses. For this system the")
say("  medium must be anoxic, sulfate-replete and low in fixed carbon, or the")
say("  result will not describe Etoliko.")
say()

# =====================================================================
rule('VERDICT')
missing = []
SOLVER_OK = bool(GO) or R_SOLVES or bool(filled)
if not SOLVER_OK:
    missing.append('a working MILP solver')
elif not GO and (R_SOLVES or filled):
    say("  NOTE. No PYTHON solver, but a MILP is solvable here in R (or has")
    say("  demonstrably been solved). gapseq is therefore the route with the")
    say("  fewest new parts. SMETANA is python and still needs its own solver.")
    say()
if not (have_carve or have_gapseq):
    missing.append('a model builder (carveme or gapseq)')
elif have_gapseq and gapseq_deps and not have_carve:
    missing.append('gapseq dependencies: ' + ', '.join(gapseq_deps))
if not have_smet:
    missing.append('smetana')
if not faas:
    missing.append('proteomes (09_annotation/*/*.faa)')

if probe_broken:
    missing.insert(0, 'a trustworthy solver probe (the negative control failed)')

if not missing:
    verdict = 'GO'
    say("  GO. Solver, model builder and simulator are all present and the")
    say("  solver has been shown to return the correct MILP optimum. The")
    say("  remaining decisions are scientific, not technical: which genomes")
    say("  enter the community, which subcommunities SMETANA is run on, and")
    say("  what medium is used.")
elif not faas or (not SOLVER_OK and not writable):
    verdict = 'NO-GO'
    say("  NO-GO. " + '; '.join(missing) + '.')
    if not faas:
        say("  The blocking fact is the INPUT, not the tooling: there are no")
        say("  proteomes under 09_annotation/*/*.faa, and a GEM is built from a")
        say("  proteome. No install fixes this — step 09 has to have produced its")
        say("  output first. Check that $WS2 is set (an unset $WS2 makes every")
        say("  path start at /) before concluding the annotations are missing.")
    else:
        say("  With no writable location for a new environment and no solver,")
        say("  step 38 cannot be made to run here.")
    say("  Step 38 is legitimately optional; the chapter's cross-feeding argument")
    say("  currently rests on marker presence, and that should be stated as the")
    say("  limitation it is rather than papered over.")
else:
    verdict = 'GO-AFTER-INSTALL'
    say("  GO AFTER INSTALL. Missing: " + '; '.join(missing) + '.')
    say("  Everything missing can go into a NEW environment in a writable")
    say("  location — nothing shared is touched. Install, re-run this")
    say("  pre-flight, and only start step 38 when it says GO.")
say()
say("  This pre-flight builds no model and simulates nothing. Its only claim is")
say("  whether the work can start.")
dump(verdict)
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "pre-flight completed" "$RC" || OK=1
check "report written" "$([ -s "$D/PREFLIGHT.txt" ] && echo 0 || echo 1)" || OK=1
ST=$(grep -m1 '^STEP_STATUS:' "$D/PREFLIGHT.txt" 2>/dev/null | cut -d' ' -f2-)
echo "  verdict: ${ST:-<none>}"
case "${ST:-MISSING}" in
  GO|GO-AFTER-INSTALL|NO-GO) echo "  [PASS] a verdict was reached" ;;
  *) echo "  [FAIL] no verdict — the pre-flight itself did not complete"; OK=1 ;;
esac
[ "$OK" -eq 0 ] && runlog "38_preflight verdict: $ST" || runlog "**38_preflight FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
