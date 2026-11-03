#!/usr/bin/env bash
# 74_cys_control_panel.sh -- the control's failures ARE the finding.
#
#   bash scripts/run.sh scripts/74_cys_control_panel_XXXXXXX.sh
#
# 🔴 STEP 73 ABORTED ON AN ARBITRARY NUMBER, AND THE NUMBER IT ABORTED ON IS THE
# RESULT.
#
#     control K01883 ... 20 of 29 genome(s) pass.
#     !! the control still fails. Nothing below may be believed.
#
# My gate was `npass < len(best) * 0.7`. 0,7 x 29 = 20,3. Twenty passed. The run
# stopped one tenth of a genome short of a threshold I invented, and printed
# nothing about cysteine.
#
# Worse: 20/29 is not a broken parse. 1.509.145 lines parsed, 29 genome prefixes,
# rule-115 clean, and a housekeeping enzyme found in twenty genomes. A parse that
# were broken would find it in none -- which is exactly what step 72 showed, and
# why that abort was right and this one was not.
#
# AND K01883 IS THE WRONG SINGLE CONTROL FOR THIS PARTICULAR QUESTION. Canonical
# cysteinyl-tRNA synthetase is precisely the enzyme that organisms using the
# INDIRECT route can do without: they charge tRNA(Cys) with O-phosphoserine using
# SepRS and convert it in place with SepCysS. So "no K01883" is a candidate
# SIGNAL, not a candidate failure -- and I used it as the gate that suppressed the
# signal.
#
#   > Rule 116. A control is per item, not per run. "70% must pass" discards the
#   > items that did pass, and the items that failed may be the result.
#
#   > Rule 117. Never use a gene that is part of the hypothesis as the control for
#   > the search that tests the hypothesis.
#
# THIS STEP THEREFORE:
#   A. controls with a PANEL of four aminoacyl-tRNA synthetases that no route
#      around exists for, and reports the verdict PER GENOME. A genome that passes
#      three of four has a working parse whatever its neighbours do.
#   B. reports the cysteine panel for every genome that passed its own control,
#      including BOTH CysRS classes, and says plainly which genomes were excluded.
#   C. cross-tabulates the thing that matters: among archaea, does absence of
#      canonical CysRS coincide with presence of SepRS + SepCysS? That pairing is
#      the signature of the indirect route, and it is testable here.
#
# READ-ONLY except 74_cys/. One pass over all_hits.tbl (408 MB), a minute or two.
# Rule 101 honoured: all_hits.tbl only, never the chunks as well.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/74_cys"; mkdir -p "$OUT"
exec > >(tee "$OUT/74.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
ok(){ printf '  [OK]   %s\n' "$*"; }
no(){ printf '  [--]   %s\n' "$*"; }

echo "74_cys_control_panel -- per-genome control, then the cysteine routes"
echo "WS2 = $WS2"
echo "started $(date -Is)"

AH="$WS2/45_kofam_full/all_hits.tbl"
KOLIST="$WS2/36_gapclose/ko_list"
MASTER="$WS2/47_master/MASTER.tsv"
[ -s "$AH" ] || { no "missing $AH"; echo "STEP_STATUS: FAILED"; exit 1; }

python3 - "$AH" "$KOLIST" "$MASTER" "$OUT" <<'PYEOF'
import sys, os, re, csv, collections
ah, kolist, master, out = sys.argv[1:5]

# ---------------- what we look for, and why each one is here
CONTROL = [
    ("K01875", "seryl-tRNA synthetase"),        # universal, no bypass
    ("K01869", "leucyl-tRNA synthetase"),       # universal
    ("K01872", "alanyl-tRNA synthetase"),       # universal
    ("K01885", "glutamyl-tRNA synthetase"),     # universal
]
CYSRS = [("K01883", "cysteinyl-tRNA synthetase"),
         ("K01884", "cysteinyl-tRNA synthetase, unknown class")]
BACT = [("K01738", "cysK, cysteine synthase"),
        ("K00640", "cysE, serine O-acetyltransferase"),
        ("K10150", "cysM, cysteine synthase B")]
ARCH = [("K07587", "sepS, O-phosphoseryl-tRNA synthetase"),
        ("K06868", "sepCysS, Sep-tRNA:Cys-tRNA synthetase")]
ALL = CONTROL + CYSRS + BACT + ARCH
WANT = {k for k, _ in ALL}

thr, defn = {}, {}
with open(kolist) as fh:
    head = fh.readline().rstrip("\n").split("\t")
    it = head.index("threshold") if "threshold" in head else 1
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) <= it or p[0] not in WANT:
            continue
        try: thr[p[0]] = float(p[it])
        except ValueError: thr[p[0]] = None
        defn[p[0]] = p[-1]
missing = [k for k in WANT if k not in thr]
if missing:
    print("    !! not in ko_list, so untestable (rule 56): %s" % ", ".join(sorted(missing)))

M, dom = {}, {}
for r in csv.DictReader(open(master), delimiter="\t"):
    wn = r.get("working_name")
    if wn:
        M[wn] = (r.get("mag_id", "") or "").replace("ETO_", "")
        dom[wn] = r.get("domain", "")

KO_I, SC_I = 2, 5
best = collections.defaultdict(dict)
n = 0
with open(ah, errors="replace") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.split()
        if len(f) <= SC_I:
            continue
        ko = f[KO_I]
        if ko not in WANT:
            continue
        n += 1
        g = f[0]
        if "|" not in g:
            continue
        b = g.split("|", 1)[0]
        try: sc = float(f[SC_I])
        except ValueError: continue
        if sc > best[b].get(ko, -1.0):
            best[b][ko] = sc
print("    %d relevant hit line(s); %d genome(s) seen" % (n, len(best)))
for b in M:
    best.setdefault(b, {})


def passes(b, k):
    t = thr.get(k)
    return t is not None and best[b].get(k, -1.0) >= t


# ---------------------------------------------------------------- A
print()
print("================================================================")
print("A. THE CONTROL PANEL, PER GENOME")
print("================================================================")
for k, nm in CONTROL:
    print("      %-9s %-46s thr %s" % (k, nm[:46], thr.get(k)))
print()
hdr = "      %-20s %-9s" % ("genome", "domain")
for k, _ in CONTROL:
    hdr += "%9s" % k
hdr += "%8s  %s" % ("passed", "verdict")
print(hdr)
good, bad = [], []
order = sorted(best, key=lambda b: (dom.get(b, "") != "Archaea", M.get(b, b)))
for b in order:
    p = sum(1 for k, _ in CONTROL if passes(b, k))
    line = "      %-20s %-9s" % (M.get(b, b), dom.get(b, "?")[:9])
    for k, _ in CONTROL:
        sc = best[b].get(k, -1.0)
        line += "%9s" % (("%.0f" % sc) if passes(b, k) else
                         ("(%.0f)" % sc if sc > 0 else "-"))
    v = "usable" if p >= 3 else ("weak" if p == 2 else "EXCLUDED")
    line += "%8d  %s" % (p, v)
    print(line)
    (good if p >= 3 else bad).append(b)
print()
print("    %d genome(s) usable, %d excluded." % (len(good), len(bad)))
print("    This is a PER-GENOME verdict. Step 73 asked whether 70% of genomes")
print("    passed one control and stopped when 20 of 29 did. Twenty genomes with a")
print("    working parse are twenty answers, and the nine without are a question.")

# ---------------------------------------------------------------- B
print()
print("================================================================")
print("B. CYSTEINE: CHARGING, BACTERIAL SYNTHESIS, ARCHAEAL INDIRECT ROUTE")
print("================================================================")
hdr = "      %-20s %-9s" % ("genome", "domain")
for k, _ in CYSRS + BACT + ARCH:
    hdr += "%9s" % k
print(hdr)
rows = []
for b in order:
    if b not in good:
        continue
    line = "      %-20s %-9s" % (M.get(b, b), dom.get(b, "?")[:9])
    rec = {"mag": M.get(b, b), "domain": dom.get(b, "?"), "bin": b}
    for k, _ in CYSRS + BACT + ARCH:
        sc = best[b].get(k, -1.0)
        line += "%9s" % (("%.0f" % sc) if passes(b, k) else
                         ("(%.0f)" % sc if sc > 0 else "-"))
        rec[k] = "%.1f" % sc if sc > 0 else ""
        rec[k + "_pass"] = "yes" if passes(b, k) else "no"
    rows.append(rec)
    print(line)
print()
print("    bare = passed its ko_list threshold; (brackets) = below it; dash = no hit")
if bad:
    print("    excluded for failing their own control: %s"
          % ", ".join(M.get(b, b) for b in bad))

# ---------------------------------------------------------------- C
print()
print("================================================================")
print("C. THE PAIRING THAT DECIDES IT")
print("================================================================")
print("    An organism using the INDIRECT route charges tRNA(Cys) with")
print("    O-phosphoserine (SepRS, K07587) and converts it in place (SepCysS,")
print("    K06868). Such an organism can lack canonical CysRS entirely. So the")
print("    signature is: no K01883/K01884, but K07587 AND K06868 present.")
print()
print("      %-20s %-9s %-14s %-14s %s"
      % ("genome", "domain", "CysRS", "SepRS+SepCysS", "reading"))
for r in rows:
    crs = any(r[k + "_pass"] == "yes" for k, _ in CYSRS)
    sep = [k for k, _ in ARCH if r[k + "_pass"] == "yes"]
    syn = any(r[k + "_pass"] == "yes" for k, _ in BACT)
    if not crs and len(sep) == 2:
        read = "indirect route -- CysRS not needed"
    elif crs and syn:
        read = "direct: charges and synthesises"
    elif crs and not syn:
        read = "charges; synthesis route not detected"
    elif not crs and sep:
        read = "partial indirect (%d of 2)" % len(sep)
    else:
        read = "neither route detected"
    print("      %-20s %-9s %-14s %-14s %s"
          % (r["mag"], r["domain"][:9], "yes" if crs else "NO",
             "%d of 2" % len(sep), read))

arc = [r for r in rows if str(r["domain"]).startswith("Arch")]
bac = [r for r in rows if not str(r["domain"]).startswith("Arch")]
print()
for k, nm in CYSRS + BACT + ARCH:
    a = sum(1 for r in arc if r[k + "_pass"] == "yes")
    b_ = sum(1 for r in bac if r[k + "_pass"] == "yes")
    print("      %-9s %-40s archaea %2d/%-3d bacteria %2d/%-3d"
          % (k, nm[:40], a, len(arc), b_, len(bac)))

mtp = next((r for r in rows if r["mag"] == "LIQ_MAG_003"), None)
print()
if mtp is None:
    print("    MTP4 is not in the usable set -- see section A for why.")
else:
    crs = any(mtp[k + "_pass"] == "yes" for k, _ in CYSRS)
    sep = sum(1 for k, _ in ARCH if mtp[k + "_pass"] == "yes")
    syn = any(mtp[k + "_pass"] == "yes" for k, _ in BACT)
    print("    MTP4: CysRS %s | SepRS+SepCysS %d of 2 | bacterial synthesis %s"
          % ("present" if crs else "ABSENT", sep, "present" if syn else "absent"))
    if sep == 2 and not syn:
        print("    >>> the indirect archaeal route is complete and the bacterial one")
        print("        is absent. 19b's Cys = 0,00 is a MARKER-SET GAP: routes.tsv")
        print("        knows only cysK/cysE/cysM. Genome and model agree, and")
        print("        fig06's empty MTP4 column is a result.")
    elif syn:
        print("    >>> the bacterial route is present after all, so 19b's zero is a")
        print("        scoring question rather than a route question.")
    else:
        print("    >>> no complete route of either kind. Cysteine biosynthesis in")
        print("        MTP4 stays unresolved and fig06's empty column rests on the")
        print("        model alone. Say exactly that in the text.")

with open(os.path.join(out, "cys_panel.tsv"), "w", newline="") as fh:
    cols = (["mag", "domain", "bin"] + [k for k, _ in CYSRS + BACT + ARCH]
            + [k + "_pass" for k, _ in CYSRS + BACT + ARCH])
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow(r)
print()
print("    wrote %s" % os.path.join(out, "cys_panel.tsv"))
PYEOF

say "D. HOUSEKEEPING NOTICED WHILE READING THAT DIRECTORY"
cat <<'EOT'
  45_kofam_full holds all_proteins.faa at 432 MB and all_hits.tbl at 408 MB, plus
  a 24-file chunks/ directory that duplicates all_hits.tbl line for line. The
  chunks are what rule 101 exists for. all_proteins.faa is the 1,3 million-protein
  concatenation flagged for deletion several volumes ago and still there. Deleting
  the two would return about 840 MB and remove the double-count trap at source.
  Not this script's business, but it costs nothing to say while looking at it.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
