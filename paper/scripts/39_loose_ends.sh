#!/usr/bin/env bash
# =============================================================================
# 39_loose_ends.sh  --  three open items that need measurement, not new tools
#
#   A. `all_mags_contigs.fna` header provenance   (open since step 18)
#   B. archaeal secretome uncertainty, MEASURED   (open since step 23)
#   C. culturomics workspace inventory            (never opened)
#
# Each is a question this project has been carrying as "unexplained" or
# "approximate". None needs a download or an install; all three need somebody
# to look, and to write down what is actually there rather than what is assumed.
#
# ---------------------------------------------------------------------------
# A. all_mags_contigs.fna — WHY DO ITS HEADERS NOT RESOLVE?
#
#   Volume III: "unexplained; nothing from that file may go in the thesis."
#   geNomad called 38 proviruses on it and assigned ZERO to a MAG; 254 contig
#   names were discarded as ambiguous.
#
#   This does not guess. It compares that file's headers against the real MAG
#   contig headers and tries to RECOVER the transformation: exact match, case
#   folding, pipe-prefix stripping, trailing-index stripping, and the substring
#   relation in both directions. Whatever fraction each rule explains is
#   reported. If one rule explains >=95%, the provenance is solved and step 18's
#   provirus calls become usable; if none does, the file stays excluded — but
#   with a measured reason, and the ASSEMBLY hypothesis is then tested directly
#   against 04_polish/ and 03_assembly/.
#
#   DETERMINISM: candidate matches are searched over a SORTED list, never over
#   a set. Iterating a set of strings gives a different first match between runs
#   under a different hash seed — the defect found in 16b's genus column.
#
# ---------------------------------------------------------------------------
# B. ARCHAEAL SECRETOME — measuring the size of a known approximation
#
#   DeepSig has no archaeal model, so all 8 archaeal genomes were scored with
#   the GRAM-POSITIVE one and the notebook calls those percentages "approximate
#   and currently uncorrectable". Two things are wrong with leaving it there.
#   "Approximate" is not a number, and "uncorrectable" was never tested.
#
#   Archaeal Sec signal peptides genuinely resemble Gram-positive ones — neither
#   cell type has an outer membrane, both use long weakly-charged n-regions —
#   which is WHY Gram-positive was chosen. So the question is not "is it wrong"
#   but "by how much". That is measurable from what is already on disk: the
#   three sequence properties every signal-peptide predictor keys on, computed
#   for archaeal vs bacterial proteomes in THIS dataset.
#
#     n-region net charge (residues 1-8)
#     h-region hydrophobicity (max Kyte-Doolittle over a full 9-residue window)
#     AXA-like cleavage motif in the 15-35 window
#
#   CONTROL — metric sensitivity. A difference of zero is only meaningful if the
#   metric could have detected one. Each N-terminus is re-measured after a
#   composition-preserving shuffle (fixed seed). If real N-termini do not
#   separate from their own shuffles, the metric has no resolution and the
#   comparison is reported as UNINFORMATIVE rather than as agreement.
#
#   LIMIT, stated in the output: this measures the INPUT distribution, not
#   predictor accuracy. It cannot validate a per-protein call. It replaces an
#   unquantified caveat with a bounded one.
#
# ---------------------------------------------------------------------------
# C. CULTUROMICS WORKSPACE — an inventory, not an analysis
#
#   `culturomics_workspace.txt` has been "still unopened" since Volume II. This
#   locates it, classifies what it describes, and reports the two things the
#   metagenome arm has been waiting for: isolate genomes that could be ANI-
#   matched to the 29 MAGs (step 12), and the Arc_Otu22 representative sequence
#   step 09 asked for by name and never received (audit finding F23). It draws
#   no conclusions; the point is to replace "unopened" with a list.
#
# OUTPUT  39_loose_ends/{PROVENANCE,SECRETOME,CULTUROMICS}.txt · LOOSE_ENDS.txt
# RUNTIME 2-10 min
# USAGE   ws2 && bash scripts/39_loose_ends.sh 2>&1 | tee logs/39.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="39_loose_ends"
A="$WS2/09_annotation"; F="$WS2/07_mags"
D="$WS2/39_loose_ends"; mkdir -p "$D"

echo "######## 39 — THREE OPEN ITEMS ########"
echo "  started: $(date -Is)"
# delete first: a crashed run must never be validated against the last run's files
rm -f "$D/LOOSE_ENDS.txt" "$D/PROVENANCE.txt" "$D/SECRETOME.txt" "$D/CULTUROMICS.txt"

# the two files that are not where any config points
AMC="$(find "$WS2" -maxdepth 3 -name 'all_mags_contigs.fna' 2>/dev/null | sort | head -1)"
# -path exclusion matters: this step's OWN output directory is inside $WS2 and
# CULTUROMICS.txt matches *culturomic*. Without it the script finds itself.
CWS="$(find "$WS2" "$HOME" "$(dirname "$WS2")" -maxdepth 4 -type f \
       -name 'culturomics_workspace*' -not -path "$D/*" 2>/dev/null | sort | head -1)"
[ -z "$CWS" ] && CWS="$(find "$WS2" "$HOME" "$(dirname "$WS2")" -maxdepth 4 -type f \
       -iname '*culturomic*' -not -path "$D/*" 2>/dev/null | sort | head -1)"
echo "  all_mags_contigs.fna : ${AMC:-<not found>}"
echo "  culturomics workspace: ${CWS:-<not found>}"

echo
python3 - "$D" "$A" "$F" "$WS2" "${AMC:-}" "${CWS:-}" <<'PY'
import sys, os, re, glob, csv, random, traceback
from collections import Counter, defaultdict

D, A, F, WS2, AMC, CWS = sys.argv[1:7]
W = 100
rep = []
STATUS = ['INCOMPLETE']
FAILED = []


def say(s=''):
    print(s)
    rep.append(s)


def rule(t=''):
    say('=' * W)
    if t:
        say('  ' + t)
        say('=' * W)


def dump():
    """status is written ONCE, last, only from here. Never before completion."""
    with open(os.path.join(D, 'LOOSE_ENDS.txt'), 'w') as fh:
        fh.write('\n'.join(rep) + '\n')
        fh.write(f"STEP_STATUS: {STATUS[0]}\n")
    print(f"\n  written: {os.path.join(D, 'LOOSE_ENDS.txt')}")


def section(letter, title, fn, outfile):
    """Contain a failure to its own section: a crash in B must not erase A."""
    mark = len(rep)
    rule(title)
    try:
        fn()
    except Exception as e:
        FAILED.append(letter)
        say()
        say(f"  !! SECTION {letter} FAILED: {type(e).__name__}: {e}")
        for ln in traceback.format_exc().strip().split('\n')[-4:]:
            say(f"     {ln}")
    with open(os.path.join(D, outfile), 'w') as fh:
        fh.write('\n'.join(rep[mark:]) + '\n')


def heads(path, cap=200000):
    out = []
    try:
        with open(path, errors='replace') as fh:
            for line in fh:
                if line.startswith('>'):
                    out.append(line[1:].rstrip('\n'))
                    if len(out) >= cap:
                        break
    except Exception as e:
        say(f"  !! cannot read {path}: {e}")
    return out


def mag_fastas():
    """Per-MAG FASTAs only. The query file MUST NOT enter its own reference set:
    all_mags_contigs.fna lives in 07_mags/ and a bare *.fna glob swallows it,
    which would make every rule match the file against itself and report SOLVED
    no matter what. Aggregate files are excluded by name as well as by path."""
    got = []
    for sub in ('bins', 'final', 'MAGs', ''):
        hits = []
        for e in ('fa', 'fna', 'fasta'):
            hits += glob.glob(os.path.join(F, sub, '*.' + e) if sub
                              else os.path.join(F, '*.' + e))
        if hits:
            got = hits
            break                      # first populated location wins; do not mix
    bad = os.path.realpath(AMC) if AMC else None
    out = []
    for p in got:
        b = os.path.basename(p)
        if bad and os.path.realpath(p) == bad:
            continue
        if b.startswith('all_') or 'all_mags' in b or b.startswith('concat'):
            continue
        out.append(p)
    return sorted(set(out))


# =====================================================================  A
def sec_a():
    if not AMC or not os.path.exists(AMC):
        say("  all_mags_contigs.fna not found under $WS2 (maxdepth 3).")
        say("  Recorded as NOT LOCATED — not as resolved, and not as absent.")
        return
    q = heads(AMC)
    say(f"  {os.path.relpath(AMC, WS2)}: {len(q):,} headers")
    for h in q[:3]:
        say(f"    > {h[:88]}")
    say()

    # the truth: contig names inside the frozen MAG FASTAs
    ref = {}
    fas = mag_fastas()
    for fa in fas:
        mag = re.sub(r'\.(fa|fna|fasta)$', '', os.path.basename(fa))
        for h in heads(fa):
            ref[h.split()[0] if h.split() else h] = mag
    say(f"  MAG contig headers on record: {len(ref):,} from {len(fas)} FASTA(s)")
    if not ref:
        say(f"  !! no MAG FASTAs found under {os.path.relpath(F, WS2)} —")
        say("     cannot compare. The provenance question stays open for want of")
        say("     the reference, not for want of a test.")
        return

    # ---- deterministic lookup structures: sorted lists, never sets
    R = set(ref)
    Rsorted = sorted(R)
    Rlow = {}
    for r in Rsorted:                       # first in sorted order wins
        Rlow.setdefault(r.lower(), r)

    def strip_pref(x):
        return re.sub(r'^[^|]*\|', '', x)

    def strip_suf(x):
        return re.sub(r'[_.]\d+$', '', x)

    def first_field(x):
        p = x.split()
        return p[0] if p else x

    def sub_qr(x):                           # query inside a MAG contig name
        return next((r for r in Rsorted if x and x in r), None)

    def sub_rq(x):                           # MAG contig name inside query
        return next((r for r in Rsorted if r and r in x), None)

    # (name, fn, is_quadratic)
    rules = [
        ('exact match',                        lambda x: x if x in R else None, 0),
        ('first whitespace field',             lambda x: first_field(x) if first_field(x) in R else None, 0),
        ('case-folded',                        lambda x: Rlow.get(x.lower()), 0),
        ('drop text before a pipe',            lambda x: strip_pref(x) if strip_pref(x) in R else None, 0),
        ('drop a trailing _N / .N',            lambda x: strip_suf(x) if strip_suf(x) in R else None, 0),
        ('query is a SUBSTRING of a contig',   sub_qr, 1),
        ('a contig is a substring of query',   sub_rq, 1),
    ]

    fast_n = min(len(q), 4000)
    slow_n = min(len(q), 300)                # substring rules are O(n*m)
    say()
    say(f"  {'transformation':<40}{'n tested':>10}{'resolved':>10}{'%':>10}")
    say('-' * W)
    best = (0.0, None, 0)
    for name, fn, quad in rules:
        probe = q[:slow_n] if quad else q[:fast_n]
        n = 0
        for x in probe:
            try:
                if fn(x):
                    n += 1
            except Exception:
                pass
        pct = 100.0 * n / max(1, len(probe))
        say(f"  {name:<40}{len(probe):>10,}{n:>10,}{pct:>9.1f}%")
        if pct > best[0]:
            best = (pct, name, len(probe))
    say()
    say(f"  (non-quadratic rules on the first {fast_n:,} headers; the two")
    say(f"   substring rules on {slow_n:,}, since each is O(n x {len(Rsorted):,}).")
    say("   The substring rules are AMBIGUOUS by construction — contig_1 is inside")
    say("   contig_11 — so a high score from one of them is a lead, not a mapping.)")
    say()

    if best[0] >= 95.0:
        say(f"  >> SOLVED: '{best[1]}' explains {best[0]:.1f}% of the headers")
        say(f"     (n = {best[2]:,}).")
        say("     The file's contigs CAN be mapped back to MAGs, so step 18's 38")
        say("     provirus calls are recoverable and the '0 assigned to a MAG'")
        say("     result was a name-space failure, not a biological absence.")
        say("     NEXT: re-run the step-18 assignment through this rule before")
        say("     anything from that file enters the thesis.")
        return

    if best[0] > 0:
        say(f"  >> PARTIAL: the best rule ('{best[1]}') explains only {best[0]:.1f}%.")
        say("     Not enough to trust — a partial name-space mapping would assign")
        say("     some proviruses and silently drop the rest, which is worse than")
        say("     assigning none. Headers matching nothing:")
        unres = [x for x in q[:fast_n] if x not in R and first_field(x) not in R][:5]
        for x in unres:
            say(f"       {x[:84]}")
    else:
        say("  >> NO RULE MATCHES A SINGLE HEADER. The file does not contain MAG")
        say("     contigs under any name this test can construct.")

    # the alternative hypothesis, tested rather than asserted
    asm = []
    for pat in ('04_polish/*.fasta', '04_polish/*.fa',
                '03_assembly/*/assembly.fasta', '03_assembly/*/*/assembly.fasta'):
        asm += glob.glob(os.path.join(WS2, pat))
    asm = sorted(set(asm))
    say()
    if not asm:
        say("  assembly hypothesis: NOT TESTED — no assembly FASTA found under")
        say("  04_polish/ or 03_assembly/. Absence of the test, not of the match.")
    else:
        say(f"  ASSEMBLY HYPOTHESIS — were these headers taken from the assemblies")
        say(f"  rather than from the MAG FASTAs? Testing {len(asm)} file(s):")
        probe = q[:fast_n]
        union = set()
        for a2 in asm[:4]:
            ah = {(h.split()[0] if h.split() else h) for h in heads(a2, 200000)}
            union |= ah
            n = sum(1 for x in probe if first_field(x) in ah)
            say(f"    {os.path.relpath(a2, WS2)[:58]:<60}"
                f"{100.0 * n / max(1, len(probe)):>8.1f}% match")
        # the file plausibly spans BOTH assemblies, so the union is the
        # hypothesis, not the best single file
        nu = sum(1 for x in probe if first_field(x) in union)
        hit_any = 100.0 * nu / max(1, len(probe))
        say(f"    {'UNION of all assemblies':<60}{hit_any:>8.1f}% match")
        say()
        if hit_any >= 95.0:
            say(f"  >> CONFIRMED at {hit_any:.1f}%: the file was built from the")
            say("     ASSEMBLIES, not from the MAG FASTAs. That fully explains step")
            say("     18's zero — an assembly contig has no MAG to be assigned to")
            say("     unless it was binned, and 254 were not. The 38 proviruses are")
            say("     real calls on real contigs; they are simply community-level,")
            say("     and may be reported as such, per assembly, not per MAG.")
        elif hit_any > 0:
            say(f"  >> the assemblies explain {hit_any:.1f}% — also not enough.")
            say("     The file's provenance remains unexplained and nothing from it")
            say("     may enter the thesis.")
        else:
            say("  >> the assemblies explain none of it either. The file came from")
            say("     somewhere this workspace does not contain. It stays excluded.")


# =====================================================================  B
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
POS = set('KR')
NEG = set('DE')


def nterm_stats(s):
    """the three properties every signal-peptide predictor keys on.
    Windows are FULL LENGTH only — a short tail must not be averaged over 9."""
    if len(s) < 25:
        return None
    n = s[:8]
    charge = sum(1 for c in n if c in POS) - sum(1 for c in n if c in NEG)
    win = [sum(KD.get(c, 0.0) for c in s[i:i + 9]) / 9.0
           for i in range(2, min(23, len(s) - 8))]
    if not win:
        return None
    hmax = max(win)
    axa = 1 if re.search(r'A.A', s[15:35]) else 0
    return charge, hmax, axa


def load_domains():
    dom = {}
    tp = os.path.join(F, 'MAG_TAXONOMY.tsv')
    if not os.path.exists(tp):
        return dom
    with open(tp, errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            blob = ' '.join(str(v or '') for v in r.values())
            d = next((str(r[c]).strip() for c in sorted(k for k in r if k)
                      if c.strip().lower() == 'domain'), '')
            if not d:
                m = re.search(r'd__([A-Za-z]+)', blob)
                d = m.group(1) if m else ''
            if not d:
                continue
            # a genome id in this project has an underscore AND a digit;
            # a genus name has neither. Row-local, one hop, never transitive.
            for v in r.values():
                v = str(v or '').strip()
                if v and '_' in v and any(ch.isdigit() for ch in v):
                    dom.setdefault(v, d)
    return dom


def sec_b():
    say("  DeepSig has no archaeal model, so the archaeal genomes were scored")
    say("  with the GRAM-POSITIVE one and the result called 'approximate and")
    say("  uncorrectable'. Archaeal Sec signal peptides do resemble Gram-positive")
    say("  ones — neither cell type has an outer membrane — which is why that")
    say("  model was chosen. The question is how far the N-termini in THIS")
    say("  dataset sit from that assumption.")
    say()
    dom = load_domains()
    if not dom:
        say("  !! no MAG_TAXONOMY.tsv / no domain column — cannot split by domain.")
        return

    rng = random.Random(20260904)     # fixed seed: the shuffle control is reproducible
    acc = defaultdict(lambda: {'n': 0, 'ch': 0.0, 'kd': 0.0, 'ax': 0,
                               'ch2': 0.0, 'kd2': 0.0, 'sk': 0.0, 'g': set()})

    def add(grp, gid, st, shuf_kd):
        p = acc[grp]
        p['n'] += 1
        p['ch'] += st[0]
        p['kd'] += st[1]
        p['ax'] += st[2]
        p['ch2'] += st[0] ** 2
        p['kd2'] += st[1] ** 2
        p['sk'] += shuf_kd
        p['g'].add(gid)

    for m in (sorted(os.listdir(A)) if os.path.isdir(A) else []):
        faa = os.path.join(A, m, m + '.faa')
        if not os.path.exists(faa):
            continue
        d = dom.get(m, '')
        if not d:
            continue
        grp = 'Archaea' if d.lower().startswith('archae') else 'Bacteria'
        seq, have = [], False

        def flush():
            if not have:
                return
            s = ''.join(seq)[:40]
            st = nterm_stats(s)
            if not st:
                return
            ch = list(s)
            rng.shuffle(ch)
            sh = nterm_stats(''.join(ch))
            add(grp, m, st, sh[1] if sh else st[1])

        with open(faa, errors='replace') as fh:
            for line in fh:
                if line.startswith('>'):
                    flush()
                    seq, have = [], True
                else:
                    seq.append(line.strip())
        flush()

    if not acc:
        say("  !! no proteomes matched a genome with a domain assignment.")
        return

    say(f"  {'domain':<11}{'MAGs':>6}{'proteins':>10}{'n-charge':>16}"
        f"{'h-region KD':>18}{'AXA':>9}{'shuffled KD':>14}")
    say('-' * W)
    for grp in ('Bacteria', 'Archaea'):
        p = acc.get(grp)
        if not p or not p['n']:
            continue
        n = p['n']
        sd_ch = max(0.0, p['ch2'] / n - (p['ch'] / n) ** 2) ** 0.5
        sd_kd = max(0.0, p['kd2'] / n - (p['kd'] / n) ** 2) ** 0.5
        say(f"  {grp:<11}{len(p['g']):>6}{n:>10,}"
            f"{p['ch'] / n:>10.2f}+-{sd_ch:<4.2f}"
            f"{p['kd'] / n:>12.3f}+-{sd_kd:<4.3f}"
            f"{100.0 * p['ax'] / n:>8.1f}%{p['sk'] / n:>14.3f}")
    say()

    b, a = acc.get('Bacteria'), acc.get('Archaea')

    # ---- CONTROL first. A null result is only meaningful from a metric with
    #      resolution. Real vs composition-preserving shuffle, same pipeline.
    sens = []
    for grp in ('Bacteria', 'Archaea'):
        p = acc.get(grp)
        if p and p['n']:
            sens.append(abs(p['kd'] / p['n'] - p['sk'] / p['n']))
    ctrl = min(sens) if sens else 0.0
    say("  CONTROL — metric sensitivity (real vs composition-shuffled N-termini,")
    say("  seed 20260904; the shuffle keeps amino-acid composition and destroys")
    say(f"  only POSITION): smallest real-vs-shuffle h-region gap = {ctrl:.3f} KD units")
    if ctrl < 0.05:
        say("  >> [FAIL] the measurement cannot separate a real N-terminus from its")
        say("     own shuffle, so it has no resolution and CANNOT bound the")
        say("     approximation. Reported as UNINFORMATIVE, not as agreement.")
        say("     The archaeal secretome caveat stands unquantified.")
        return
    say("  >> [PASS] the metric responds to N-terminal positional structure, so a")
    say("     small archaea-vs-bacteria difference is evidence of similarity")
    say("     rather than of a blunt instrument.")
    say()

    if not (b and a and b['n'] and a['n']):
        say("  only one domain is represented in the annotated proteomes —")
        say("  no cross-domain comparison possible. NOT MEASURED.")
        return

    dch = abs(b['ch'] / b['n'] - a['ch'] / a['n'])
    dkd = abs(b['kd'] / b['n'] - a['kd'] / a['n'])
    dax = abs(b['ax'] / b['n'] - a['ax'] / a['n']) * 100
    say(f"  archaea vs bacteria: charge {dch:.3f} | hydrophobicity {dkd:.3f}"
        f" | AXA {dax:.1f} points")
    say(f"  for scale, the control's real-vs-shuffle gap is {ctrl:.3f} KD units")
    say()
    if dkd < ctrl and dch < 0.25 and dax < 5:
        say("  >> The archaea-vs-bacteria difference is SMALLER than the metric's")
        say("     own real-vs-shuffle response, on a metric shown above to have")
        say("     resolution. On all three axes a signal-peptide model keys on,")
        say("     the archaeal N-termini in this dataset sit within the bacterial")
        say("     spread. The Gram-positive model is a SMALL approximation here.")
        say("     WRITE IT AS: 'archaeal secretome fractions were obtained with")
        say("     DeepSig's Gram-positive model; archaeal N-terminal charge,")
        say(f"     hydrophobicity and cleavage-motif frequency differ from bacterial")
        say(f"     by {dch:.2f}, {dkd:.2f} and {dax:.1f} points respectively' —")
        say("     instead of the word 'approximate'.")
    else:
        say("  >> The distributions DIFFER materially on at least one axis. The")
        say("     Gram-positive model is a MATERIAL approximation for these")
        say("     archaea. Report the archaeal secretome percentages with the")
        say("     three numbers above as the size of the caveat, and do not")
        say("     compare archaeal and bacterial secretome fractions to each")
        say("     other without it.")
    say()
    say("  LIMITS. This measures the INPUT distribution, not predictor accuracy;")
    say("  it cannot tell you any single call is right, and agreement of these")
    say("  three properties does not prove the model transfers. It replaces an")
    say("  unquantified caveat with a bounded one. The real remedy is SignalP 6.0,")
    say("  which has an archaeal mode and needs academic registration.")


# =====================================================================  C
def sec_c():
    if not CWS or not os.path.exists(CWS):
        say("  No file matching *culturomic* found under $WS2, $HOME, or the")
        say("  parent of $WS2 (maxdepth 4).")
        say("  Recorded as NOT LOCATED. It has been listed as 'still unopened'")
        say("  since Volume II; it is now also known not to be where it was")
        say("  expected. If it exists it is on another machine or under another")
        say("  name — that is a question for the culturomics arm, not a finding.")
        return
    lines = [l.rstrip('\n') for l in open(CWS, errors='replace')]
    shown = os.path.relpath(CWS, WS2) if CWS.startswith(WS2) else CWS
    say(f"  {shown}")
    say(f"  {len(lines):,} lines")
    say()
    ext = Counter()
    for l in lines:
        m = re.search(r'\.([A-Za-z0-9]{1,6})$', l.strip())
        if m:
            ext[m.group(1).lower()] += 1
    if ext:
        say("  file types described:")
        for e, n in sorted(ext.items(), key=lambda kv: (-kv[1], kv[0]))[:15]:
            say(f"    .{e:<12}{n:>7}")
        say()

    iso = [l for l in lines
           if re.search(r'\.(fa|fna|fasta|gbk|gbff|gff)$', l.strip(), re.I)]
    otu = [l for l in lines if re.search(r'otu|rep[_-]?seq|16s', l, re.I)]
    ani = [l for l in lines if re.search(r'\bani\b|fastani|skani|dddh|gtdb', l, re.I)]

    say(f"  genome / annotation files listed          : {len(iso):,}")
    for l in iso[:10]:
        say(f"    {l.strip()[:84]}")
    say()
    say(f"  lines mentioning OTU / 16S / rep sequences : {len(otu):,}")
    for l in otu[:10]:
        say(f"    {l.strip()[:84]}")
    say()
    say(f"  lines mentioning ANI / dDDH / GTDB         : {len(ani):,}")
    for l in ani[:6]:
        say(f"    {l.strip()[:84]}")
    say()
    say("  WHY THESE THREE. Step 12 flagged 'ANI: liquid MAGs vs culturomics")
    say("  isolates' as the link between the two chapters — that needs the")
    say("  isolate genome files. Step 09 asked for the Arc_Otu22 representative")
    say("  FASTA by name and never received it — that sequence is what closes")
    say("  audit finding F23 (16S identity reported as both 75.8% and 93.1%,")
    say("  two conflicting values on which a candidate-new-phylum claim rests).")
    say()
    if not iso:
        say("  >> NO genome files are listed. The ANI link cannot be made from")
        say("     this inventory alone; the isolate assemblies live elsewhere.")
    else:
        say(f"  >> {len(iso)} genome/annotation paths are named. The ANI comparison")
        say("     against the 29 MAGs is now a concrete, runnable job (skani or")
        say("     fastANI, both cheap) rather than an open item.")


# =====================================================================
section('A', 'A. all_mags_contigs.fna — RECOVERING THE HEADER TRANSFORMATION',
        sec_a, 'PROVENANCE.txt')
section('B', 'B. ARCHAEAL SECRETOME — HOW BIG IS THE APPROXIMATION?',
        sec_b, 'SECRETOME.txt')
section('C', 'C. CULTUROMICS WORKSPACE — WHAT IS ACTUALLY IN IT?',
        sec_c, 'CULTUROMICS.txt')

rule('WHAT TO REPORT')
say("  A. Either the header transformation is recovered — in which case step 18's")
say("     MAG-level provirus calls become usable and the '0 assigned' was a")
say("     name-space failure — or the assemblies explain it, in which case the 38")
say("     proviruses are community-level calls and may be reported per assembly.")
say("     If neither, the file stays excluded with a measured reason.")
say("  B. Report the archaeal secretome caveat with the three property")
say("     differences, NOT the word 'approximate' — but only if the sensitivity")
say("     control passed. If it did not, the caveat stays unquantified and that")
say("     is the honest statement.")
say("  C. The culturomics workspace is inventoried. Whether the isolate genomes")
say("     and the Arc_Otu22 sequence exist is now a fact rather than an")
say("     assumption, and the ANI job is either runnable or known not to be.")

if FAILED:
    STATUS[0] = 'PARTIAL (section ' + ','.join(FAILED) + ' failed)'
else:
    STATUS[0] = 'OK'
dump()
PY
RC=$?

echo
echo "######## VALIDATION ########"
OK=0
check "python block exited cleanly" "$RC" || OK=1
check "report written" "$([ -s "$D/LOOSE_ENDS.txt" ] && echo 0 || echo 1)" || OK=1
for f in PROVENANCE SECRETOME CULTUROMICS; do
  check "$f.txt written" "$([ -s "$D/$f.txt" ] && echo 0 || echo 1)" || OK=1
done
ST=$(grep -m1 '^STEP_STATUS:' "$D/LOOSE_ENDS.txt" 2>/dev/null | cut -d' ' -f2-)
echo "  status line: ${ST:-<none>}"
[ "${ST:-MISSING}" = "OK" ] || { echo "  [STOP] ${ST:-no status line}"; OK=1; }
[ "$OK" -eq 0 ] && runlog "39_loose_ends complete" || runlog "**39_loose_ends FAILED**"
echo "  ended: $(date -Is)"
exit "$OK"
