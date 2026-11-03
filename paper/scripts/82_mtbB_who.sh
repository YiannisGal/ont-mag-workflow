#!/usr/bin/env bash
# 82_mtbB_who.sh -- which organisms hold mtbB's reads, and does the ceiling
# even apply to them?
#
# MY EXPLANATION FOR mtbB WAS WRONG AND STEP 81 REFUTED IT. I predicted its
# 283 reads were MttB-superfamily sequence landing on mtbB references by
# accident. The runner-up test says otherwise:
#
#     label         reads   no other   gap<10%
#     mtmB             18       100%        0%     <- the calibrated marker
#     DECOYmtbB       282       100%        0%     <- indistinguishable
#
# Those reads see nothing outside the mtbB label. They are not crossovers.
#
# WHAT IS ACTUALLY ODD. The excess sits in two or three individual references:
# TQS84672.1 holds 110 of 282 sediment reads (19.9x its expected share) and
# AAM03976.1 holds 129 of 296 culture reads (22.2x). Rule 109 says uneven
# distribution is promiscuity only when the references are INTERCHANGEABLE --
# and these reads have no outside hits, so the references are not
# interchangeable and the unevenness is not promiscuity. Something real is
# being detected; the question is what.
#
# THE ASSUMPTION THE CEILING RESTS ON, WHICH NOBODY HAS CHECKED. "mtbB cannot
# exceed mcrA" is true only if every mtbB carrier is a methanogen. It is not a
# law about the gene, it is an inference from who has it. MttB-superfamily
# methyltransferases occur in BACTERIA -- this project already found that the
# hard way when an NCBI protein-name query imported seven bacterial proteins
# into the mtmB target set, two of them Desulfosporosinus meridiei, the genus
# of this sediment's own unbinned sulfate reducer. A bacterium carrying an
# mtbB-like gene has no mcrA, so its reads are not bounded by mcrA at all.
#
# Step 63 removed those seven by NAME. A name query finds what the database
# calls a name (rule 107); it does not find a bacterial protein filed under a
# different description. So the mtbB label may still hold bacteria, and this
# step asks the only question that settles it: WHO are these references?
#
# No network needed -- the organism is in the fasta header.
set -u
WS2="${WS2}"
P61="$WS2/61_mtmB"
P63="$WS2/63_mtmB_relabel"
OUT="$WS2/82_mtbB_who"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

"$PY" - "$P61" "$P63" "$OUT" <<'PYEOF'
import sys, os, re, collections
P61, P63, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MINLEN = 30
SHOW = ("mtmB", "DECOYmtaB", "DECOYmtbB", "mcrA")

bp = os.path.join(P63, "bacterial.txt")
bact = {l.strip() for l in open(bp)} if os.path.exists(bp) else set()
if not bact:
    sys.exit("   bacterial.txt missing -- refusing")


def lab(i):
    return "DECOYbact" if i.split("|")[-1] in bact else i.split("|")[0]


# --- full headers, so the organism travels with the id
hdr = {}
for f in sorted(os.listdir(P61)):
    if not f.endswith((".faa", ".fa", ".fasta")):
        continue
    for ln in open(os.path.join(P61, f), errors="replace"):
        if ln.startswith(">"):
            i = ln[1:].split()[0]
            hdr.setdefault(i, ln[1:].rstrip("\n"))

# organism: NCBI puts it in [brackets]; fall back to the trailing words
def organism(h):
    m = re.search(r"\[([^\]]+)\]\s*$", h)
    if m:
        return m.group(1)
    m = re.search(r"OS=([A-Za-z0-9 ._-]+?)(?: [A-Z]{2}=|$)", h)
    if m:
        return m.group(1).strip()
    return "(no organism in header)"


# A methanogen is what the ceiling assumes. Every known methanogen genus is
# Methano*, plus a short list that is not. Anything else is NOT covered by the
# ceiling -- and is printed so the call can be checked by eye rather than
# trusted.
METHANOGEN = re.compile(
    r"^(Methano|Methermicoccus|Methanomassiliicoccus|Candidatus Methano)", re.I)


def best(tag):
    top = {}
    for ln in open(os.path.join(P61, "%s_hits.tsv" % tag), errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        try:
            L, s = int(f[3]), float(f[5])
        except ValueError:
            continue
        if L < MINLEN:
            continue
        c = top.get(f[0])
        if c is None or s > c[1]:
            top[f[0]] = (f[1], s)
    return top


out = open(os.path.join(OUT, "reference_taxonomy.tsv"), "w")
out.write("sample\tlabel\treference\treads\torganism\tmethanogen\n")

for tag in ("SED", "LIQ"):
    top = best(tag)
    byref = collections.Counter(r for r, _ in top.values())
    bylab = collections.defaultdict(collections.Counter)
    for r, n in byref.items():
        bylab[lab(r)][r] = n
    mcra = sum(bylab.get("mcrA", {}).values())

    print("=" * 68)
    print("%s   (mcrA = %d reads)" % (tag, mcra))
    print("=" * 68)
    for L in SHOW:
        cnt = bylab.get(L)
        if not cnt:
            continue
        tot = sum(cnt.values())
        meth = non = 0
        print("   %s -- %d reads across %d references" % (L, tot, len(cnt)))
        print("      %5s  %-22s %s" % ("reads", "reference", "organism"))
        for r, c in cnt.most_common(12):
            h = hdr.get(r, "")
            org = organism(h)
            ism = bool(METHANOGEN.match(org))
            if ism:
                meth += c
            else:
                non += c
            print("      %5d  %-22s %s%s"
                  % (c, r.split("|")[-1][:22], org[:44],
                     "" if ism else "   <-- NOT a methanogen"))
            out.write("%s\t%s\t%s\t%d\t%s\t%s\n"
                      % (tag, L, r, c, org, "yes" if ism else "no"))
        # the tail, counted but not printed
        for r, c in cnt.most_common()[12:]:
            h = hdr.get(r, "")
            org = organism(h)
            ism = bool(METHANOGEN.match(org))
            meth += c if ism else 0
            non += 0 if ism else c
            out.write("%s\t%s\t%s\t%d\t%s\t%s\n"
                      % (tag, L, r, c, org, "yes" if ism else "no"))
        if len(cnt) > 12:
            print("      ... %d further reference(s), counted below"
                  % (len(cnt) - 12))
        print("      methanogen references hold %d reads; others hold %d"
              % (meth, non))
        if mcra and L != "mcrA":
            print("      ratio to mcrA: all %.2f   methanogen-only %.2f"
                  % (tot / float(mcra), meth / float(mcra)))
        print()
out.close()
print("   wrote %s" % os.path.join(OUT, "reference_taxonomy.tsv"))
print()
print("=" * 68)
print("WHAT DECIDES IT")
print("=" * 68)
print("   The ceiling is not a property of the gene. It is the inference that")
print("   every carrier also carries mcrA, which holds for methanogens and for")
print("   nobody else. So:")
print()
print("   - If the references holding mtbB's excess are METHANOGENS, the count")
print("     genuinely exceeds what the catalogue can explain and stays")
print("     unreportable -- the reads come from methanogen lineages with no")
print("     genome in this catalogue, which is a finding about coverage.")
print()
print("   - If they are NOT methanogens, the ceiling never applied to them.")
print("     The methanogen-only ratio above is then the number to judge, and")
print("     step 63's bacterial removal was incomplete because it searched by")
print("     protein NAME rather than by organism (rule 107).")
PYEOF
