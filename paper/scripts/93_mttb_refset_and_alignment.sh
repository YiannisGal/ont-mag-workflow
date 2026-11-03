#!/usr/bin/env bash
# 93_mttb_refset_and_alignment.sh -- the same audit for mttB, and the question
# of whether the reference sets are what makes these alignments short.
#
# WHY mttB IS NOW SUSPECT. Its verdict table says
#
#   every one of the five queries      CONFIRMED, k_nearest = TMA_MttB:10
#   closest reference, all five        TMA_MttB|CAQ5842095.1
#   distance to it                     1,2171 to 1,2830
#
# Three things are wrong with that picture at once. A closest relative at more
# than 1,2 substitutions per site is not a close relative. The SAME accession
# being closest to all five queries is what happens when nothing in the panel
# is near any of them and one sequence sits marginally nearer the middle. And
# step 92 printed the panels: the TMA_MttB set is mostly BACTERIAL
# (Anaerovoracaceae, Desulfosporosinus meridiei x10, Desulfobacula,
# Aquamicrobium, Magnetomoraceae, Bilophila) with only a handful of archaeal
# entries, while the BETAINE_SF "decoy" set is 39 "MttB family protein"
# sequences from Alphaproteobacteria plus one genuine glycine betaine
# methyltransferase. Our five queries are archaeal.
#
# So the contrast the mttB tree is supposed to draw -- true trimethylamine
# methyltransferase versus the betaine superfamily -- may be being drawn
# between two sets of bacterial proteins, neither of which is near the queries.
#
# THE CONTROL THAT DECIDES IT, AND IT IS NOT AN OPINION. For each panel,
# measure how far every REFERENCE is from its own nearest neighbour in the same
# panel. That gives the scale the alignment works on. If the references sit
# 1,2 apart from each other too, then our queries at 1,25 are unremarkable and
# the verdict stands. If the references sit 0,3 apart and only our queries are
# at 1,25, then the queries are outliers to the whole reference set and
# "CONFIRMED" means only "closer to this panel than to that one", which is a
# much weaker statement than the report currently makes.
#
# Section D asks the same of mcrA, where the answer should be the opposite --
# our sequences at 0,014-0,097 should sit INSIDE the reference spacing.
#
# Section E asks where the 139 sites went, by counting columns before and after
# trimming and measuring per-column occupancy on the untrimmed alignment.
set -u
WS2="${WS2}"
OUT="$WS2/93_mttb_audit"
mkdir -p "$OUT" || exit 1
PY=$(command -v python3 || command -v python) || { echo "no python"; exit 1; }

"$PY" - "$WS2" "$OUT" <<'PYEOF'
import sys, os, re, collections, statistics
WS2, OUT = sys.argv[1], sys.argv[2]
M = os.path.join(WS2, "15_marker_trees")

ARCH = ("methano", "methermicoccus", "archaeoglob", "halo", "thermoplasma",
        "candidatus methano", "nitrososphaer", "sulfolobus", "thermococc",
        "pyrococc", "bathyarchae", "asgard", "lokiarchae")


def is_arch(org):
    o = (org or "").lower()
    return any(o.startswith(a) or a in o for a in ARCH)


def deflines(d):
    out = {}
    rd = os.path.join(d, "refs")
    if not os.path.isdir(rd):
        return out
    for fn in sorted(os.listdir(rd)):
        if not fn.endswith((".faa", ".fa", ".fasta")):
            continue
        for ln in open(os.path.join(rd, fn), errors="replace"):
            if not ln.startswith(">"):
                continue
            h = ln[1:].rstrip("\n")
            a = h.split()[0]
            # sp|Q24SP7.1|MTGB_DESHY style ids keep their whole first token
            out[a] = h[len(a):].strip()
            if "|" in a:
                out[a.split("|")[-1]] = out[a]
                out[a.split("|")[1] if a.count("|") > 1 else a] = out[a]
    return out


def org_of(desc):
    m = re.search(r"\[([^\]]+)\]", desc or "")
    if m:
        return m.group(1)
    m = re.search(r"_([A-Z]{4,6})\b", desc or "")
    return ""


def mldist(path):
    lines = [l.rstrip("\n") for l in open(path, errors="replace") if l.strip()]
    names, M_ = [], []
    for l in lines[1:]:
        f = l.split()
        if len(f) < 2:
            continue
        names.append(f[0]); M_.append([float(x) for x in f[1:]])
    return names, M_


def panel(t):
    return t.split("|")[0] if "|" in t else "?"


def accn(t):
    return t.split("|")[-1]


for marker in ("mttB", "mcrA"):
    d = os.path.join(M, marker)
    md = os.path.join(d, "%s.mldist" % marker)
    print("=" * 62)
    print("%s" % marker)
    print("=" * 62)
    if not os.path.exists(md):
        print("   MISSING %s" % md); continue
    DESC = deflines(d)
    names, MM = mldist(md)
    print("   %d sequence(s) in the distance matrix, %d defline(s) on disk"
          % (len(names), len(DESC)))

    # ---- composition by domain, per panel
    pan = collections.defaultdict(list)
    for t in names:
        pan[panel(t)].append(t)
    for pname in sorted(pan):
        if pname == "QUERY":
            continue
        arc = bac = unk = 0
        for t in pan[pname]:
            dsc = DESC.get(accn(t)) or DESC.get(t) or ""
            o = org_of(dsc)
            if not o:
                unk += 1
            elif is_arch(o):
                arc += 1
            else:
                bac += 1
        print("   panel %-12s n=%3d   archaeal %3d   non-archaeal %3d   unknown %3d"
              % (pname, len(pan[pname]), arc, bac, unk))

    # ---- THE CONTROL: how far is each reference from its own nearest neighbour
    idx = {n: i for i, n in enumerate(names)}
    print()
    print("   how far each REFERENCE sits from its nearest neighbour in the")
    print("   same panel -- this is the scale the alignment works on:")
    for pname in sorted(pan):
        if pname == "QUERY" or len(pan[pname]) < 2:
            continue
        nn = []
        for t in pan[pname]:
            i = idx[t]
            v = [MM[i][idx[u]] for u in pan[pname] if u != t]
            if v:
                nn.append(min(v))
        print("      %-12s min %.4f  median %.4f  max %.4f"
              % (pname, min(nn), statistics.median(nn), max(nn)))

    qs = [t for t in names if panel(t) == "QUERY"]
    if not qs:
        print("   no QUERY sequences"); print(); continue
    print()
    print("   and our sequences, against each panel:")
    print("   %-24s %-14s %10s   %s" % ("query", "panel", "nearest", "which"))
    flag = []
    for q in qs:
        i = idx[q]
        f = q.split("|")
        lab = (f[1] + " " + f[3]) if len(f) > 3 else q
        for pname in sorted(pan):
            if pname == "QUERY":
                continue
            c = sorted((MM[i][idx[u]], u) for u in pan[pname])
            dsc = DESC.get(accn(c[0][1])) or ""
            o = org_of(dsc) or accn(c[0][1])
            print("   %-24s %-14s %10.4f   %s" % (lab.replace("ETO_", "")[:24],
                                                  pname, c[0][0], o[:38]))
        # the same thing restricted to ARCHAEAL references of any panel
        arcs = [u for u in names if panel(u) != "QUERY"
                and is_arch(org_of(DESC.get(accn(u)) or ""))]
        if arcs:
            c = sorted((MM[i][idx[u]], u) for u in arcs)[0]
            dsc = DESC.get(accn(c[1])) or ""
            print("   %-24s %-14s %10.4f   %s"
                  % ("", "ARCHAEAL refs", c[0], org_of(dsc)[:38]))
            flag.append(c[0])
    print()
    if flag:
        print("   nearest ARCHAEAL reference, across all queries: %.4f to %.4f"
              % (min(flag), max(flag)))
    print()
PYEOF

echo
echo "=============================================================="
echo "E. where the alignment columns went"
echo "=============================================================="
echo "   139 sites for 105 sequences was called a thin alignment. Whether it is"
echo "   a trimming choice or a consequence of aligning proteins that are not"
echo "   homologous is decided by counting columns before and after, and by"
echo "   measuring how well each column is filled BEFORE trimming."
echo
"$PY" - "$WS2" <<'PYEOF'
import sys, os, collections
WS2 = sys.argv[1]
M = os.path.join(WS2, "15_marker_trees")


def read_fa(p):
    seqs, name, buf = {}, None, []
    for ln in open(p, errors="replace"):
        if ln.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = ln[1:].strip(), []
        else:
            buf.append(ln.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


for marker in ("mcrA", "mttB"):
    d = os.path.join(M, marker)
    aln = os.path.join(d, "%s_aln.faa" % marker)
    trm = os.path.join(d, "%s_trim.faa" % marker)
    print("   --- %s ---" % marker)
    if not os.path.exists(aln):
        print("      MISSING %s" % aln); continue
    A = read_fa(aln)
    L = len(next(iter(A.values())))
    print("      untrimmed: %d sequence(s) x %d columns" % (len(A), L))
    if os.path.exists(trm):
        Tr = read_fa(trm)
        LT = len(next(iter(Tr.values())))
        print("      trimmed:   %d sequence(s) x %d columns  (%.1f%% of columns kept)"
              % (len(Tr), LT, 100.0 * LT / L))
    # occupancy per column, and how it splits by panel
    def occ(keys):
        if not keys:
            return []
        return [sum(1 for k in keys if A[k][j] != "-") / float(len(keys))
                for j in range(L)]
    groups = collections.defaultdict(list)
    for k in A:
        groups[k.split("|")[0]].append(k)
    allo = occ(list(A))
    print("      columns filled by >=90%% of ALL sequences: %d"
          % sum(1 for x in allo if x >= 0.9))
    for g in sorted(groups):
        o = occ(groups[g])
        print("      %-12s n=%3d   columns filled by >=90%% of this group: %d"
              % (g, len(groups[g]), sum(1 for x in o if x >= 0.9)))
    # what happens if the known-bad sequences are dropped
    bad = [k for k in A if "WP_484489188" in k or "WP_471601478" in k]
    if bad:
        keep = [k for k in A if k not in bad]
        o = occ(keep)
        print("      dropping the %d mis-imported sequence(s): columns filled by"
              " >=90%% would be %d" % (len(bad), sum(1 for x in o if x >= 0.9)))
    print()
PYEOF

echo
echo "=============================================================="
echo "WHAT THIS DECIDES"
echo "=============================================================="
echo "   If the mttB references sit ~1,2 apart from each other, our queries at"
echo "   1,22-1,28 are normal and fig15 can be drawn like fig14. If they sit"
echo "   much closer, the mttB reference set does not contain the queries'"
echo "   relatives, the CONFIRMED verdicts mean only 'nearer this panel than"
echo "   that one', and FINDINGS_REPORT section 5's MttB sentence has to say so."
echo "   E says whether rebuilding the mcrA reference set would buy back"
echo "   alignment columns, which is the only reason to rerun that tree."
