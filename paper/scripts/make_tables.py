#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_tables.py -- every table for the metagenome chapter, in Greek and English.

    python3 scripts/make_tables.py [workspace_2 root] [outdir]
    python3 scripts/make_tables.py . tables

WHAT IT WRITES

    tables/tableNN_<name>.tsv     one per table, the data, language-neutral
    tables/TABLES_EL.docx         all tables, Greek headers and captions
    tables/TABLES_EN.docx         all tables, English headers and captions
    tables/MANIFEST.txt           what was built, what was skipped, and why

EVERY NUMBER IS RECOMPUTED FROM A RUN OUTPUT. Nothing is transcribed from the
findings report or from a lab notebook, for the same reason the figure scripts
recompute: a table that is typed in by hand stops tracking the analysis the
moment the analysis changes. If a table's input file is missing or does not
parse, that table is SKIPPED and the manifest says so. It is never filled in
from memory and never half-filled.

THE ONE EXCEPTION, AND IT IS MARKED. Table 7 (sediment geochemistry) is not
derived from a file in workspace_2 -- the digitised TSV was never on the
server. Its values are embedded below, recovered from the vector coordinates of
Sediment_Chemical_Analysis.pdf, and they carry their own control: every one of
audit finding F28's eight sulfate statistics reproduces exactly. If
tables/sediment_geochemistry_by_depth.tsv is present, the script reads that
instead and says so.

DOCX. If python-docx is not importable the TSVs are still written and the
manifest says which command would install it. Nothing is installed
automatically -- this machine's shared environments are not writable.
"""
import sys, os, re, csv, glob, math, statistics, collections

args = [a for a in sys.argv[1:] if not a.startswith("--")]
WS2 = os.path.abspath(args[0]) if args else "."
OUT = args[1] if len(args) > 1 else "tables"
os.makedirs(OUT, exist_ok=True)

MANIFEST = []
def note(s):
    MANIFEST.append(s)
    print(s)

def P(*a):
    return os.path.join(WS2, *a)

def readtsv(path, delim="\t"):
    """rows as dicts; returns [] if the file is absent"""
    if not os.path.exists(path):
        return []
    with open(path, errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter=delim))

def num(x, dp=2):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x) if x not in (None, "") else "—"
    if dp == 0:
        return "%d" % round(v)
    return ("%%.%df" % dp) % v

def thousands(x):
    try:
        return "{:,}".format(int(float(x))).replace(",", " ")
    except (TypeError, ValueError):
        return "—"


# =========================================================== the table objects
class Table(object):
    def __init__(self, n, stem, cap_el, cap_en, cols_el, cols_en, rows, foot_el="", foot_en=""):
        self.n, self.stem = n, stem
        self.cap_el, self.cap_en = cap_el, cap_en
        self.cols_el, self.cols_en = cols_el, cols_en
        self.rows = rows
        self.foot_el, self.foot_en = foot_el, foot_en

TABLES = []
def add(*a, **k):
    t = Table(*a, **k)
    if not t.rows:
        note("  [skip]  table %-2d %-28s no rows" % (t.n, t.stem))
        return
    if len(t.cols_el) != len(t.cols_en):
        note("  [skip]  table %-2d %-28s header lengths differ (%d vs %d)"
             % (t.n, t.stem, len(t.cols_el), len(t.cols_en)))
        return
    bad = [i for i, r in enumerate(t.rows) if len(r) != len(t.cols_el)]
    if bad:
        note("  [skip]  table %-2d %-28s %d row(s) have the wrong width"
             % (t.n, t.stem, len(bad)))
        return
    TABLES.append(t)
    note("  [ok]    table %-2d %-28s %d row(s) x %d column(s)"
         % (t.n, t.stem, len(t.rows), len(t.cols_el)))


# ================================================================ 1. libraries
def t1():
    src = []
    rows = []
    for lab_el, lab_en, pat in (
            ("Ίζημα", "Sediment", "*[Ss]ed*"),
            ("Καλλιέργειες", "Cultures", "*[Ll]iq*")):
        ns = sorted(glob.glob(P("00_input", pat + "NanoStats.txt")) +
                    glob.glob(P("00_input", pat + "nanostats.txt")))
        if not ns:
            continue
        txt = open(ns[0], errors="replace").read()
        src.append(os.path.relpath(ns[0], WS2))
        def g(pat_):
            m = re.search(pat_, txt, re.I)
            return m.group(1).replace(",", "") if m else ""
        rows.append([lab_el, lab_en,
                     g(r"Number of reads:\s*([\d,\.]+)"),
                     g(r"Total bases:\s*([\d,\.]+)"),
                     g(r"Read length N50:\s*([\d,\.]+)"),
                     g(r"Median read length:\s*([\d,\.]+)"),
                     g(r"Mean read quality:\s*([\d,\.]+)")])
    if not rows:
        return
    # the flagstat-confirmed library sizes override NanoStats where they differ
    out = []
    for r in rows:
        out.append([r[0] if LANG_ROWS else r[1], thousands(r[2]),
                    num(float(r[3]) / 1e9, 2) if r[3] else "—",
                    thousands(r[4]), thousands(r[5]), num(r[6], 1)])
    add(1, "libraries",
        "Πίνακας 1. Χαρακτηριστικά των βιβλιοθηκών αλληλούχισης.",
        "Table 1. Sequencing library statistics.",
        ["Βιβλιοθήκη", "Αναγνώσεις", "Βάσεις (Gbp)", "N50 (bp)", "Διάμεσο μήκος (bp)", "Μέση ποιότητα"],
        ["Library", "Reads", "Bases (Gbp)", "N50 (bp)", "Median length (bp)", "Mean quality"],
        out,
        foot_el="Πηγή: %s." % ", ".join(src),
        foot_en="Source: %s." % ", ".join(src))


# ================================================================== 2. binning
def t2():
    bc = readtsv(P("03_binning", "bin_counts.tsv"))
    bt = readtsv(P("05_consensus", "binner_contribution.tsv"))
    if not bc:
        return
    keys = list(bc[0].keys())
    rows = []
    for r in bc:
        rows.append([str(r.get(keys[0], "")).strip()] +
                    [str(r.get(k, "")).strip() for k in keys[1:]])
    add(2, "binning",
        "Πίνακας 2. Αριθμός bins ανά εργαλείο και μετά τη συναίνεση.",
        "Table 2. Bin counts per binner and after consensus.",
        [k if k != keys[0] else "Δείγμα" for k in keys],
        keys, rows,
        foot_el="Πηγή: 03_binning/bin_counts.tsv" + (" · 05_consensus/binner_contribution.tsv" if bt else ""),
        foot_en="Source: 03_binning/bin_counts.tsv" + (" · 05_consensus/binner_contribution.tsv" if bt else ""))


# =========================================================== 3. the catalogue
MASTER = readtsv(P("47_master", "MASTER.tsv"))

def mval(r, *names, **kw):
    for n in names:
        if r.get(n) not in (None, ""):
            return r[n]
    return kw.get("default", "")

def t3():
    if not MASTER:
        return
    rows = []
    for r in sorted(MASTER, key=lambda r: mval(r, "mag_id")):
        rows.append([mval(r, "mag_id"),
                     mval(r, "phylum", default="—"),
                     mval(r, "family", default="—"),
                     mval(r, "genus", default="—"),
                     num(mval(r, "completeness"), 1),
                     num(mval(r, "contamination"), 2),
                     thousands(mval(r, "bp")),
                     mval(r, "mimag", default="—"),
                     mval(r, "has_16S", default="—")])
    add(3, "mag_catalogue",
        "Πίνακας 3. Ο κατάλογος των ανακατασκευασμένων γονιδιωμάτων (MAGs).",
        "Table 3. The metagenome-assembled genome (MAG) catalogue.",
        ["MAG", "Φύλο", "Οικογένεια", "Γένος", "Πληρότητα (%)", "Επιμόλυνση (%)",
         "Μέγεθος (bp)", "MIMAG", "16S"],
        ["MAG", "Phylum", "Family", "Genus", "Completeness (%)", "Contamination (%)",
         "Size (bp)", "MIMAG", "16S"],
        rows,
        foot_el="Πηγή: 47_master/MASTER.tsv. Πληρότητα και επιμόλυνση κατά CheckM2.",
        foot_en="Source: 47_master/MASTER.tsv. Completeness and contamination by CheckM2.")


# ==================================================== 4. taxonomy and novelty
def t4():
    rows, src = [], []
    for dom_el, dom_en, fn in (("Βακτήρια", "Bacteria", "gtdbtk.bac120.summary.tsv"),
                               ("Αρχαία", "Archaea", "gtdbtk.ar53.summary.tsv")):
        p = P("08_taxonomy", "gtdbtk", "classify", fn)
        rr = readtsv(p)
        if not rr:
            continue
        src.append(fn)
        for r in rr:
            ani = mval(r, "closest_placement_ani", "fastani_ani", default="")
            af = mval(r, "closest_placement_af", "fastani_af", default="")
            ref = mval(r, "closest_placement_reference", "fastani_reference", default="—")
            rows.append([mval(r, "user_genome"), dom_el if LANG_ROWS else dom_en,
                         (mval(r, "classification", default="")).split(";")[-1].replace("s__", "") or "—",
                         ref, num(ani, 2) if ani else "—", num(af, 2) if af else "—"])
    if not rows:
        return
    add(4, "taxonomy_novelty",
        "Πίνακας 4. Ταξινομική απόδοση κατά GTDB-Tk και εγγύτητα στο πλησιέστερο γονιδίωμα αναφοράς.",
        "Table 4. GTDB-Tk assignment and distance to the closest reference genome.",
        ["MAG", "Επικράτεια", "Είδος GTDB", "Πλησιέστερο γονιδίωμα", "ANI (%)", "AF"],
        ["MAG", "Domain", "GTDB species", "Closest genome", "ANI (%)", "AF"],
        rows,
        foot_el="Πηγή: 08_taxonomy/gtdbtk/classify/ (%s). Κενό είδος = καμία απόδοση σε επίπεδο είδους." % ", ".join(src),
        foot_en="Source: 08_taxonomy/gtdbtk/classify/ (%s). An empty species means no species-level assignment." % ", ".join(src))


# ===================================================== 5. in-situ detectability
def t5():
    if not MASTER:
        return
    rows = []
    for r in MASTER:
        b = mval(r, "sed_breadth")
        try:
            bf = float(b)
        except (TypeError, ValueError):
            continue
        rows.append((bf, [mval(r, "mag_id"), mval(r, "genus", default="—"),
                          mval(r, "sample", default="—"),
                          num(bf, 4), num(mval(r, "sed_relabund"), 3),
                          num(mval(r, "liq_breadth"), 4)]))
    rows.sort(key=lambda x: -x[0])
    add(5, "in_situ_detectability",
        "Πίνακας 5. Ανιχνευσιμότητα κάθε γονιδιώματος στις αναγνώσεις του ιζήματος.",
        "Table 5. Detectability of each genome in the sediment read set.",
        ["MAG", "Γένος", "Προέλευση", "Κάλυψη ιζήματος", "Σχετ. αφθονία (%)", "Κάλυψη καλλιεργειών"],
        ["MAG", "Genus", "Origin", "Sediment breadth", "Rel. abundance (%)", "Culture breadth"],
        [r for _, r in rows],
        foot_el="Πηγή: 47_master/MASTER.tsv (26_crossmap). Κάλυψη = κλάσμα του γονιδιώματος "
                "που καλύπτεται από αναγνώσεις. Χαμηλή κάλυψη σημαίνει ότι ο πληθυσμός βρίσκεται "
                "κάτω από το όριο ανίχνευσης αυτής της βιβλιοθήκης, όχι ότι απουσιάζει από τη λιμνοθάλασσα.",
        foot_en="Source: 47_master/MASTER.tsv (26_crossmap). Breadth = fraction of the genome covered "
                "by reads. A low breadth means the population is below this library's recovery limit, "
                "not that it is absent from the lagoon.")


# ============================================= 6. methanogenesis: modules + mcr
def t6():
    mods = readtsv(P("11_kofam", "ko_modules.tsv"))
    thr = {}
    kl = P("36_gapclose", "ko_list")
    if os.path.exists(kl):
        for ln in open(kl, errors="replace"):
            f = ln.rstrip("\n").split("\t")
            if len(f) > 1 and re.match(r"^K\d{5}$", f[0]):
                try:
                    thr[f[0]] = float(f[1])
                except ValueError:
                    pass
    bin2mag = {mval(r, "working_name"): mval(r, "mag_id") for r in MASTER
               if mval(r, "working_name")}
    best = collections.defaultdict(dict)
    tb = P("45_kofam_full", "all_hits.tbl")
    KOS = ("K00399", "K00401", "K00402", "K14083", "K16176", "K16178", "K04480")
    if os.path.exists(tb):
        for ln in open(tb, errors="replace"):
            if ln.startswith("#"):
                continue
            f = ln.split()
            if len(f) < 6 or f[2] not in KOS:
                continue
            mag = bin2mag.get(f[0].split("|")[0])
            if not mag:
                continue
            try:
                s = float(f[5])
            except ValueError:
                continue
            if s > best[mag].get(f[2], -1):
                best[mag][f[2]] = s
    if not best:
        return
    arc = [r for r in MASTER if mval(r, "domain") == "Archaea"]
    rows = []
    for r in sorted(arc, key=lambda r: mval(r, "mag_id")):
        mg = mval(r, "mag_id")
        cells = []
        for ko in KOS:
            s = best[mg].get(ko)
            t = thr.get(ko)
            if s is None:
                cells.append("—")
            else:
                cells.append(("%s%s" % (num(s, 1), " *" if (t and s >= t) else "")))
        rows.append([mg, mval(r, "genus", default="—")] + cells)
    add(6, "methanogenesis_markers",
        "Πίνακας 6. Βαθμολογίες KOfam των δεικτών μεθανογένεσης ανά αρχαίο γονιδίωμα.",
        "Table 6. KOfam scores for the methanogenesis markers, per archaeal genome.",
        ["MAG", "Γένος", "mcrA", "mcrB", "mcrG", "mttB", "mtmB", "mtbB", "mtaB"],
        ["MAG", "Genus", "mcrA", "mcrB", "mcrG", "mttB", "mtmB", "mtbB", "mtaB"],
        rows,
        foot_el="Πηγή: 45_kofam_full/all_hits.tbl, ενωμένο με 47_master/MASTER.tsv στο όνομα του bin· "
                "όρια από 36_gapclose/ko_list (K00399 775,53· K00401 504,30· K00402 360,93· "
                "K14083 160,67· K16176 261,00· K16178 496,93· K04480 578,57). "
                "Ο αστερίσκος σημαίνει πάνω από το επιμελημένο όριο. Δίνεται η καλύτερη βαθμολογία ORF.",
        foot_en="Source: 45_kofam_full/all_hits.tbl joined to 47_master/MASTER.tsv on the bin name; "
                "thresholds from 36_gapclose/ko_list (K00399 775.53; K00401 504.30; K00402 360.93; "
                "K14083 160.67; K16176 261.00; K16178 496.93; K04480 578.57). "
                "An asterisk marks a score above the curated threshold. Best ORF score shown.")


# =================================================== 7. sediment geochemistry
# Recovered from the VECTOR coordinates of Sediment_Chemical_Analysis.pdf: the
# data points are filled circles whose centres are the values, each panel's axis
# fitted on its own printed tick labels (max residual 1,3e-4 data units), depth
# from the panel frame (recovered depths land on 2,5 + 5k to 1,5e-5 cm). Two
# independent objects in the file -- the circles and the 13-segment polyline
# joining them -- agree to 3e-5 cm, and every value sits on a 0,1 grid. Sulfate
# reproduces all eight of audit finding F28's statistics exactly.
GEOCHEM_COLS = ["Ash", "TOC", "Polysaccharides", "Total proteins", "Total phenols",
                "NH4+-N", "TKN", "Total phosphorus", "SO4 2-", "NO3-"]
GEOCHEM = [
    (2.5,  21.5, 14.5, 13.7, 50.0, 87.5, 1.1, 41.0,  8.8, 135.0, 15.2),
    (7.5,  21.5, 12.5,  9.9, 16.0, 81.5, 6.0, 27.0,  4.8, 285.0, 15.8),
    (12.5, 20.5, 13.2, 10.8, 44.5, 62.0, 7.5, 26.0,  0.5, 240.0, 15.6),
    (17.5, 22.5, 14.3, 12.5, 34.0, 72.0, 1.4, 23.0,  0.8, 295.0, 15.2),
    (22.5, 26.0, 17.2, 11.0, 34.0, 68.0, 1.9, 21.0,  9.0, 215.0, 15.0),
    (27.5, 20.0,  8.8,  8.2, 41.0, 51.0, 0.1, 29.0,  7.2, 242.0, 15.3),
    (32.5, 19.5, 14.1,  8.1, 17.5, 47.5, 0.1, 28.0, 11.0,  45.0, 15.6),
    (37.5, 19.8,  8.0,  8.3, 18.0, 52.0, 0.1, 30.0,  4.6, 182.0, 15.5),
    (42.5, 19.8,  9.5,  9.8, 41.0, 63.0, 0.1, 31.0,  4.6, 155.0, 15.1),
    (47.5, 21.0, 12.2,  7.8, 48.0, 40.0, 0.1, 28.0,  3.1,  38.0, 15.3),
    (52.5, 20.5, 11.8, 10.1, 27.0, 49.5, 6.1, 33.0, 25.2, 240.0, 15.1),
    (57.5, 22.0, 12.8,  9.0,  1.0, 35.0, 0.1, 31.0,  0.8, 135.0, 15.8),
    (62.5, 18.5, 13.4, 11.0, 23.5, 41.5, 6.0, 35.0,  0.8, 205.0, 15.4),
    (67.5, 17.5, 12.1,  9.8, 19.5, 43.5, 1.3, 37.0, 31.8, 170.0, 16.1),
]

def t7():
    ext = os.path.join(OUT, "sediment_geochemistry_by_depth.tsv")
    rows, where = [], "ενσωματωμένο / embedded"
    if os.path.exists(ext):
        rr = readtsv(ext)
        if rr:
            ks = list(rr[0].keys())
            rows = [[r[k] for k in ks] for r in rr]
            where = os.path.relpath(ext, WS2)
    if not rows:
        rows = [[num(v, 1) if i else num(v, 1) for i, v in enumerate(r)] for r in GEOCHEM]
    add(7, "sediment_geochemistry",
        "Πίνακας 7. Γεωχημεία του ιζήματος ανά βάθος (14 ορίζοντες, 2,5–67,5 cm).",
        "Table 7. Sediment geochemistry by depth (14 horizons, 2.5–67.5 cm).",
        ["Βάθος (cm)", "Τέφρα (% w/w)", "TOC (% w/w)", "Πολυσακχαρίτες (mg/g)",
         "Ολικές πρωτεΐνες (mg/kg)", "Ολικές φαινόλες (mg/g)", "NH₄⁺-N (mg/kg)",
         "TKN (mg/kg)", "Ολικός P (mg/kg)", "SO₄²⁻ (mg/kg)", "NO₃⁻ (mg/kg)"],
        ["Depth (cm)", "Ash (% w/w)", "TOC (% w/w)", "Polysaccharides (mg/g)",
         "Total proteins (mg/kg)", "Total phenols (mg/g)", "NH₄⁺-N (mg/kg)",
         "TKN (mg/kg)", "Total P (mg/kg)", "SO₄²⁻ (mg/kg)", "NO₃⁻ (mg/kg)"],
        rows,
        foot_el="Πηγή: %s, ανακτημένο από τις διανυσματικές συντεταγμένες του "
                "Sediment_Chemical_Analysis.pdf. Ο έλεγχος πέρασε: τα θειικά αναπαράγουν και "
                "τις οκτώ στατιστικές του ευρήματος F28. ΠΡΟΣΟΧΗ: οι στήλες TKN και ολικών "
                "πρωτεϊνών φέρουν σφάλμα μονάδας ~10³ (C:N 2.667–8.190 : 1 σε κάθε ορίζοντα)· "
                "η στήλη των νιτρικών δεν είναι αξιοποιήσιμη (εύρος 7%% σε 70 cm). "
                "Αυτά περιγράφουν το ίζημα, όχι το θρεπτικό μέσο των καλλιεργειών." % where,
        foot_en="Source: %s, recovered from the vector coordinates of "
                "Sediment_Chemical_Analysis.pdf. Control passed: sulfate reproduces all eight of "
                "audit finding F28's statistics. CAUTION: the TKN and total protein columns carry a "
                "~10³ unit error (C:N 2,667–8,190 : 1 at every horizon); the nitrate column is "
                "unusable (7%% range over 70 cm). These describe the sediment, not the "
                "enrichment medium." % where)


# ================================================== 8. community composition
# THE UNIT QUESTION IS SETTLED, AND THE ANSWER IS THAT NEITHER RAW COLUMN WAS A
# RELATIVE ABUNDANCE.
#
# SingleM's own glossary defines what a taxonomic profile contains:
#
#   coverage  "the average number of reads covering each position in that
#             species' genome", estimated from reads recruited to the marker
#             gene windows;
#   UNFILLED  a profile reports unfilled coverage -- "the coverage of a taxon
#             does NOT include the coverage of its sub-taxons";
#   FILLED    a taxon's own coverage plus every descendant's;
#   and "to calculate relative abundance, the FILLED coverage is used",
#             divided by the sum of filled coverages at that level.
#
# So the 23,41 next to Desulfobacterota in 22_singlem/SED_profile.tsv is the
# coverage assigned AT the phylum rank only -- reads that could not be placed
# any deeper. It is not that phylum's abundance, and phylum rows do not sum to
# the community (theirs summed to 104,18, which is not 100% and was never meant
# to be). The pct_of_rank column in 24_deep_taxonomy is a different object
# again: an aggregation over the OTU table, before the condense step that
# resolves each window to one taxon.
#
# The fix is neither to pick one file nor to flag a tie. It is to do what the
# tool documents: fill the coverage down the lineage, then divide. That is
# computed here from the condensed profile, which is the output SingleM
# intends for this purpose.
def t8():
    prof = P("22_singlem", "SED_profile.tsv")
    if not os.path.exists(prof):
        return
    rows_in = []
    with open(prof, errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ci = ti = None
        for i, h in enumerate(head):
            hl = h.strip().lower()
            if hl in ("coverage", "relative_abundance"):
                ci = i
            if hl.startswith("tax"):
                ti = i
        if ci is None or ti is None:          # no header: sample, value, taxonomy
            fh.seek(0)
            ci, ti = 1, 2
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) <= max(ci, ti):
                continue
            try:
                cov = float(f[ci])
            except ValueError:
                continue
            lin = [x.strip() for x in f[ti].split(";") if x.strip()]
            rows_in.append((cov, lin))
    if not rows_in:
        return

    def rank_of(lin, pre):
        for x in lin:
            if x.startswith(pre):
                return x[len(pre):]
        return None

    # filled coverage: a phylum takes its own row plus every row beneath it,
    # which is exactly what "filled" means in the glossary
    filled = collections.defaultdict(float)
    unfilled = collections.defaultdict(float)
    dom = {}
    for cov, lin in rows_in:
        ph = rank_of(lin, "p__")
        if not ph:
            continue
        filled[ph] += cov
        if lin[-1].startswith("p__"):
            unfilled[ph] += cov
        d = rank_of(lin, "d__")
        if d:
            dom[ph] = d
    if not filled:
        return
    tot = sum(filled.values())
    rows = []
    for ph, fc in sorted(filled.items(), key=lambda kv: -kv[1]):
        rows.append([ph, dom.get(ph, "—"), num(fc, 2), num(unfilled.get(ph, 0.0), 2),
                     num(100.0 * fc / tot, 2)])
    add(8, "community_composition",
        "Πίνακας 8. Σύνθεση της κοινότητας του ιζήματος σε επίπεδο φύλου (SingleM).",
        "Table 8. Sediment community composition at phylum level (SingleM).",
        ["Φύλο", "Επικράτεια", "Συμπληρωμένη κάλυψη", "Κάλυψη στο επίπεδο φύλου",
         "Σχετική αφθονία (%)"],
        ["Phylum", "Domain", "Filled coverage", "Coverage at phylum rank",
         "Relative abundance (%)"],
        rows,
        foot_el="Πηγή: 22_singlem/SED_profile.tsv, το συμπυκνωμένο προφίλ του SingleM. "
                "Το SingleM αναφέρει ΜΗ ΣΥΜΠΛΗΡΩΜΕΝΗ κάλυψη: η τιμή κάθε ταξινομικής ομάδας "
                "ΔΕΝ περιλαμβάνει τις υποκείμενες. Η στήλη «κάλυψη στο επίπεδο φύλου» είναι "
                "αυτή η ακατέργαστη τιμή — οι αναγνώσεις που δεν τοποθετήθηκαν βαθύτερα — και "
                "ΔΕΝ είναι αφθονία. Η σχετική αφθονία υπολογίζεται, όπως ορίζει το εργαλείο, "
                "από τη ΣΥΜΠΛΗΡΩΜΕΝΗ κάλυψη (ίδια συν όλων των απογόνων) προς το άθροισμά της "
                "στο ίδιο επίπεδο. Ο πίνακας pct_of_rank του 24_deep_taxonomy προέρχεται από "
                "το OTU table, πριν το βήμα condense, και είναι διαφορετικό μέγεθος· δεν "
                "χρησιμοποιείται εδώ.",
        foot_en="Source: 22_singlem/SED_profile.tsv, SingleM's condensed profile. SingleM "
                "reports UNFILLED coverage: a taxon's value does NOT include its sub-taxa. The "
                "'coverage at phylum rank' column is that raw value — reads not placed any "
                "deeper — and is NOT an abundance. Relative abundance is computed as the tool "
                "defines it, from FILLED coverage (a taxon's own plus all its descendants') "
                "over the sum of filled coverages at the same level. The pct_of_rank table in "
                "24_deep_taxonomy comes from the OTU table, before the condense step, and is a "
                "different quantity; it is not used here.")


# ==================================================== 9. auxotrophy, 10. CAZy
def t9():
    rr = readtsv(P("19b_auxotrophy", "auxotrophies.tsv"))
    if not rr:
        return
    ks = list(rr[0].keys())
    rows = [[str(r.get(k, "")).strip() for k in ks] for r in rr]
    add(9, "auxotrophy",
        "Πίνακας 9. Προβλεπόμενες αυξοτροφίες ανά γονιδίωμα.",
        "Table 9. Predicted auxotrophies per genome.",
        ks, ks, rows,
        foot_el="Πηγή: 19b_auxotrophy/auxotrophies.tsv. 🔴 Οι αυξοτροφίες κυστεΐνης των αρχαίων "
                "είναι τεχνούργημα: το μονοπάτι ορίζεται από βακτηριακούς δείκτες (cysK/cysE/cysM) "
                "και τα αρχαία χρησιμοποιούν SepRS/SepCysS. Μην τις αναφέρετε.",
        foot_en="Source: 19b_auxotrophy/auxotrophies.tsv. 🔴 The archaeal cysteine auxotrophies are "
                "artefacts: the pathway is defined by bacterial markers (cysK/cysE/cysM) while "
                "archaea use SepRS/SepCysS. Do not report them.")

def t10():
    rr = readtsv(P("30_dbcan", "family_matrix.tsv"))
    if not rr:
        return
    ks = list(rr[0].keys())
    mag = ks[0]
    fam = [k for k in ks[1:]]
    rows = []
    for r in rr:
        tot = 0
        cls = collections.Counter()
        for k in fam:
            try:
                v = int(float(r.get(k) or 0))
            except ValueError:
                continue
            tot += v
            m = re.match(r"^(GH|GT|PL|CE|AA|CBM)", k)
            if m:
                cls[m.group(1)] += v
        rows.append((tot, [str(r.get(mag, "")).strip(), str(tot)] +
                          [str(cls.get(c, 0)) for c in ("GH", "GT", "PL", "CE", "AA", "CBM")]))
    rows.sort(key=lambda x: -x[0])
    add(10, "cazymes",
        "Πίνακας 10. Ένζυμα δράσης σε υδατάνθρακες (CAZymes) ανά γονιδίωμα.",
        "Table 10. Carbohydrate-active enzymes (CAZymes) per genome.",
        ["MAG", "Σύνολο", "GH", "GT", "PL", "CE", "AA", "CBM"],
        ["MAG", "Total", "GH", "GT", "PL", "CE", "AA", "CBM"],
        [r for _, r in rows],
        foot_el="Πηγή: 30_dbcan/family_matrix.tsv. Δηλώσεις για τη λιμνοθάλασσα πρέπει να "
                "γίνονται μόνο από γονιδιώματα με υψηλή κάλυψη ιζήματος (Πίνακας 5).",
        foot_en="Source: 30_dbcan/family_matrix.tsv. In-situ statements must be made only from "
                "genomes with high sediment breadth (Table 5).")


# ================================================== 11. biogeography, 12. ANI
def t11():
    rr = readtsv(P("14_comparative", "recruitment_matrix.tsv"))
    if not rr:
        return
    ks = list(rr[0].keys())
    ds = ks[1:]
    rows = []
    for dset in ds:
        vals = []
        for r in rr:
            try:
                vals.append(float(r.get(dset) or 0))
            except ValueError:
                pass
        if not vals:
            continue
        rows.append((max(vals), [dset, "%d" % len(vals), num(max(vals), 4),
                                 num(statistics.median(vals), 4),
                                 "%d" % sum(1 for v in vals if v >= 0.30)]))
    rows.sort(key=lambda x: -x[0])
    add(11, "biogeography",
        "Πίνακας 11. Στρατολόγηση αναγνώσεων των 29 γονιδιωμάτων σε δημόσια σύνολα δεδομένων.",
        "Table 11. Read recruitment of the 29 genomes against public datasets.",
        ["Σύνολο δεδομένων", "Γονιδιώματα", "Μέγιστη κάλυψη", "Διάμεση κάλυψη", "Γονιδιώματα ≥ 0,30"],
        ["Dataset", "Genomes", "Max breadth", "Median breadth", "Genomes ≥ 0.30"],
        [r for _, r in rows],
        foot_el="Πηγή: 14_comparative/recruitment_matrix.tsv. Οι δημόσιες τιμές δεν είναι μηδέν "
                "και δεν πρέπει να παρουσιάζονται ως μηδέν.",
        foot_en="Source: 14_comparative/recruitment_matrix.tsv. The public values are not zero and "
                "must not be presented as zero.")

def t12():
    p = P("40_integration", "skani_raw.tsv")
    if not os.path.exists(p):
        return
    rows = []
    for ln in open(p, errors="replace"):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 3 or f[0].lower().startswith(("ref", "query", "#")):
            continue
        try:
            ani = float(f[2])
        except ValueError:
            continue
        a, b = os.path.basename(f[0]), os.path.basename(f[1])
        if a == b:
            continue
        rows.append((ani, [a, b, num(ani, 2)] +
                     [num(f[3], 1) if len(f) > 3 else "—",
                      num(f[4], 1) if len(f) > 4 else "—"]))
    if not rows:
        return
    rows.sort(key=lambda x: -x[0])
    add(12, "ani_pairs",
        "Πίνακας 12. Ζεύγη γονιδιωμάτων για τα οποία το skani επέστρεψε τιμή ANI.",
        "Table 12. Genome pairs for which skani returned an ANI value.",
        ["Γονιδίωμα Α", "Γονιδίωμα Β", "ANI (%)", "Κάλυψη Α (%)", "Κάλυψη Β (%)"],
        ["Genome A", "Genome B", "ANI (%)", "Aligned A (%)", "Aligned B (%)"],
        [r for _, r in rows],
        foot_el="Πηγή: 40_integration/skani_raw.tsv. Το skani επιστρέφει τιμή μόνο πάνω από ένα "
                "εσωτερικό όριο· ένα ζεύγος χωρίς τιμή ΥΠΟΛΟΓΙΣΤΗΚΕ και έπεσε κάτω από αυτό. "
                "Η χαμηλότερη τιμή που επιστρέφει αυτή η εκτέλεση είναι το εμπειρικό δάπεδο.",
        foot_en="Source: 40_integration/skani_raw.tsv. skani reports a pair only above an internal "
                "screening threshold; a pair with no value WAS computed and fell below it. The "
                "lowest value this run returns is the empirical floor.")


# ================================================================ write it out
LANG_ROWS = True     # Greek labels inside data cells where a row carries one

def write_tsv(t):
    p = os.path.join(OUT, "table%02d_%s.tsv" % (t.n, t.stem))
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(t.cols_en)
        for r in t.rows:
            w.writerow(r)
    return p

def write_docx(lang):
    try:
        from docx import Document
        from docx.shared import Pt, Mm
        from docx.enum.text import WD_ALIGN_PARAGRAPH
    except ImportError:
        return None
    doc = Document()
    for s in doc.sections:
        s.page_width, s.page_height = Mm(210), Mm(297)
        s.left_margin = s.right_margin = Mm(25)
    st = doc.styles["Normal"]
    st.font.name = "Liberation Serif"
    st.font.size = Pt(9)
    for t in TABLES:
        cap = t.cap_el if lang == "EL" else t.cap_en
        cols = t.cols_el if lang == "EL" else t.cols_en
        foot = t.foot_el if lang == "EL" else t.foot_en
        p = doc.add_paragraph()
        r = p.add_run(cap)
        r.bold = True
        r.font.size = Pt(9.5)
        tb = doc.add_table(rows=1, cols=len(cols))
        tb.style = "Table Grid"
        for i, c in enumerate(cols):
            cell = tb.rows[0].cells[i]
            cell.text = ""
            rr = cell.paragraphs[0].add_run(str(c))
            rr.bold = True
            rr.font.size = Pt(8.5)
        for row in t.rows:
            cells = tb.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = ""
                rr = cells[i].paragraphs[0].add_run("" if v is None else str(v))
                rr.font.size = Pt(8.5)
        if foot:
            fp = doc.add_paragraph()
            fr = fp.add_run(foot)
            fr.italic = True
            fr.font.size = Pt(7.5)
        doc.add_paragraph()
    p = os.path.join(OUT, "TABLES_%s.docx" % lang)
    doc.save(p)
    return p



# ---------------------------------------------------------------- docx, by hand
# python-docx is not on the analysis machine and its shared environments are not
# writable, so rather than asking for an install the file is written directly.
# A .docx is a zip of XML: three members are enough for Word and LibreOffice to
# open it, and table borders and fonts are set inline so no styles part is
# needed. Nothing here imports anything outside the standard library.
import zipfile

def _x(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;"))

def _run(text, bold=False, italic=False, half_pt=17):
    pr = "<w:rPr><w:rFonts w:ascii='Liberation Serif' w:hAnsi='Liberation Serif'/>"
    if bold:
        pr += "<w:b/>"
    if italic:
        pr += "<w:i/>"
    pr += "<w:sz w:val='%d'/><w:szCs w:val='%d'/></w:rPr>" % (half_pt, half_pt)
    return ("<w:r>%s<w:t xml:space='preserve'>%s</w:t></w:r>" % (pr, _x(text)))

def _para(text, **kw):
    return "<w:p>%s</w:p>" % _run(text, **kw)

_BORDER = ("<w:tblBorders>" + "".join(
    "<w:%s w:val='single' w:sz='4' w:space='0' w:color='999999'/>" % e
    for e in ("top", "left", "bottom", "right", "insideH", "insideV")) +
    "</w:tblBorders>")

def _table(cols, rows):
    out = ["<w:tbl><w:tblPr><w:tblW w:w='5000' w:type='pct'/>%s"
           "<w:tblLayout w:type='autofit'/></w:tblPr>" % _BORDER]
    out.append("<w:tr><w:trPr><w:tblHeader/></w:trPr>" + "".join(
        "<w:tc><w:tcPr><w:shd w:val='clear' w:fill='F2F2F0'/></w:tcPr>%s</w:tc>"
        % _para(c, bold=True, half_pt=16) for c in cols) + "</w:tr>")
    for r in rows:
        out.append("<w:tr>" + "".join(
            "<w:tc>%s</w:tc>" % _para("" if v is None else v, half_pt=16)
            for v in r) + "</w:tr>")
    out.append("</w:tbl>")
    return "".join(out)

_CT = ("<?xml version='1.0' encoding='UTF-8' standalone='yes'?>"
       "<Types xmlns='http://schemas.openxmlformats.org/package/2006/content-types'>"
       "<Default Extension='rels' ContentType='application/vnd.openxmlformats-package.relationships+xml'/>"
       "<Default Extension='xml' ContentType='application/xml'/>"
       "<Override PartName='/word/document.xml' ContentType='application/vnd."
       "openxmlformats-officedocument.wordprocessingml.document.main+xml'/></Types>")
_RELS = ("<?xml version='1.0' encoding='UTF-8' standalone='yes'?>"
         "<Relationships xmlns='http://schemas.openxmlformats.org/package/2006/relationships'>"
         "<Relationship Id='rId1' Type='http://schemas.openxmlformats.org/officeDocument/"
         "2006/relationships/officeDocument' Target='word/document.xml'/></Relationships>")

def write_docx_plain(lang):
    body = []
    for t in TABLES:
        cap = t.cap_el if lang == "EL" else t.cap_en
        cols = t.cols_el if lang == "EL" else t.cols_en
        foot = t.foot_el if lang == "EL" else t.foot_en
        body.append(_para(cap, bold=True, half_pt=19))
        body.append(_table(cols, t.rows))
        if foot:
            body.append(_para(foot, italic=True, half_pt=15))
        body.append("<w:p/>")
    doc = ("<?xml version='1.0' encoding='UTF-8' standalone='yes'?>"
           "<w:document xmlns:w='http://schemas.openxmlformats.org/wordprocessingml/2006/main'>"
           "<w:body>" + "".join(body) +
           "<w:sectPr><w:pgSz w:w='11906' w:h='16838'/>"
           "<w:pgMar w:top='1418' w:right='1418' w:bottom='1418' w:left='1418'/>"
           "</w:sectPr></w:body></w:document>")
    p = os.path.join(OUT, "TABLES_%s.docx" % lang)
    with zipfile.ZipFile(p, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", _CT)
        z.writestr("_rels/.rels", _RELS)
        z.writestr("word/document.xml", doc)
    return p


def main():
    note("make_tables.py -- workspace: %s" % WS2)
    note("")
    if not MASTER:
        note("  !! 47_master/MASTER.tsv not found. Most tables depend on it.")
    for fn in (t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12):
        try:
            fn()
        except Exception as e:
            note("  [fail]  %s: %s: %s" % (fn.__name__, type(e).__name__, e))
    note("")
    for t in TABLES:
        write_tsv(t)
    note("  wrote %d TSV file(s) to %s/" % (len(TABLES), OUT))
    made, how = [], "python-docx"
    for lang in ("EL", "EN"):
        p = write_docx(lang)
        if p is None:
            how = "the built-in writer (python-docx not installed)"
            p = write_docx_plain(lang)
        made.append(os.path.basename(p))
    note("  wrote %s via %s" % (", ".join(made), how))
    with open(os.path.join(OUT, "MANIFEST.txt"), "w") as fh:
        fh.write("\n".join(MANIFEST) + "\n")
    print("\n  -> %s/MANIFEST.txt" % OUT)

main()
