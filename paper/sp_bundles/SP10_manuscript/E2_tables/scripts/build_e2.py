#!/usr/bin/env python3
"""SP10 E2: export the supplementary tables. Re-reads existing outputs only; no analysis is re-run.

Every input is logged (path + md5) to ../../inputs_md5.tsv; items that cannot be filled are
written as NA and listed in E2_notes.tsv. Run from anywhere:  python3 build_e2.py
"""
import csv, gzip, hashlib, json, os, re, sys
from collections import Counter, OrderedDict, defaultdict

WS2 = "${WS2}"
WS1 = "${WS1}"
M = os.path.join(WS2, "99_thesis_audit/SP10_manuscript")
OUT = os.path.join(M, "E2_tables")
INPUTS = OrderedDict()
NOTES = []

def P(rel):
    p = rel if rel.startswith("/") else os.path.join(WS2, rel)
    if p not in INPUTS:
        h = hashlib.md5()
        with open(p, "rb") as f:
            for b in iter(lambda: f.read(1 << 20), b""):
                h.update(b)
        INPUTS[p] = h.hexdigest()
    return p

def rows(rel, header=True):
    with open(P(rel), newline="") as f:
        r = csv.reader(f, delimiter="\t")
        data = [x for x in r if x]
    if not header:
        return data
    h = data[0]
    return [dict(zip(h, x)) for x in data[1:]]

def write(name, header, data):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in data:
            w.writerow(["NA" if (v is None or v == "") else v for v in r])
    print("  %-30s %4d rows" % (name, len(data)))

def note(table, item, why):
    NOTES.append((table, item, why))

# ------------------------------------------------------------------ ID mapping
CW = rows("07_mags/MAG_CROSSWALK.tsv")
W2M = {r["working_name"]: r["mag_id"] for r in CW}
W2M.update({r["qc_bin_name"]: r["mag_id"] for r in CW})
W2M["LIQ_das_unknown_000b"] = "ETO_LIQ_MAG_025"   # CheckM2 alias (SP10 NCBI)
MAGS = sorted(r["mag_id"] for r in CW)

def mid(x):
    """working name / short name / bin~contig -> final MAG ID"""
    x = x.split("~")[0]
    if x in W2M: return W2M[x]
    if re.fullmatch(r"(LIQ|SED)_MAG_\d{3}", x): return "ETO_" + x
    return x

# ================================================================== T1 versions
def t1():
    out, seen = [], set()
    def add(tool, ver, db, src, note_=""):
        k = (tool.strip().lower(), ver.strip())
        if k in seen: return
        seen.add(k); out.append([tool.strip(), ver.strip(), db, src, note_.strip()])
    for r in rows("versions.tsv", header=False)[1:]:
        r += [""] * 3
        add(r[0], r[1], "", "workspace_2/versions.tsv", "env: " + r[2] if r[2] else "")
    for sp in ("SP01_markers", "SP02_recruitment", "SP03_16S_links", "SP04_archaeal_orders",
               "SP05_dsr_mcr", "SP09_final"):
        rel = "99_thesis_audit/%s/versions.tsv" % sp
        for r in rows(rel, header=False)[1:]:
            r += [""] * 3
            dbish = any(w in r[0].lower() for w in ("db", "database", "gtdb", "release", "kofam", "silva", "reference"))
            add(r[0], r[1] if not dbish else "", r[1] if dbish else "", sp + "/versions.tsv", r[2])
    write("T_versions.tsv", ["tool_or_item", "version", "database_or_release", "source_file", "note"], out)
    note("T_versions", "duplicates", "rows merged on identical (tool, version); different versions of the same tool in different steps are kept as separate rows")

# ================================================================== T2 read QC
def t2():
    out = []
    mq = {r[0]: r for r in rows("99_thesis_audit/SP10_ncbi_submission/metadata/read_mean_quality.tsv", header=False)}
    for lib, fn, key in (("SED", "01_qc/SED_readstats_prepost.tsv", "ETO_SED_shotgun.fastq.gz"),
                         ("LIQ", "01_qc/LIQ_readstats_prepost.tsv", "ETO_LIQ_shotgun.fastq.gz")):
        rr = rows(fn)
        for i, r in enumerate(rr):
            stage = "before >=1 kb filter (deduplicated)" if i == 0 else "after >=1 kb filter (Flye input)"
            meanq = mq[key][2] if i == 1 and key in mq else ""
            out.append([lib, stage, os.path.relpath(r["file"], WS2), r["num_seqs"], r["sum_len"],
                        r.get("N50", ""), r.get("Q2", ""), meanq, "", r.get("GC(%)", ""), "", ""])
    for lib, d in (("SED", "sediment"), ("LIQ", "liquid_cultures")):
        f = os.path.join(WS1, "provenance_harvest/02_qc", d, "NanoStats.txt")
        t = open(P(f)).read()
        g = lambda k: re.search(k + r":\s+([\d,\.]+)", t).group(1).replace(",", "")
        q15 = re.search(r">Q15:\s+\d+ \(([\d\.]+)%\)", t).group(1)
        out.append([lib, "raw, all runs pooled (NanoPlot; SED includes the duplicated ligation file)",
                    os.path.relpath(f, os.path.dirname(WS2)), g("Number of reads"), g("Total bases"),
                    g("Read length N50"), g("Median read length"), g("Mean read quality"), q15, "", "", ""])
    note("T_read_qc", "mean Q / % > Q15 before filter", "not recorded for the deduplicated pre-filter file; only NanoPlot on the raw pooled reads (row 'raw')")
    note("T_read_qc", "% > Q15 after filter", "not recorded for *_min1000 files")
    note("T_read_qc", "GC of raw", "NanoPlot does not report GC")
    # per-run (sediment runs + LIQ runs): seqkit stats -a computed for SP10 on the original per-run FASTQs
    st = {os.path.basename(r["file"]): r for r in rows(os.path.join(OUT, "scripts/per_run_seqkit_stats.tsv"))}
    kit = {
        "ligation_new.fastq": ("Ligation (protocol_group_id GALIATSOS_LIGATION; LSK114 not written)", "no", "in assembly set"),
        "Ligation_promithion.fastq": ("Ligation (byte-identical copy of ligation_new.fastq)", "no", "duplicate; removed by dedup"),
        "Rapid_PCR_promithion.fastq": ("Rapid PCR Barcoding, RPB (protocol_group_id GALIATSOS_RAPID_PCR_BARCODING)", "yes", "in assembly set"),
        "TSIAMIS_EXP-PBC096.fastq": ("PCR Barcoding expansion EXP-PBC096 (file name)", "yes", "in assembly set"),
        "TSIAMIS_SQK-NBD114-24.fastq": ("Native Barcoding SQK-NBD114.24 (file name)", "no", "in assembly set"),
        "metagenome_LSK_test.fastq": ("Ligation LSK (file name)", "no", "in assembly set"),
        "metagenome_RBK_test.fastq": ("Rapid Barcoding RBK (file name); 16S amplicon library", "amplicon", "EXCLUDED (step 13f)"),
        "metagenome_2.fastq": ("not recorded", "unknown", "in assembly set"),
        "metagenome_3.fastq": ("not recorded", "unknown", "in assembly set"),
        "metagenome_new.fastq": ("not recorded", "unknown", "in assembly set"),
    }
    cnt = {os.path.basename(r[0]): r[1] for r in rows(os.path.join(WS1, "provenance_harvest/02_qc/per_run_read_counts.tsv"), header=False)}
    for fn in sorted(st):
        s = st[fn]
        lib = "SED" if fn in kit else "LIQ"
        k, pcr, use = kit.get(fn, ("not recorded", "unknown", "in assembly set"))
        out.append([lib, "per run: " + fn, "workspace/liquid_cultures/profile/fastqs/" + ("Sediment_Metagenome/" if lib == "SED" else "Liquid_Cultures_Metagenome/") + fn,
                    s["num_seqs"], s["sum_len"], s.get("N50", ""), s.get("Q2", ""), "", "", s.get("GC(%)", ""), k, pcr + "; " + use])
        if cnt.get(fn) and cnt[fn] != s["num_seqs"]:
            note("T_read_qc", fn, "seqkit reads %s != per_run_read_counts.tsv %s" % (s["num_seqs"], cnt[fn]))
    write("T_read_qc.tsv", ["library", "stage", "file", "reads", "bases_bp", "N50_bp", "median_length_bp",
                            "mean_Q", "pct_reads_gt_Q15", "GC_pct", "library_kit", "PCR_based_and_use"], out)

# ================================================================== T3 assembly
def fai_stats(p):
    L = sorted((int(r[1]) for r in rows(p, header=False)), reverse=True)
    tot, acc, n50 = sum(L), 0, 0
    for x in L:
        acc += x
        if acc >= tot / 2: n50 = x; break
    return len(L), tot, n50, L[0]

def t3():
    out = []
    for lib, d in (("SED", "sediment"), ("LIQ", "liquid_cultures")):
        info = rows("03_assembly/%s/assembly_info.txt" % d, header=False)
        circ = sum(1 for r in info[1:] if r[3] == "Y")
        n, tot, n50, mx = fai_stats("03_assembly/%s/assembly.fasta.fai" % d)
        out.append([lib, "Flye 2.9.3 (--nano-hq --meta)", "03_assembly/%s/assembly.fasta" % d, n, tot, n50, mx, circ])
        n, tot, n50, mx = fai_stats("04_polish/%s_asm_final.fasta.%s" % (lib, "fai" if lib == "LIQ" else "seqkit.fai"))
        out.append([lib, "Medaka 2.1.1 polished (binning input)", "04_polish/%s_asm_final.fasta" % lib, n, tot, n50, mx, circ])
    write("T_assembly.tsv", ["library", "assembly", "file", "contigs", "total_length_bp", "N50_bp", "largest_contig_bp",
                             "circular_contigs_Flye"], out)
    note("T_assembly", "QUAST", "QUAST was never run on these assemblies (only 00_smoketest); values read from Flye assembly_info.txt and samtools/seqkit .fai; circular count is Flye's circ.=Y (polishing does not change topology)")

# ================================================================== T4 binning
def t4():
    bc = rows("03_binning/bin_counts.tsv")
    das = rows("05_consensus/DAS_CROSSWALK.tsv")
    dasn = Counter((r["sample"], r["source_binner"]) for r in das)
    fin = Counter((r["sample"], r["binner"]) for r in CW)
    out = []
    for r in bc:
        k = (r["sample"], r["binner"])
        out.append([r["sample"], r["binner"], r["n_bins"], r["total_Mbp"], r["largest_Mb"], r["median_Mb"],
                    dasn.get(k, 0), fin.get(k, 0)])
    out.append(["all", "all", sum(int(r["n_bins"]) for r in bc), "", "", "", len(das), len(CW)])
    write("T_binning.tsv", ["sample", "binner", "bins_before_DAS_Tool", "total_binned_Mbp", "largest_bin_Mb", "median_bin_Mb",
                            "bins_kept_by_DAS_Tool", "final_MAGs_contributed"], out)
    note("T_binning", "31 vs 29", "DAS Tool kept 31 bins (05_consensus/DAS_CROSSWALK.tsv); 29 passed to the MAG catalogue (07_mags/MAG_CROSSWALK.tsv)")

# ================================================================== T5 catalogue
def fasta_stats(p):
    L, gc, at = [], 0, 0
    with open(P(p)) as f:
        cur = 0
        for ln in f:
            if ln.startswith(">"):
                if cur: L.append(cur)
                cur = 0
            else:
                s = ln.strip().upper(); cur += len(s)
                gc += s.count("G") + s.count("C"); at += s.count("A") + s.count("T")
        if cur: L.append(cur)
    L.sort(reverse=True); tot, acc = sum(L), 0
    for x in L:
        acc += x
        if acc >= tot / 2: return round(100 * gc / (gc + at), 2), x
    return "", ""

def t5():
    gt = {}
    for f in ("08_taxonomy/gtdbtk/classify/gtdbtk.bac120.summary.tsv", "08_taxonomy/gtdbtk/classify/gtdbtk.ar53.summary.tsv"):
        for r in rows(f): gt[mid(r["user_genome"])] = r
    ab = {r["mag_id"]: r for r in rows("26_crossmap/cross_abundance.tsv")}
    an = {r["mag_id"]: r for r in rows("09_annotation/annotation_summary.tsv")}
    sp2 = defaultdict(dict)
    for r in rows("99_thesis_audit/SP02_recruitment/isolate_calls.tsv"):
        if r["index"] == "A": sp2[r["genome"]][r["sample"]] = r
    org = {r["isolate"]: r["organism"] for r in rows("NCBI_submission/metadata/biosample_MAGs.tsv")}
    # circular: Flye circ flag via polish rename map
    circ = {}
    for lib, d in (("SED", "sediment"), ("LIQ", "liquid_cultures")):
        ci = {r[0]: r[3] for r in rows("03_assembly/%s/assembly_info.txt" % d, header=False)[1:]}
        for r in rows("04_polish/%s_contig_rename.tsv" % lib, header=False):
            circ[r[-1]] = ci.get(r[0], ci.get(r[1], "")) if len(r) > 1 else ""
    cwcols = list(CW[0].keys())
    hdr = cwcols + ["gtdb_r232_classification", "gtdb_red_value", "gtdb_closest_reference", "gtdb_closest_ani_pct",
                    "gtdb_closest_af", "GC_pct", "N50_bp", "CDS_count", "n_circular_contigs",
                    "sed_breadth", "sed_depth_x", "sed_relabund_pct", "liq_breadth", "liq_depth_x", "liq_relabund_pct",
                    "SP02_SED_identity_mode_pct", "SP02_SED_call_calibrated", "SP02_LIQ_identity_mode_pct",
                    "SP02_LIQ_call_calibrated", "NCBI_organism_name", "NCBI_isolate_name"]
    out = []
    for c in CW:
        m = c["mag_id"]; g = gt.get(m, {}); a = ab[m]
        gc, n50 = fasta_stats("07_mags/bins/%s.fa" % c["working_name"])
        contigs = [ln[1:].split()[0] for ln in open(os.path.join(WS2, "07_mags/bins/%s.fa" % c["working_name"])) if ln.startswith(">")]
        nc = sum(1 for x in contigs if circ.get(x) == "Y") if circ else ""
        s, l = sp2[m].get("SED", {}), sp2[m].get("LIQ", {})
        out.append([c[k] for k in cwcols] + [g.get("classification"), g.get("red_value"), g.get("closest_genome_reference"),
                    g.get("closest_genome_ani"), g.get("closest_genome_af"), gc, n50, an.get(m, {}).get("n_cds"), nc,
                    a["sed_breadth"], a["sed_depth"], a["sed_relabund"], a["liq_breadth"], a["liq_depth"], a["liq_relabund"],
                    s.get("hist_mode"), s.get("call_calibrated"), l.get("hist_mode"), l.get("call_calibrated"), org.get(m), m])
    write("T_mag_catalogue_full.tsv", hdr, out)
    note("T_mag_catalogue_full", "cross_abundance", "re-read from the raw 26_crossmap/cross_abundance.tsv (values keep their leading 0.)")
    note("T_mag_catalogue_full", "GC / N50", "computed from 07_mags/bins/<working_name>.fa (sequence re-read, not re-assembled)")
    note("T_mag_catalogue_full", "RED", "GTDB-Tk reports red_value only for genomes placed by RED (novel at some rank); N/A otherwise")

# ================================================================== T6/T7 public runs
def t6_t7():
    proj = {r[0]: r for r in rows("14_comparative/meta/projects.tsv", header=False)}
    rm = rows("14_comparative/recruitment_matrix.tsv")
    cols = [c for c in rm[0] if c.endswith("_breadth") and not c.startswith("ETOLIKO")]
    ri = {}
    with open(P("14_comparative/meta/runinfo_all.csv")) as f:
        for r in csv.DictReader(f): ri[r["Run"]] = r
    out = []
    for c in cols:
        pj, run = c.split("_")[0], c.split("_")[1]
        r = ri.get(run, {})
        out.append([pj, run, proj.get(pj, ["", "", ""])[1], proj.get(pj, ["", "", ""])[2], r.get("BioSample"),
                    r.get("Platform"), r.get("Model"), r.get("LibraryStrategy"), r.get("spots"), 2000000,
                    "searched"])
    for pj in proj:
        if not any(c.startswith(pj + "_") for c in cols):
            why = "RNA metatranscriptome, deliberately excluded" if proj[pj][2] == "RNA" else "no runs retrievable (meta/%s_runs.txt empty)" % pj
            out.append([pj, "", proj[pj][1], proj[pj][2], "", "", "", "", "", "", "NOT SEARCHED: " + why])
    write("T_public_runs.tsv", ["BioProject", "SRA_run", "environment", "nucleic_acid", "BioSample", "platform", "model",
                                "library_strategy", "spots_in_run", "spots_used_fastq_dump_X", "status"], out)
    note("T_public_runs", "country", "not in runinfo/projects.tsv; country is only implicit in the environment text")
    note("T_public_runs", "spots used", "2,000,000 per run (SPOTS default in scripts/14_comparative.sh; RUNLOG 2026-08-19 'spots=2000000 runs=2'); fastq-dump -X takes the FIRST N spots")
    # T7
    hdr = ["mag_id", "domain", "genus"] + [c.replace("_breadth", "") + "_breadth" for c in cols]
    write("T_recruitment_public.tsv", hdr, [[mid(r["mag_id"]), r["domain"], r["genus"]] + [r[c] for c in cols] for r in rm])

# ================================================================== T8 KOfam panels
def t8():
    fp = rows("45_kofam_full/tables/function_panels.tsv")
    gs = [g for g in fp[0] if g not in ("panel", "n_KO_in_panel")]
    # KO lists per panel, if recorded in genome_ko_summary/json
    out = []
    for r in fp:
        out.append([r["panel"], r["n_KO_in_panel"]] + [r[g] for g in gs])
    write("T_kofam_panels.tsv", ["panel", "n_KO_in_panel"] + [mid(g) + "_KOs_present" for g in gs], out)
    note("T_kofam_panels", "panels", "45_kofam_full/tables/function_panels.tsv holds 12 panels (values = number of panel KOs above the curated threshold). Panels requested but not in that file (methylotrophy split, B12-dependent, ECF/BtuCDF, oxidative stress, osmolytes, metal resistance, sporulation) are in T_b12_ecf / T_env_stress where available; others NA")
    note("T_kofam_panels", "KO list column", "the panel -> KO definition is not stored in any table under 45_kofam_full (it lived in the step-45 script); left NA rather than re-derived")

# ================================================================== T9 CAZy
def t9():
    fam = rows("30_dbcan/cazyme_families.tsv")
    cnt = defaultdict(Counter)
    for r in fam: cnt[r["mag_id"]][r["primary_family"]] += 1
    fams = sorted({f for c in cnt.values() for f in c}, key=lambda x: (re.sub(r"\d.*", "", x), int(re.sub(r"\D", "", x) or 0), x))
    fm = {r["mag_id"]: r for r in rows("30_dbcan/family_matrix.tsv")}
    out = [[m, fm[m]["total"]] + [cnt[m].get(f, 0) for f in fams] + ["", ""] for m in MAGS]
    write("T_cazy_family_matrix.tsv", ["mag_id", "total_CAZymes"] + fams + ["PUL_count", "T9SS_cargo_count"], out)
    note("T_cazy_family_matrix", "PUL / T9SS", "no PUL or T9SS table exists under 30_dbcan (only gaoshiqia_nonpul.tsv for one MAG); NA")

# ================================================================== T10 pyrrolysine
def t10():
    mk = rows("34b_methylamine/markers.tsv")
    out = [[r["ko"], r["gene"], r["role"], r["source"], r["contig"], r["segment"], r["score"], r["in_frame_UAG"],
            r["internal_UAG"], r["rel_positions"]] for r in mk]
    write("T_pyrrolysine.tsv", ["KO", "gene", "class", "source", "contig", "segment", "hmm_score_bits", "in_frame_UAG",
                                "internal_UAG", "UAG_relative_positions"], out)
    # class summary straight from VERDICT.txt
    t = open(P("34b_methylamine/VERDICT.txt")).read()
    summ = re.findall(r"^\s+(PYL|NONPYL_CTRL|CORRINOID_DECOY|PYL_MACHINERY)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)%\s+([\d\.]+)$", t, re.M)
    fisher = re.search(r"Fisher exact, one-sided\s*:\s*p = (\S+)", t).group(1)
    write("T_pyrrolysine_summary.tsv", ["class", "segments", "any_UAG", "terminal_UAG", "internal_UAG", "pct_internal", "internal_per_1000_aa", "Fisher_one_sided_p_PYL_vs_controls"],
          [list(s) + [fisher if s[0] == "PYL" else ""] for s in summ])

# ================================================================== T11 read-level markers
def t11():
    out = []
    v = open(P("31_readlevel/VERDICT.txt")).read()
    for r in rows("79_mtbB_mtaB/counts.tsv"):
        out.append(["79_mtbB_mtaB/counts.tsv", r["label"], r["SED"], r["LIQ"], ""])
    for r in rows("34_unbinned/marker_levels.tsv"):
        out.append(["34_unbinned/marker_levels.tsv", r["gene"] + " (" + r["ko"] + ")", r["reads_sed"], r["reads_liq"],
                    "sed/rpoB=%s; level SED=%s; level LIQ=%s" % (r["sed_over_rpoB"], r["level"], r["level_liq"])])
    for r in rows("80_mtbB_finish/ko_counts.tsv", header=False):
        out.append(["80_mtbB_finish/ko_counts.tsv", r[1] + " (" + r[0] + ") reference-panel size", "", "", "panel sequences: " + r[2]])
    write("T_readlevel_markers.tsv", ["source", "marker", "reads_SED", "reads_LIQ", "panel_composition_or_ratio"], out)
    rt = rows("82_mtbB_who/reference_taxonomy.tsv")
    write("T_readlevel_mtbB_besthit.tsv", list(rt[0].keys()), [list(r.values()) for r in rt])
    note("T_readlevel_markers", "83_who_really/reference_taxonomy.tsv", "file does not exist; the only reference_taxonomy.tsv is 82_mtbB_who/ (used). 83_who_really holds own_tax.tsv / ext_tax.tsv instead")
    note("T_readlevel_markers", "EF-1alpha, tuf, glnA, mtmB(61/63)", "only markers present in the listed tables are exported; 61_mtmB/63 hit files are per-read and were not re-scored")

# ================================================================== T12 auxotrophy
def t12():
    pc = rows("19b_auxotrophy/pathway_completeness.tsv")
    aux = {(r["mag_id"], r["pathway"]) for r in rows("19b_auxotrophy/auxotrophies.tsv")}
    cys = {mid(r["mag"]): r for r in rows("74_cys/cys_panel.tsv")}
    paths = [c for c in pc[0] if c not in ("mag_id", "working_name", "sample", "domain", "genus", "completeness")]
    out = []
    for r in pc:
        for p in paths:
            call = "auxotroph" if (r["mag_id"], p) in aux else "prototroph"
            flag = ""
            if p == "Cys" and r["domain"] == "Archaea":
                flag = "ARTEFACT: archaeal Cys via SepRS/SepCysS (74_cys)"
                c = cys.get(r["mag_id"])
                if c: flag += "; SepRS K07587 pass=%s, SepCysS K06868 pass=%s" % (c["K07587_pass"], c["K06868_pass"])
            out.append([r["mag_id"], r["domain"], r["genus"], r["completeness"], p, r[p], call, flag])
    write("T_auxotrophy_full.tsv", ["mag_id", "domain", "genus", "CheckM2_completeness_pct", "pathway", "route_completeness_fraction", "call", "flag"], out)

# ================================================================== T13 B12 / ECF
def t13():
    js = json.load(open(P("45_kofam_full/b12_genomewide.json")))
    ecf = {mid(r["mag"]): r for r in rows("36_gapclose/ecf/ecf_matrix.tsv")}
    out = []
    for g in js:
        m = mid(g["genome"]); e = ecf.get(m, {})
        out.append([m, len(g.get("prod", [])), ",".join(g.get("prod", [])), len(g.get("cons", [])), ",".join(g.get("cons", [])),
                    ",".join(g.get("tran", [])), e.get("module_kos"), e.get("scomponent_kos"), e.get("complete")])
    write("T_b12_ecf.tsv", ["mag_id", "n_B12_biosynthesis_KOs", "B12_biosynthesis_KOs", "n_B12_dependent_enzyme_KOs",
                            "B12_dependent_enzyme_KOs", "B12_transport_KOs_BtuCDF_etc", "ECF_module_KOs", "ECF_S_component_KOs", "ECF_complete"], sorted(out))
    note("T_b12_ecf", "source path", "b12_genomewide.json is at 45_kofam_full/, not 45_kofam_full/tables/ (single candidate)")

# ================================================================== T14 SMETANA
def t14():
    pr = rows("43_smetana/run_43q/pairs/pairs_detailed.tsv")
    vi = json.load(open(P("43_smetana/run_43q/viable_pairs.json")))
    cls = {}
    cn = os.path.join(WS2, "67_network/compound_names.tsv")
    if os.path.exists(cn):
        for r in rows("67_network/compound_names.tsv", header=False)[1:]:
            cls[r[0]] = r[-1]
    keys = list(pr[0].keys())
    out = []
    for r in pr:
        vals = [mid(r[k]) if k in ("donor", "receiver") else r[k] for k in keys]
        out.append(vals)
    write("T_smetana_pairs.tsv", keys, out)
    write("T_smetana_viable_pairs.tsv", ["viable_pairs_json"], [[json.dumps(vi)]])
    note("T_smetana_pairs", "compound class / positive controls", "exported as stored in pairs_detailed.tsv; class column only if present there; viable pairs exported verbatim from viable_pairs.json")

# ================================================================== T15 env stress
def t15():
    em = rows("16_environment/env_matrix.tsv")
    tg = {r[0]: r for r in rows("16_environment/targets.tsv", header=False)}
    kos = [c for c in em[0] if re.fullmatch(r"K\d{5}", c)]
    out = []
    for k in kos:
        t = tg.get(k, [k, "", "", "", ""])
        out.append([k, t[1], t[2], t[3], t[4] if len(t) > 4 else ""] + [r[k] for r in em])
    write("T_env_stress.tsv", ["KO", "gene", "category", "subcategory", "definition"] + [mid(r["mag"]) for r in em], out)
    note("T_env_stress", "threshold", "per-KO thresholds are KOfam ko_list values; the step-16 ko_list is not stored in 16_environment, so the column is omitted (NA)")

# ================================================================== T16 mobilome
def t16():
    write("T_mobilome_assembly.tsv", list(rows("18_viruses/assembly_virus_summary.tsv")[0].keys()),
          [list(r.values()) for r in rows("18_viruses/assembly_virus_summary.tsv")])
    pv = rows("18_viruses/provirus_host_links.tsv")
    out = []
    for r in pv:
        sn = r["seq_name"].split("~")[-1]
        m_ = re.search(r"\|provirus_(\d+)_(\d+)", sn)
        out.append([r["mag_id"], r["genus"], sn.split("|")[0], m_.group(1) if m_ else "", m_.group(2) if m_ else "",
                    r["length"], r["virus_score"], r["n_hallmarks"], r["taxonomy"]])
    write("T_mobilome_proviruses.tsv", ["mag_id", "genus", "contig", "start", "end", "length_bp", "geNomad_virus_score", "n_hallmarks", "taxonomy"], out)
    ds = rows("28_crispr/defence_systems.tsv")
    write("T_mobilome_defence.tsv", list(ds[0].keys()), [list(r.values()) for r in ds])
    hv = rows("28_crispr/hits_verified.tsv")
    write("T_mobilome_crispr_hits.tsv", list(hv[0].keys()), [list(r.values()) for r in hv])
    note("T_mobilome", "split", "exported as 4 TSV (assembly, proviruses, defence, CRISPR spacer hits with flanking-k-mer control verdict); defence-system tool name is not recorded in defence_systems.tsv")

# ================================================================== T17 antiSMASH
def t17():
    mh = rows("33_antismash/mibig_hits.tsv")
    bg = rows("33_antismash/bgc_regions.tsv")
    # antiSMASH numbers regions per contig in start order; join on (MAG, contig, region)
    byc = defaultdict(list)
    for b in bg: byc[(mid(b["genome"]), b["contig"])].append(b)
    for k in byc: byc[k].sort(key=lambda b: int(b["start"]))
    out = []
    for r in mh:
        b = byc[(r["mag_id"], r["contig"])][int(r["region"]) - 1]
        if set(b["products"].split(";")) != set(r["products"].split(";")):
            note("T_antismash", "%s %s r%s" % (r["mag_id"], r["contig"], r["region"]), "product mismatch on join: %s vs %s" % (b["products"], r["products"]))
        out.append([r["mag_id"], r["contig"], r["region"], b["start"], b["end"], r["products"], r["best_accession"],
                    r["best_compound"], r["best_type"], r["genes_matched"], r["query_genes"],
                    "" if r["frac_matched"] in ("", "NA") else round(100 * float(r["frac_matched"]), 1)])
    write("T_antismash.tsv", ["mag_id", "contig", "region", "start", "end", "class", "MIBiG_best_hit", "MIBiG_compound",
                              "MIBiG_type", "genes_matched", "query_genes", "pct_genes_matched"], out)

# ================================================================== T18 strain
def t18():
    fa = rows("32_strain/fixed_sites_annotated.tsv")
    write("T_strain_divergence.tsv", list(fa[0].keys()), [list(r.values()) for r in fa])
    note("T_strain_divergence", "summary numbers", "per-site table exported; sites compared / per kb / Ti/Tv summary are in 32_strain/DIVERGENCE.txt and ANATOMY.txt (copied to E2_tables/text/)")

# ================================================================== T19 profilers
def t19():
    out = []
    for lib in ("SED", "LIQ"):
        rep = rows("13_profile/kraken/%s_reads.report" % lib, header=False)
        unc = next((r for r in rep if r[3].strip() == "U"), None)
        tot = int(unc[1]) + sum(int(r[1]) for r in rep if r[3].strip() == "R")
        byrank = Counter()
        for r in rep:
            byrank[r[3].strip()[0]] += int(r[2])
        out.append(["Kraken2", lib, "classified_pct", round(100 - float(unc[0]), 2)])
        for k, n in (("D", "domain"), ("P", "phylum"), ("C", "class"), ("O", "order"), ("F", "family"), ("G", "genus"), ("S", "species")):
            out.append(["Kraken2", lib, "reads_assigned_at_%s" % n, byrank.get(k, 0)])
        t = open(P("35_kaiju/ASSIGNMENT_RANK.txt")).read()
        sec = re.search(lib + r" — RANK OF THE TAXON.*?classified\s*:\s*([\d,]+)\s*\(([\d\.]+)%\)(.*?)\n={10}", t, re.S)
        if sec:
            out.append(["Kaiju", lib, "classified_pct", sec.group(2)])
            for m_ in re.finditer(r"^\s{2}(\w[\w ]*?)\s{2,}([\d,]+)\s+([\d\.]+)%\s+([\d\.]+)%", sec.group(3), re.M):
                out.append(["Kaiju", lib, "reads_assigned_at_" + m_.group(1).strip().replace(" ", "_"), m_.group(2).replace(",", "")])
    for lib in ("SED", "LIQ"):
        prof = rows("22_singlem/%s_profile.tsv" % lib)
        depth = Counter()
        for r in prof:
            last = r["taxonomy"].split("; ")[-1][:3]
            depth[last] += float(r["coverage"])
        for k in ("d__", "p__", "c__", "o__", "f__", "g__", "s__"):
            out.append(["SingleM", lib, "coverage_with_deepest_rank_" + k.rstrip("_"), round(depth.get(k, 0), 2)])
        for r in prof:
            if "o__Methanomassiliicoccales" in r["taxonomy"]:
                out.append(["SingleM", lib, "Methanomassiliicoccales: " + r["taxonomy"].split("; ")[-1], r["coverage"]])
    write("T_profilers.tsv", ["profiler", "library", "metric", "value"], out)
    note("T_profilers", "SingleM filled phylum table / % of filled coverage", "singlem_filled_phylum.py not on the server; not computed (see E4 #2)")
    note("T_profilers", "SingleM unnamed fractions", "no stored definition of 'unnamed'; not computed")

# ================================================================== T20, T21
def t20():
    out = []
    for lib in ("SED", "LIQ"):
        for r in rows("27_rrna_corrected/%s_corrected_profile.tsv" % lib):
            out.append([lib, r["genus"], r["source_mag"], r["operons"], r["raw_reads"], r["raw_pct_of_subset"],
                        r["corrected_pct_of_subset"], r["change_pct_points"]])
    write("T_rrna_copy_correction.tsv", ["library", "genus", "source_mag", "rRNA_operons", "raw_reads", "raw_pct_of_subset",
                                         "corrected_pct_of_subset", "change_pct_points"], out)

def t21():
    ms = {r["bam"]: r for r in rows("02_mapping/mapping_summary.tsv")}
    out = []
    for lib in ("SED", "LIQ"):
        fs = open(P("02_mapping/%s_map_self.flagstat.txt" % lib)).read()
        prim = re.search(r"primary mapped \(([\d\.]+)%", fs).group(1)
        asm = sum(int(r[1]) for r in rows("04_polish/%s_asm_final.fasta.%s" % (lib, "fai" if lib == "LIQ" else "seqkit.fai"), header=False))
        mag = sum(int(c["bp"]) for c in CW if c["sample"] == lib)
        out.append([lib, ms["%s_map_self" % lib]["pct_mapped"], prim, asm, mag, round(100 * mag / asm, 2)])
    write("T_unbinned_fraction.tsv", ["library", "pct_reads_mapped_to_assembly_mapping_summary", "pct_primary_reads_mapped_flagstat",
                                      "assembly_bp", "bp_in_MAGs", "pct_assembly_bases_in_MAGs"], out)
    ml = rows("34_unbinned/marker_levels.tsv")
    write("T_unbinned_marker_levels.tsv", list(ml[0].keys()), [list(r.values()) for r in ml])
    note("T_unbinned_fraction", "% reads mapped to MAGs", "not stored per library in any table; not computed (needs a re-mapping)")

if __name__ == "__main__":
    os.makedirs(os.path.join(OUT, "text"), exist_ok=True)
    for fn in (t1, t2, t3, t4, t5, t6_t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21):
        try:
            fn()
        except Exception as e:
            print("  !! %s failed: %r" % (fn.__name__, e)); NOTES.append((fn.__name__, "FAILED", repr(e)))
    for f in ("32_strain/DIVERGENCE.txt", "32_strain/ANATOMY.txt", "34b_methylamine/VERDICT.txt", "31_readlevel/VERDICT.txt", "35_kaiju/ASSIGNMENT_RANK.txt"):
        import shutil; shutil.copy(P(f), os.path.join(OUT, "text", f.replace("/", "__")))
    write("E2_notes.tsv", ["table", "item", "note"], NOTES)
    with open(os.path.join(M, "inputs_md5_E2.tsv"), "w") as f:
        f.write("part\tpath\tmd5\n")
        for p, h in INPUTS.items(): f.write("E2\t%s\t%s\n" % (p, h))
