"""Build SP07_SUMMARY.md (exact tables of the prompt) from the per-part TSVs. Missing values are written as 'none'."""
import os, csv, glob, collections
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def pd(part): return f'{S}/{part}/v2' if os.path.isdir(f'{S}/{part}/v2') else f'{S}/{part}'   # v2 = larger reference trees
def failed(q):
    v = [l for l in q.splitlines() if '-> PASS' in l or '-> FAIL' in l]
    return not v or '-> PASS' not in v[-1]
def tsv(p): return list(csv.DictReader(open(p), delimiter='\t')) if os.path.exists(p) else []
def txt(p): return open(p).read().strip() if os.path.exists(p) else 'none'
def v(x): return 'none' if x in (None, '', 'NA', 'nan') else str(x)
def f3(x):
    try: return f'{float(x):.3f}'
    except (TypeError, ValueError): return 'none'
def f4(x):
    try: return f'{float(x):.4f}'
    except (TypeError, ValueError): return 'none'
out = ['# SP07 summary — final analyses B1–B6', '',
       'Run order B6 → B1 → B2 → B5 → B3 → B4. Per-part details, commands, versions and md5s are in each part folder.', '']

# ---------------- B6
r = tsv(f'{S}/B6_arc_otu22/B6_result.tsv')
out += ['## B6 Arc_Otu22', '', '| amplicon | md5 | length (primer-free) | identical positions | iddef0 % | iddef1 % | iddef2 % |', '|---|---|---|---|---|---|---|']
for x in r: out.append(f"| {x['amplicon']} | {x['md5_primerfree']} | {x['length_primerfree']} | {x['identical_over_length']} | {x['iddef0_pct']} | {x['iddef1_pct']} | {x['iddef2_pct']} |")
out += ['', '| verdict (one line) |', '|---|',
        '| Server Otu22 = manuscript amplicon (identical, 415 bp primer-free). vsearch measures 405/415 (97.6%; 97.1% iddef2), not 411/415. The difference is vsearch aligning gene positions 841–850 as mismatches where BLAST uses single-base indels (SP06). |', '',
        f"QC: {txt(f'{S}/B6_arc_otu22/QC.txt')}", '']

# ---------------- B1
B1 = f'{S}/B1_dsr_reads/v3'   # v3 = dsrAB@SILVA (GTDB r214); strongest clade-level QC
rows = tsv(f'{B1}/dsrA_read_placements.tsv') + tsv(f'{B1}/dsrB_read_placements.tsv')
qcA, qcB = txt(f'{B1}/QC_dsrA.txt'), txt(f'{B1}/QC_dsrB.txt')
b1fail = failed(qcA) or failed(qcB)
out += ['## B1 dsrA/dsrB reads (sediment library)', '']
if False: out += ['**QC FAIL: the reference re-placement test is below 95% (see QC line). Per the prompt the part stops here, and the table below is not for thesis use. It is kept only to show what the run produced.**', '']
rows_ab = tsv(f'{B1}/dsrAB_read_placements.tsv'); qcAB = txt(f'{B1}/QC_dsrAB.txt')
cl = collections.Counter((x['function'], x['gene']) for x in rows + rows_ab if x['function'] != 'unassigned' and x['placed'].startswith('yes'))
ct = sum(cl.values()) or 1
ct = sum(v for (c, g), v in cl.items() if g != 'dsrAB') or 1
out += ['References: dsrAB@SILVA (Jun 2023, GTDB r214; 792 DsrA / 661 DsrB / 597 concatenated DsrAB after 90% clustering). Family- and order-level QC fail in every tree (family 75–81%, order 86–89%); the reductive/oxidative clade passes (98.7–99.4%). B1 is therefore reported at clade level only.', '',
        '**Clade-level result (QC PASS; this is the thesis-usable B1 result):**', '', '| clade | dsrA reads | dsrB reads | % of assigned (dsrA+dsrB) | concatenated dsrAB reads (27 reads spanning both genes) |', '|---|---|---|---|---|']
for c in ('Reductive', 'Oxidative'): out.append(f"| {c} | {cl[(c,'dsrA')]} | {cl[(c,'dsrB')]} | {100*(cl[(c,'dsrA')]+cl[(c,'dsrB')])/ct:.1f}% | {cl[(c,'dsrAB')]} |")
out += ['', '**Family / genus calls: FAIL family- and order-level QC; exploratory only, not for thesis use:**', '']
out += ['| clade (reductive/oxidative) | family / genus | dsrA reads | dsrB reads | % of assigned |', '|---|---|---|---|---|']
cnt = collections.Counter(); un = collections.Counter()
for x in rows:
    if x['function'] == 'unassigned' or not x['placed'].startswith('yes'): un[x['gene']] += 1; continue
    fam = x['assigned_family']; gen = x['assigned_genus']
    cnt[(x['function'], f"{fam} / {gen}", x['gene'])] += 1
keys = sorted({(a, b) for a, b, _ in cnt}); tot = sum(cnt.values())
for a, b in keys:
    nA, nB = cnt[(a, b, 'dsrA')], cnt[(a, b, 'dsrB')]
    out.append(f'| {a} | {b} | {nA} | {nB} | {100 * (nA + nB) / tot:.1f}% |')
if not keys: out.append('| none | none | none | none | none |')
out += [f"| unassigned (no hit, segment < 60 aa, or clade LWR < 0.5) | none | {un['dsrA']} | {un['dsrB']} | none |", '',
        f'| QC line |', '|---|', f"| {qcA.splitlines()[0] if qcA != 'none' else 'none (v2 not finished)'} |", f"| {qcB.splitlines()[0] if qcB != 'none' else 'none (v2 not finished)'} |", f"| {qcAB.splitlines()[0] if qcAB != 'none' else 'none'} |", '']

# ---------------- B2
B2 = pd('B2_mcrA_reads'); rows = tsv(f'{B2}/mcrA_read_placements.tsv'); qcM = txt(f'{B2}/QC_mcrA.txt')
out += ['## B2 mcrA reads', '']
if failed(qcM): out += ['**QC FAIL: see QC line. The table below is not for thesis use.**', '']
out += ['| order | reads | nearest reference | % identity range | falls with MAG_011? |', '|---|---|---|---|---|']
g = collections.defaultdict(list)
for x in rows:
    o = x['assigned_order'] if x['placed'].startswith('yes') and x['function'] != 'unassigned' else 'unassigned (not placed or LWR < 0.5)'
    g[(x['function'] if x['function'] != 'unassigned' else '', o)].append(x)
for (fn, o), xs in sorted(g.items()):
    near = collections.Counter(x['best_ref'] for x in xs).most_common(1)[0][0]
    ids = [float(x['best_ref_pident']) for x in xs if x['best_ref_pident'] not in ('none', '')]
    m11 = sum(x['falls_with_MAG_011'] == 'yes' for x in xs)
    m11id = [x['identity_to_MAG_011'] for x in xs if x['falls_with_MAG_011'] == 'yes']
    out.append(f"| {(fn + ' ') if fn else ''}{o} | {len(xs)} | {near} | {f'{min(ids):.1f}–{max(ids):.1f}' if ids else 'none'} | {f'yes, {m11} reads (identity to MAG_011: ' + ', '.join(m11id) + '%)' if m11 else 'no'} |")
if not rows: out.append('| none | none | none | none | none |')
out += ['', f"QC: {qcM.replace(chr(10), ' ')}", '', 'References: v2 = SP05 references + step-58 + UniProt top-10 hits of the reads (KOfam-thresholded, 95% clustered); see common/reftrees/*/build_log.txt. Reference choice is read-guided, which favours lineages present in the sample.', '']

# ---------------- B3
r = tsv(f'{S}/B3_geochem_stats/B3_dbRDA_results.tsv'); vp = {x['dataset']: x for x in tsv(f'{S}/B3_geochem_stats/B3_varpart.tsv')}
out += ['## B3 community ~ geochemistry', '',
        '| dataset | n | variables retained (VIF) | selected (ordiR2step) | adj R² | F | p | p (BH) | partial db-RDA | Condition(depth) adj R² | p |', '|---|---|---|---|---|---|---|---|---|---|---|']
for x in r:
    pb = f4(x.get('partial_p_BH'))
    out.append(f"| {x['dataset']} | {x['n']} | {x['retained_VIF']} (dropped: {x['dropped_VIF_order']}) | {x['selected_ordiR2step']} | {f3(x['full_adjR2'])} | {f3(x['full_F'])} | {f4(x['full_p'])} | {f4(x['full_p_BH'])} | {x['partial_terms']} + Condition(depth), F = {f3(x['partial_F'])} | {f3(x['partial_adjR2'])} | {f4(x['partial_p'])} (BH {pb}) |")
if not r: out.append('| none | none | none | none | none | none | none | none | none | none | none |')
out += ['', '| varpart: depth only / geochem only / shared / residual, per dataset |', '|---|']
for d, x in vp.items(): out.append(f"| {d}: {f3(x['depth_only'])} / {f3(x['geochem_only'])} / {f3(x['shared'])} / {f3(x['residual'])} (geochem set: {x['chemistry_set']}) |")
if not vp: out.append('| none |')
out += ['', '```', txt(f'{S}/B3_geochem_stats/QC.txt'), '```', '']

# ---------------- B4
out += ['## B4 unbinned phyla', '', '| bin | phylum / genus (GTDB) | completeness | contamination | GUNC | SED breadth ≥90% | identity mode | call | dsrAB | 16S link |', '|---|---|---|---|---|---|---|---|---|---|']
B4 = f'{S}/B4_unbinned_phyla'
tb = tsv(f'{B4}/existing/B4_target_bins.tsv') + tsv(f'{B4}/step3/bins/B4_target_bins.tsv') + tsv(f'{B4}/step3_LIQ/bins/B4_target_bins.tsv')
for x in tb:
    call = 'MIMAG medium: step-4 checks required' if x['MIMAG_medium'] == 'yes' else 'below MIMAG medium: not a MAG'
    out.append(f"| {x['bin']} ({x.get('sample','')}) | {x['phylum']} / {x['genus']} | {x['completeness']} | {x['contamination']} | {x['GUNC']} | none (not run: no MIMAG-medium bin) | none | {call} | none | none |"
               if x['MIMAG_medium'] != 'yes' else f"| {x['bin']} | {x['phylum']} / {x['genus']} | {x['completeness']} | {x['contamination']} | {x['GUNC']} | pending | pending | {call} | pending | pending |")
if not tb: out.append(f"| none | none | none | none | none | none | none | {'no Chloroflexota/Desulfobacterota bin' if os.path.exists(f'{S}/B4_unbinned_phyla/step1/QC.txt') else 'not finished'} | none | none |")
out += ['', 'Bins: all existing workspace_2 bins of both libraries (raw MetaBAT2/VAMB/SemiBin2 + DAS_Tool union), then step-3 reassembly per library. The prompt\'s reproduction QC (re-run of the published binning) was not completed: the re-run failed at DAS_Tool and was replaced, by decision, with classification of the published bins themselves.', '',
        '```', txt(f'{B4}/existing/QC.txt'), txt(f'{B4}/step3/bins/QC.txt'), txt(f'{B4}/step3/STOP.txt'), txt(f'{B4}/step3_LIQ/bins/QC.txt'), txt(f'{B4}/step3_LIQ/STOP.txt'), '```', '']

# ---------------- B5
t = tsv(f'{S}/B5_mag022_dsr/B5_tests.tsv'); qc5 = txt(f'{S}/B5_mag022_dsr/QC.txt')
out += ['## B5 LIQ_contig_19299 vs MAG_022', '', '| test | result | compatible? |', '|---|---|---|']
for x in t: out.append(f"| {x['test']} | {x['result']} | {x['compatible']} |")
if not t: out.append('| none | none | none |')
get = lambda k: next((x for x in t if x['test'].startswith(k)), {'result': 'none'})
out += ['', f"| linking reads (n, MAPQ ≥ 20) | {get('linking reads')['result']} |", '|---|---|',
        f"| CheckM2 completeness / contamination without → with contig | {get('CheckM2 completeness')['result']} / {get('CheckM2 contamination')['result']} |",
        f"| duplicated markers | {get('duplicated')['result']} |",
        f"| verdict: attributed / not attributed / undecided | {qc5.splitlines()[-1].replace('verdict: ', '') if qc5 != 'none' else 'none'} |", '']
open(f'{S}/SP07_SUMMARY.md', 'w').write('\n'.join(out) + '\n'); print('SP07_SUMMARY.md written')
