"""Full catalogue of every workspace_2 bin (SED + LIQ; raw MetaBAT2/VAMB/SemiBin2 + DAS_Tool union) with assembly statistics,
CheckM2, GUNC, GTDB-Tk (R232) and the final-MAG crosswalk, plus a two-panel thesis figure.
Outputs: B4_unbinned_phyla/catalogue/ALL_BINS_catalogue.tsv (+ .xlsx), Fig_all_bins.pdf/.png, catalogue_notes.txt"""
import os, csv, glob, re, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); E = f'{P}/existing'; O = f'{P}/catalogue'; os.makedirs(O, exist_ok=True)
WS2 = os.path.abspath(f'{P}/../../..')
def tsv(p): return list(csv.DictReader(open(p), delimiter='\t'))
if not glob.glob(f'{E}/gtdbtk/gtdbtk.bac120.summary.tsv'): raise SystemExit('bac120 summary missing: catalogue not built')
cm = {r['Name']: r for r in tsv(f'{E}/checkm2/quality_report.tsv')}
gu = {r['genome']: r for p in glob.glob(f'{E}/gunc/GUNC*maxCSS*.tsv') for r in tsv(p)}
gt = {r['user_genome']: r for p in glob.glob(f'{E}/gtdbtk/gtdbtk.*.summary.tsv') for r in tsv(p)}
cw = {r['qc_bin_name']: r for r in tsv(f'{WS2}/07_mags/MAG_CROSSWALK.tsv')}
cw_raw = {r['raw_binner_file'][:-3]: r for r in tsv(f'{WS2}/07_mags/MAG_CROSSWALK.tsv')}   # raw bin a final MAG was refined from
def stats(p):
    L, gc, n = [], 0, 0
    s = []
    for l in open(p):
        if l.startswith('>'):
            if s: x = ''.join(s); L.append(len(x)); gc += x.count('G') + x.count('C'); s = []
        else: s.append(l.strip().upper())
    if s: x = ''.join(s); L.append(len(x)); gc += x.count('G') + x.count('C')
    L.sort(reverse=True); tot = sum(L); acc = 0; n50 = 0
    for x in L:
        acc += x
        if acc >= tot / 2: n50 = x; break
    return len(L), tot, n50, L[0] if L else 0, gc / tot if tot else 0
rank = lambda c, k: next((x[3:] for x in c.split(';') if x.startswith(k + '__')), '') if c else ''
rows = []
for f in sorted(glob.glob(f'{E}/bins/*.fa')):
    b = os.path.basename(f)[:-3]; m = re.match(r'(SED|LIQ)_(das_)?(metabat|vamb|semibin|unknown)_(\w+)', b)
    nc, bp, n50, mx, gc = stats(f); q = cm.get(b, {}); g = gt.get(b, {}); c = g.get('classification', '')
    comp = float(q.get('Completeness', 'nan')); cont = float(q.get('Contamination', 'nan'))
    tier = 'high' if comp > 90 and cont < 5 else ('medium' if comp >= 50 and cont < 10 else 'low')
    mag = cw.get(b, {}).get('mag_id', '') or (f"source of {cw_raw[b]['mag_id']}" if b in cw_raw else '')
    rows.append(dict(bin=b, library=m.group(1), stage='DAS_Tool' if m.group(2) else 'raw', binner=m.group(3), final_MAG=mag or 'no',
        n_contigs=nc, size_bp=bp, N50_bp=n50, longest_contig_bp=mx, GC=round(gc, 4), completeness=comp, contamination=cont,
        checkm2_model=q.get('Completeness_Model_Used', ''), MIMAG_tier_checkm2=tier,
        MIMAG_final=cw.get(b, {}).get('mimag', '') or 'n/a (not a final MAG)',
        GUNC_CSS=gu.get(b, {}).get('clade_separation_score', 'NA'), GUNC_pass=gu.get(b, {}).get('pass.GUNC', 'NA'),
        domain=rank(c, 'd') or 'unclassified', phylum=rank(c, 'p') or 'unclassified', class_=rank(c, 'c'), order=rank(c, 'o'),
        family=rank(c, 'f'), genus=rank(c, 'g'), species=rank(c, 's'), gtdb_classification=c or 'unclassified',
        closest_reference=g.get('closest_genome_reference', g.get('fastani_reference', '')) or 'none',
        closest_ANI=g.get('closest_genome_ani', g.get('fastani_ani', '')) or 'none', red_value=g.get('red_value', '') or 'none'))
cols = list(rows[0])
with open(f'{O}/ALL_BINS_catalogue.tsv', 'w') as o:
    o.write('\t'.join(c.rstrip('_') for c in cols) + '\n')
    for r in rows: o.write('\t'.join(str(r[c]) for c in cols) + '\n')
try:
    import openpyxl
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = 'all_bins'; ws.append([c.rstrip('_') for c in cols])
    for r in rows: ws.append([r[c] for c in cols])
    ws.freeze_panes = 'B2'; wb.save(f'{O}/ALL_BINS_catalogue.xlsx')
except Exception as e: print('xlsx skipped:', e)

# ------------------------------------------------------------------ figure
INK, MUTED, GRID = '#1f1f1e', '#6b6a63', '#e4e3dc'
C_MAG, C_TGT, C_OTH = '#2a78d6', '#eb6834', '#b5b3a8'   # categorical slots 1-2 + neutral
TIER = {'high': '#1c5cab', 'medium': '#6da7ec', 'low': '#cde2fb'}   # one-hue sequential ramp
TGT = ('Chloroflexota', 'Desulfobacterota')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False})
def plot(rows, TITLE, FN, MAGLAB):
    fig = plt.figure(figsize=(7.2, 9.0))
    YMAX = max(10, max(r['contamination'] for r in rows) * 1.15)
    XMAXB = 0
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], hspace=0.55, wspace=0.28)
    for j, lib in enumerate(('SED', 'LIQ')):
        ax = fig.add_subplot(gs[0, j]); R = [r for r in rows if r['library'] == lib]
        for cat, col, z, lab in (('oth', C_OTH, 1, 'other bins'), ('tgt', C_TGT, 3, 'Chloroflexota / Desulfobacterota'), ('mag', C_MAG, 4, MAGLAB)):
            pts = [r for r in R if (r['final_MAG'] != 'no' if cat == 'mag' else (r['phylum'].startswith(TGT) and r['final_MAG'] == 'no' if cat == 'tgt'
                   else r['final_MAG'] == 'no' and not r['phylum'].startswith(TGT)))]
            ax.scatter([r['completeness'] for r in pts], [r['contamination'] for r in pts], s=18, c=col, edgecolors='white', linewidths=0.6, zorder=z, label=lab)
        ax.axvline(50, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=0); ax.axvline(90, color=MUTED, lw=0.8, ls=(0, (1, 2)), zorder=0)
        ax.axhline(10, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=0); ax.axhline(5, color=MUTED, lw=0.8, ls=(0, (1, 2)), zorder=0)
        ax.set_yscale('symlog', linthresh=10); ax.set_xlim(0, 101); ax.set_ylim(-0.5, YMAX)
        ax.set_yticks([t for t in (0, 2, 5, 10, 20, 50, 100, 200) if t <= YMAX]); ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'))
        ax.set_xlabel('CheckM2 completeness (%)'); ax.set_ylabel('CheckM2 contamination (%)' if j == 0 else '')
        ax.grid(axis='y', color=GRID, lw=0.6); ax.set_axisbelow(True)
        nm = sum(r['MIMAG_tier_checkm2'] != 'low' for r in R)
        ax.set_title(f"{'Sediment' if lib == 'SED' else 'Liquid cultures'}: {len(R)} bins, {nm} ≥ MIMAG medium", loc='left', fontsize=8.5, color=INK)
        if j == 0: hA, lA = ax.get_legend_handles_labels()
    # panel B: bins per phylum x quality tier, per library (top 14 phyla, rest folded to Other)
    cnt = collections.Counter(r['phylum'] for r in rows if r['phylum'] != 'unclassified'); top = [p for p, _ in cnt.most_common(14)]
    phy = lambda r: 'not classifiable (too few markers)' if r['phylum'] == 'unclassified' else (r['phylum'] if r['phylum'] in top else 'Other')
    order = sorted(set(map(phy, rows)), key=lambda p: (p == 'not classifiable (too few markers)', p == 'Other', -sum(phy(r) == p for r in rows)))
    for j, lib in enumerate(('SED', 'LIQ')):
        ax = fig.add_subplot(gs[1, j], ); R = [r for r in rows if r['library'] == lib]
        if j == 0: axB0 = ax
        y = list(range(len(order)))[::-1]; left = [0] * len(order)
        for t in ('high', 'medium', 'low'):
            v = [sum(1 for r in R if phy(r) == p and r['MIMAG_tier_checkm2'] == t) for p in order]
            ax.barh(y, v, left=left, color=TIER[t], edgecolor='white', linewidth=1, height=0.72, label=f'{t} quality')
            left = [a + b for a, b in zip(left, v)]
        for yi, p, tot in zip(y, order, left):
            if tot: ax.text(tot + max(left) * 0.015, yi, str(tot), va='center', fontsize=6.5, color=MUTED)
        ax.set_yticks(y); ax.set_yticklabels(order if j == 0 else [''] * len(order), fontsize=7)
        for lab, p in zip(ax.get_yticklabels(), order):
            if p.startswith(TGT): lab.set_color(C_TGT); lab.set_fontweight('bold')
        ax.set_xlabel('number of bins'); ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.set_title(f"{'Sediment' if lib == 'SED' else 'Liquid cultures'}: bins per GTDB phylum", loc='left', fontsize=8.5, color=INK)
        if j == 0: hB, lB = ax.get_legend_handles_labels()
    fig.text(0.01, 0.93, 'A', fontsize=12, fontweight='bold'); fig.text(0.01, 0.50, 'B', fontsize=12, fontweight='bold')
    fig.text(0.06, 0.935, TITLE, fontsize=9.5, color=INK)
    fig.legend(hA, lA, loc='center', bbox_to_anchor=(0.5, 0.535), ncol=3, frameon=False, fontsize=7, handletextpad=0.2, columnspacing=1.2, title='lines: dashed = MIMAG medium (≥50% / <10%), dotted = high (>90% / <5%)', title_fontsize=7)
    fig.legend(hB, lB, loc='center', bbox_to_anchor=(0.5, 0.045), ncol=3, frameon=False, fontsize=7, title='CheckM2 MIMAG tier', title_fontsize=7)
    for ext in ('pdf', 'png'): fig.savefig(f'{O}/{FN}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)

plot(rows, 'All workspace_2 bins (raw MetaBAT2 / VAMB / SemiBin2 + DAS_Tool union)', 'Fig_all_bins', 'final MAG or its source bin')
for bn, lab in (('metabat', 'MetaBAT2'), ('vamb', 'VAMB'), ('semibin', 'SemiBin2')):
    plot([r for r in rows if r['binner'] == bn and r['stage'] == 'raw'], f'{lab} raw bins (both libraries)', f'Fig_bins_{lab}', f'source bin of a final MAG')
n = collections.Counter((r['library'], r['stage']) for r in rows)
open(f'{O}/catalogue_notes.txt', 'w').write(
    f'{len(rows)} bins: ' + ', '.join(f'{a} {b} {v}' for (a, b), v in sorted(n.items())) + '\n'
    'MIMAG tier from CheckM2 only (high >90/<5, medium >=50/<10); MIMAG_final is the workspace_2 call incl. rRNA/tRNA for the 29 final MAGs.\n'
    'Raw bins overlap: the same genome can appear as a MetaBAT2, VAMB, SemiBin2 and DAS_Tool bin; counts are bins, not genomes.\n'
    'Panel A y-axis symlog (linear below 10%). Panel B: 14 most frequent phyla, rest folded into Other; target phyla in orange; x-scales differ between libraries.\n'
    'not classifiable = GTDB-Tk excluded the bin (<10% of MSA columns): 121 of 122 such bins are <30% complete; the exception, LIQ_vamb_029, is a 19.8 Mb / 103% contamination VAMB catch-all bin.\n')
print(open(f'{O}/catalogue_notes.txt').read())
