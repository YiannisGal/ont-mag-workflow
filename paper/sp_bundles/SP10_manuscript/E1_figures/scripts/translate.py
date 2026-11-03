#!/usr/bin/env python3
"""Switch the copied thesis figure scripts to English labels and manuscript names.

Only strings that are DRAWN are translated; caption strings go to S.caption(), which no
longer draws anything in the manuscript style module. Every replacement must match
exactly once or more, otherwise the script stops. Numbers and data handling are untouched.
"""
import re, sys, os
D = os.path.dirname(os.path.abspath(__file__))

COMMON = [('FIG_W_MM = 160.0', 'FIG_W_MM = 170.0'), (' κ.ά.', ' et al.'), ('"κ.ά. "', '"et al. "'),
          ('"υποκ./θέση"', '"subst./site"'), ('"συμπτυγμένη ομάδα αναφοράς"', '"collapsed reference clade"'),
          ('"αλληλουχίες αυτής της εργασίας"', '"sequences from this study"')]

MAP = {
 "fig01_mag_quality.py": [
  ('"Ποιότητα των 29 MAGs κατά τα κριτήρια MIMAG"', '""'),
  ('"Πληρότητα (%)"', '"Completeness (%)"'), ('"Επιμόλυνση (%)"', '"Contamination (%)"'),
  ('"≥ 90 % πληρότητα, < 5 % επιμόλυνση"', '"≥ 90% complete, < 5% contaminated"'),
  ('"μέτρια ποιότητα"', '"medium quality"'), ('"μέγεθος γονιδιώματος"', '"genome size"'),
  ('label="MIMAG υψηλής ποιότητας"', 'label="MIMAG high-quality"'),
  ('S.save(fig, "fig01_mag_quality"', 'S.save(fig, "S_MAG_quality"')],
 "fig03_function_heatmap.py": [
  ('"μεθυλαμίνη"', '"methylamine"'), ('"μεθανόλη"', '"methanol"'), ('"υδρογονάση"', '"hydrogenase"'),
  ('"B$_{12}$ σύνθεση"', '"B$_{12}$ synthesis"'), ('"B$_{12}$ εξάρτηση"', '"B$_{12}$-dependent"'),
  ('"B$_{12}$ μεταφορά"', '"B$_{12}$ transport"'),
  ('"Λειτουργικό δυναμικό των 29 MAGs"', '""'),
  ('"αριθμός γονιδίων, ως προς το μέγιστο της στήλης"', '"gene count, relative to column maximum"'),
  ('S.save(fig, "fig03_function_heatmap"', 'S.save(fig, "S_function_heatmap"')],
 "fig06_exchange_matrix.py": [
  ('"Προβλεπόμενη ανταλλαγή μεταβολιτών μεταξύ των γονιδιωμάτων"', '""'),
  ('"δέκτης"', '"receiver"'), ('"δότης"', '"donor"'), ('"άθροισμα SMETANA ανά ζεύγος"', '"summed SMETANA score per pair"'),
  ('"σύνολο\\nπροσφοράς"', '"total\\ndonated"'), ('"σύνολο λήψης"', '"total received"'),
  ('"καμία προβλεπόμενη μεταφορά"', '"no predicted transfer"'), ('"δεν ορίζεται"', '"not defined"'),
  ('"γονιδίωμα < 90 % πλήρες (°)"', '"genome < 90% complete (°)"'),
  ('S.save(fig, "fig06_exchange_matrix"', 'S.save(fig, "S_exchange_matrix"')],
 "fig07_what_is_exchanged.py": [
  ('"Τι ανταλλάσσεται: αμινοξέα, όχι υδρογόνο"', '""'),
  ('"άθροισμα βαθμών SMETANA σε όλα τα ζεύγη"', '"summed SMETANA score over all pairs"'),
  ('"αμινοξέα"', '"amino acids"'), ('"ενεργειακά και C1"', '"energy and C1"'),
  ('"άλλα (νουκλεοσίδια, οργανικά οξέα, βιταμίνες)"', '"other (nucleosides, organic acids, vitamins)"'),
  ('"τα υποστρώματα που συζητά το κεφάλαιο"', '"substrates discussed in the text"'),
  ('S.save(fig, "fig07_what_is_exchanged"', 'S.save(fig, "S_compound_classes"')],
 "fig09_biogeography.py": [
  ('"Βιογεωγραφία: τα 29 γονιδιώματα δεν ανακτώνται αλλού"', '""'),
  ('"εύρος κάλυψης (κλάσμα του γονιδιώματος)"', '"breadth of coverage (fraction of genome)"'),
  ('"αναγνώσεις Αιτωλικού"', '"Aitoliko reads"'), ('"καλύτερο δημόσιο σύνολο δεδομένων"', '"best public dataset"'),
  ('"όριο ελέγχου %s"', '"detection threshold %s"'),
  ('S.save(fig, "fig09_biogeography"', 'S.save(fig, "S_biogeography"')],
 "fig10_dsr_binning.py": [
  ('"Σε ποιο γονιδίωμα ταιριάζουν τα contigs του dsrAB"', '""'),
  ('"τετρανουκλεοτιδική απόσταση από το γονιδίωμα"', '"tetranucleotide distance to genome"'),
  ('"μέτρο του ίδιου του γονιδιώματος"', '"genome\'s own yardstick"'),
  ('"εκτός — ασύμβατο"', '"outside: incompatible"'), ('"εντός — συμβατό"', '"inside: compatible"'),
  ('"χωρίς μέτρο σύγκρισης"', '"no yardstick"'), ('"οι υπόλοιπες %d: %s – %s"', '"remaining %d: %s – %s"'),
  ('S.save(fig, "fig10_dsr_binning"', 'S.save(fig, "S_dsr_binning"')],
 "fig11_ani_isolates.py": [
  ('"Ταυτότητα γονιδιωμάτων: MAG, απομονωμένα στελέχη, γονιδιώματα αναφοράς"', '""'),
  ('"MAG × στέλεχος"', '"MAG × isolate"'), ('"MAG × αναφορά"', '"MAG × reference"'),
  ('"στέλεχος × στέλεχος"', '"isolate × isolate"'), ('"στέλεχος × αναφορά"', '"isolate × reference"'),
  ('"αναφορά × αναφορά"', '"reference × reference"'), ('"συγκρίθηκαν"', '"compared"'), ('"με τιμή"', '"with value"'),
  ('"μέγ. ANI"', '"max. ANI"'), ('"καμία"', '"none"'), ('"όριο είδους, 95 %"', '"species threshold, 95%"'),
  ('"θετικός μάρτυρας"', '"positive control"'), ('"στέλεχος – στέλεχος"', '"isolate – isolate"'),
  ('"στέλεχος – αναφορά (ταυτοποίηση)"', '"isolate – reference (identification)"'),
  ('"αναφορά – αναφορά"', '"reference – reference"'), ('"MAG – καλλιεργημένου: καμία τιμή"', '"MAG – cultured: no value"'),
  ('S.save(fig, "fig11_ani_isolates"', 'S.save(fig, "S_ani_isolates"')],
 "fig12_phylogenomics.py": [
  ('"Φυλογενωμική θέση των αρχαίων γονιδιωμάτων του καταλόγου"', '""'),
  ('"Φυλογενωμική θέση των βακτηριακών γονιδιωμάτων του καταλόγου"', '""'),
  ('"γονιδιώματα αυτής της εργασίας"', '"genomes from this study"'), ('"γονιδιώματα αναφοράς GTDB"', '"GTDB reference genomes"'),
  ('"εξωομάδα"', '"outgroup"'), ('"(χωρίς γένος στη GTDB)"', '"(no GTDB genus)"')],
 "fig14_15_mcra.py": [
  ('"Ποια αρχαία γονιδιώματα φέρουν mcrA πάνω από το επιμελημένο όριο"', '""'),
  ('"Η τοποθέτηση της mcrA ως προς την ομάδα ελέγχου των ACR"', '""'),
  ('"Κάθε εύρημα K00399 κάθε αρχαίου γονιδιώματος του καταλόγου"', '"Every K00399 hit in every archaeal MAG"'),
  ('"βαθμολογία KOfam (πλήρους αλληλουχίας)"', '"KOfam score (full sequence)"'),
  ('"επιμελημένο όριο K00399"', '"curated K00399 threshold"'),
  ('"αναγωγάση της μεθυλο-συνενζύμου Μ (MCR)"', '"methyl-coenzyme M reductase (MCR)"'),
  ('"αναγωγάση αλκυλο-συνενζύμου Μ (ACR) — ομάδα ελέγχου"', '"alkyl-coenzyme M reductase (ACR), control"'),
  ('"Δέντρο μέγιστης πιθανοφάνειας, ριζωμένο στις αλληλουχίες ACR"', '"ML tree rooted on ACR"'),
  ('"ομάδα ελέγχου"', '"control"'), ('"άλλες πρωτεΐνες, όχι έλεγχος υπομονάδας α"', '"other proteins, not alpha-subunit control"'),
  ('"ένζυμο περιορισμού· εκτός σύγκρισης"', '"restriction enzyme; excluded"'),
  ('"ACR, υπομονάδα α — η ομάδα ελέγχου"', '"ACR alpha subunit, control"'),
  ('"υπομονάδες β/γ και «component A2» του ίδιου συνόλου"', '"beta/gamma subunits and component A2"'),
  ('"εσφαλμένα εισαγμένες αλληλουχίες (McrA ένζυμο περιορισμού)"', '"mis-annotated (McrA restriction enzyme)"'),
  ('"στο δέντρο"', '"in tree"'), ('"δεν προτάθηκε"', '"not nominated"'), ('"κάλυψη ιζήματος"', '"sediment breadth"'),
  ('"SH-aLRT (%) — το αρχείο φέρει δύο τιμές ανά κόμβο"', '"SH-aLRT (%)"'),
  ('"ACRA|ACR α · %s (n = %d)"', '"ACRA|ACR alpha · %s (n = %d)"'),
  ('S.save(fig, "fig14_mcra_nomination"', 'S.save(fig, "S_mcrA_nomination"'),
  ('S.save(fig, "fig15_mcra_tree"', 'S.save(fig, "fig15_thesis_tree_unused"')],
}

def strip_comments_docstrings(src):
    out, in_doc = [], False
    for ln in src.splitlines():
        s = ln.strip()
        if s.startswith('"""') or s.endswith('"""'):
            if s.count('"""') == 1: in_doc = not in_doc
            continue
        if in_doc or s.startswith("#"): continue
        out.append(ln)
    return out

GREEK = re.compile(r"[Ͱ-Ͽἀ-῿]")
for fn, pairs in MAP.items():
    p = os.path.join(D, fn); src = open(p).read()
    for a, b in pairs + COMMON:
        if a in src:
            src = src.replace(a, b)
        elif (a, b) not in COMMON:
            sys.exit("%s: not found: %s" % (fn, a))
    if fn == "fig07_what_is_exchanged.py":   # compound names: English = the key itself
        src = re.sub(r'"([^"]+)": "[^"]*[Ͱ-Ͽ][^"]*"', lambda m: '"%s": "%s"' % (m.group(1), "H₂" if m.group(1) == "H2" else m.group(1)), src)
    open(p, "w").write(src)
    left = [l for l in strip_comments_docstrings(src) if GREEK.search(l)]
    print("%-28s greek code lines left: %d" % (fn, len(left)))
