"""Is the oxidative (rDsrAB) label set monophyletic in the tree? Reports the smallest clade containing all Oxidative leaves."""
import sys, os
from Bio import Phylo
g = sys.argv[1]; H = os.path.dirname(os.path.abspath(__file__))
tax = {l.split('\t')[0]: l.split('\t')[1].split(';')[0] for l in open(f'{H}/{g}/taxon.tsv')}
t = Phylo.read(f'{H}/{g}/{g}.treefile', 'newick')
red = [c for c in t.get_terminals() if tax.get(c.name) == 'Reductive']
t.root_with_outgroup(red[0])   # any reductive leaf as temporary root
ox = [c for c in t.get_terminals() if tax.get(c.name) == 'Oxidative']
m = t.common_ancestor(ox); inside = [c.name for c in m.get_terminals()]
nred = sum(tax.get(n) == 'Reductive' for n in inside)
print(f'{g}: {len(ox)} oxidative leaves; smallest clade containing them has {len(inside)} leaves, of which {nred} labelled reductive')
