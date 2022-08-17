from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
from sys import argv

id = argv[1]
species = argv[2]
marks = ['1']

print(f"../evoltree/codon_aln/sp_name_aln_{id}.paml")
tree = EvolTree(f"../evoltree/species_tree/SpeciesTree.txt")
tree.link_to_alignment(f"../evoltree/codon_aln/sp_name_aln_{id}.paml")

if species.lower() == 'a99':
    marks = ['1']
elif species.lower() == 'mcon':
    marks = ['5']
elif species.lower() == 'cvar':
    marks = ['3']

# mark a group of branches
tree.mark_tree (marks, ['#1'])
print(tree.write ())

tree.run_model('b_free')
tree.run_model('M0')
pval = tree.get_most_likely ('b_free','M0')
if pval < 0.05:
    print(pval, 'free-branch model most likely')
    tree.render(f'../Results/{id}.svg', layout=evol_clean_layout)
else:
    print(pval, 'dN/dS not significantly different between branches')