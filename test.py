from Bio import Phylo
import matplotlib.pyplot as plt

# Load the phylogenetic tree
tree = Phylo.read("results/phylogenetic_tree.xml", "phyloxml")

# Visualize the tree
Phylo.draw(tree)
plt.show()