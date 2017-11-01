import dendropy
from dendropy.model import birthdeath

num_leaves = 25
niter = 1000

for _ in range(niter):
    taxa = dendropy.TaxonNamespace(["l"+str(i) for i in range(num_leaves)])
    tree = birthdeath.uniform_pure_birth_tree(taxa, birth_rate=1.0)
    print(tree.as_string(schema="newick"))
    # print(tree.as_ascii_plot())

