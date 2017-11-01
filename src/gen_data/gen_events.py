import dendropy
import random
import itertools
import sys

path = "../../data/yule_trees/sp10_trees.nwk"
tree_list = dendropy.TreeList.get(
        path=path,
        schema="newick")

prop_D = 0.33
prop_L = 0.33
prop_S = abs(prop_D - prop_L)

freq = [0,0,0]
freq_rel = [0,0,0]

def event_newick(node):
    node.trasfer = False
    if len(node.child_nodes()) == 0:
        print(node.taxon.label,end='')
        return
    print("(", end='')
    for child in node.child_nodes()[:-1]:
        event_newick(child)
        print(",", end='')
    event_newick(node.child_nodes()[-1])
    print(")", end='')
    p = random.uniform(0,1)
    event = "S"
    if p <= prop_D:
        event = "D"
        freq[1] += 1
    elif p <= prop_D + prop_L:
        event = "L"
        freq[2] += 1
        transfer_dir = random.randrange(0,1)
        node.child_nodes()[transfer_dir].transfer = True
    else:
        freq[0] += 1
    node.event = event
    print(event, end='')

def is_transfer(lca_node, leaf):
    while leaf.parent_node != lca_node:
        leaf = leaf.parent_node
    return leaf.transfer

def ortho_matrix(tree, leaves):
    n = len(leaves)
    ortho = [[0 for _ in range(n)] for _ in range(n)]
    for i, j in itertools.combinations(range(n), 2):
        n1, n2 = leaves[i], leaves[j]
        assert(n1.taxon < n2.taxon)
        lca = tree.mrca(taxon_labels=[n1.taxon.label, n2.taxon.label])
        if lca.event == "S":
            freq_rel[0] += 1
            ortho[i][j] = 1
            ortho[j][i] = 1
        elif lca.event == "D":
            freq_rel[1] += 1
            ortho[i][j] = 0
            ortho[j][i] = 0
        elif is_transfer(lca, n2):
            freq_rel[2] += 1
            ortho[i][j] = 1
            ortho[j][i] = 0
        else:
            freq_rel[2] += 1
            ortho[i][j] = 0
            ortho[j][i] = 1
    return ortho



for i, tree in enumerate(tree_list):
    # print(tree.as_string(schema="newick"))
    print("#GeneFamily", i)

    leaves = [l for l in tree.leaf_node_iter()]
    leaves.sort(key = lambda x: x.taxon)
    print("#GeneNames")
    for leaf in leaves:
        print(leaf.taxon.label, end=' ')
    print()

    print("#LabeledGeneTree=", end='')
    for node in tree.nodes():
        node.transfer = False
    event_newick(tree.seed_node)
    print(";")

    print("#OrthologyMatrix")
    om = ortho_matrix(tree, leaves)
    for row in om:
        for v in row:
            print(v, end=' ')
        print()
    print()

    # print(leaves)


    # break

total = sum(freq)
sys.stderr.write("Distribution: " + str([v/total for v in freq]) + "\n")
total = sum(freq_rel)
sys.stderr.write("Distribution: " + str([v/total for v in freq_rel]) + "\n")

