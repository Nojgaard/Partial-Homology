import subprocess as sp

def triple_diff(line):
    trees = line.split(";")
    with open("org_tree.trees", "w") as f:
        f.write(trees[0])
    with open("constructed_tree.trees", "w") as f:
        f.write(trees[1])
    # print(trees)
    sp.run(['java', '-jar','../../lib/TreeCmp/bin/TreeCmp.jar', '-r', 
            'org_tree.trees', '-d', 'tt', '-N', '-i', 'constructed_tree.trees',
            '-o', 'output.dist'])
    print(trees)
    with open("output.dist", "r") as f:
        for line in f.readlines():
            line = line.strip().split()
            print(line)
            if line[0] == "1":
                if line[3] != 'N/A': return line[3] 
    return 0

with open("../../data/trees/dl50_hl50.trees", "r") as f, open("out.dat", "w") as fo:
    u,dl,hl = 0, 0, 0
    i, sum_diff = 0, 0
    fo.write("prop\tdl\thl\ttriple_diff\n")
    for line in f.readlines():
        stats = line.strip().split()
        if stats[0] == "GROUP":
            if u != 0:
                fo.write("%s\t%s\t%s\t%s\n"%(u,dl,hl,str(sum_diff/i)))
                i, sum_diff = 0, 0
            u,dl,hl = stats[1:]

        else:
            sum_diff += float(triple_diff(stats[0]))
            i += 1
            # fo.write("%s\t%s\t%s\t%s\n"%(u,dl,hl,diff))
