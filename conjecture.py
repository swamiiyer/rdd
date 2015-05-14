"""
This script validates the following conjecture for small permutation orders:

Every connected component of a graph of all the centers of permutations 
u and v contains at least one vertex on the segment [u, v].
"""

import all_centers, itertools, networkx, pickle, rdd, sys

def main(args):
    if len(args) != 2:
        print "Usage: python conjecture.py <order> <infile>"
        sys.exit(1)
    m = int(args[0])
    infile = args[1]
    graph1 = pickle.load(open(infile, "rb"))
    permutations = itertools.permutations(range(1, m + 1))
    graph2 = networkx.Graph()
    for i, permutation in enumerate(permutations):
        permutation = list(permutation)
        for j in range(1, m):
            for k in range(j + 1, m + 1):
                neighbor = rdd.linversion(m, permutation, j, k)
                graph2.add_edge(str(permutation), str(neighbor))
    permutations = itertools.permutations(range(1, m + 1))
    v = rdd.id(m)
    for u in permutations:
        u = list(u)
        seg_centers, non_seg_centers = all_centers.find_centers(graph1, u, v)
        centers = seg_centers + non_seg_centers
        graph3 = networkx.subgraph(graph2, [str(i) for i in centers])
        print "u = %s" %(u)
        for i, cc in enumerate(networkx.connected_components(graph3)):
            print "  component %d" %(i + 1)
            count, seg_count = 0, 0
            for center in centers:
                if str(center) in cc:
                    count += 1
                    if center in seg_centers:
                        seg_count += 1
                        print "    %s seg" %(center)
                    else:
                        print "    %s" %(center)
            print "  total = %d, on seg = %d" %(count, seg_count)
        print ""

if __name__ == "__main__":
    main(sys.argv[1:])
