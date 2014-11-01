"""
This script generates a search graph in which the vertices are permutations 
of order m and two permutations u and v are connected by an edge if 
there are indices 1 <= i, j <= m such that u[i] = u[j] + 1 and 
v[i] + 1 = v[j].
"""

import pickle, sys

def search_graph(m):
    """
    Return the search graph of permutations of order m.
    """

    u = range(1, m + 1)
    u.reverse()
    v = range(1, m + 1)
    q = [u]
    graph = {}
    while len(q) > 0:
        w = q.pop(0)
        children = []
        for i in range(m - 1):
            for j in range(i + 1, m):
                if w[i] == w[j] + 1:
                    child = list(w)
                    child[i], child[j] = child[j], child[i]
                    children.append(child)
                    if not child in q:
                        q.append(child)
        graph[str(w)] = children
    return graph

def main(args):
    """
    Entry point.
    """

    if len(args) != 2:
        print "Usage: python search_graph.py <order> <outfile>"
        sys.exit(1)
    m = int(args[0])
    outfile = args[1]
    pickle.dump(search_graph(m), open(outfile, "wb"))

if __name__ == "__main__":
    main(sys.argv[1:])
