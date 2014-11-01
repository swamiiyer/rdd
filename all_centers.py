"""
This script takes the search graph of order m produced by search_graph.py and 
for each each permutation u, finds all (segment and non-segment) 
centers of (u, id).
"""

import itertools, pickle, rdd, sys

def find_centers(graph, u, v):
    """
    Searches for all the segment and non-segment centers of (u, v) in the 
    specified search graph and returns them as two separate lists.
    """

    seg_centers, non_seg_centers = [], []
    q = [u]
    duv = rdd.dist(u, v)
    min_radius = duv
    while len(q) > 0:
        w = q.pop(0)
        duw, dwv = rdd.dist(u, w), rdd.dist(w, v)
        radius = rdd.radius(u, v, w)
        if radius == min_radius:
            if duv == duw + dwv:
                seg_centers.append(w)
            else:
                non_seg_centers.append(w)
        elif radius < min_radius:
            min_radius = radius
            if duv == duw + dwv:
                seg_centers = [w]
                non_seg_centers = []
            else:
                seg_centers = []
                non_seg_centers = [w]
        children = graph[str(w)]
        for child in children:
            if not child in q:
                q.append(child)
    return seg_centers, non_seg_centers

def main(args):
    """
    Entry point.
    """

    if len(args) != 2:
        print "Usage: python all_centers.py <order> <infile>"
        sys.exit(1)
    m = int(args[0])
    infile = args[1]
    graph = pickle.load(open(infile, "rb"))
    v = rdd.id(m)
    permutations = itertools.permutations(v)
    for permutation in permutations:
        u = list(permutation)
        seg_centers, non_seg_centers = find_centers(graph, u, v)
        seg_radius = rdd.radius(u, v, seg_centers[0])
        non_seg_radius = rdd.radius(u, v, non_seg_centers[0]) \
                         if len(non_seg_centers) > 0 else "-"
        print "u = %s" %(u)
        print "  ideal radius = %d" %(rdd.ideal_radius(u, v))
        print "  seg radius = %d" %(seg_radius)
        print "  non-seg radius = %s" %(non_seg_radius)
        print "  segment centers:" 
        for center in seg_centers:
            print "    %s" %(center)
        print "  non-segment centers:"
        for center in non_seg_centers:
            print "    %s" %(center)
        print ""
    
if __name__ == "__main__":
    main(sys.argv[1:])
