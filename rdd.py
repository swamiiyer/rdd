""" 
This script computes the segment center(s) and distinguishability of two 
permutations u and v.
"""

import argparse, rdd, sys

def is_valid(u):
    """
    Returns True if u is a valid permutation and False otherwise.
    """

    for i in range(1, len(u) + 1):
        if not i in u:
            return False
    return True

def id(m):
    """
    Returns the identity permutation of length m.
    """

    return range(1, m + 1)

def is_id(u):
    """
    Returns True if u is the identity permutation and False otherwise.
    """

    return u == id(len(u))

def compose(u, v):
    """
    Returns the composition of permutations u and v.
    """

    return [u[v[i - 1] - 1] for i in range(1, len(u) + 1)]

def inv(u):
    """
    Returns the inverse permutation of u.
    """
    
    uinv = id(len(u))
    for i, ui in enumerate(u):
        uinv[ui - 1] = i + 1
    return uinv

def tau(m, i, j):
    """
    Returns the identity permutation of length m with values i and j 
    transposed.
    """

    u = id(m)
    temp = u[i - 1]
    u[i - 1] = u[j - 1]
    u[j - 1] = temp
    return u

def linversion(m, u, i, j):
    """
    Returns a permutation of u with the values i and j transposed.
    """

    return rinversion(m, u, u.index(i) + 1, u.index(j) + 1)

def rinversion(m, u, i, j):
    """
    Returns a permutation of u with the values at i and j transposed.
    """

    w = list(u)
    temp = w[i - 1]
    w[i - 1] = w[j - 1]
    w[j - 1] = temp
    return w

def dist(u, v):
    """
    Returns the Spearman Footrule distance between permutations u and v.
    """

    return sum([abs(u[i - 1] - v[i - 1]) for i in range(1, len(u) + 1)])

def radius(u, v, w):
    """
    Returns the smallest radius of a ball centered at w that covers 
    permutations u and v.
    """

    return max(dist(u, w), dist(v, w))

def ideal_radius(u, v):
    """
    Returns the ideal radius for a ball that covers permutations u and v.
    """

    d = dist(u, v) / 2
    return d if d % 2 == 0 else d + 1

def composition(m, u):
    """
    Returns the canonical composition mu of the permutation u. For example, 
    mu([1, 4, 2, 3, 5]) = [2, 3], mu([1, 2, 3, 4, 5]) = [5], and 
    mu([5, 4, 3, 2, 1]) = [1, 1, 1, 1, 1].
    """
    
    mu = []
    descent = 0
    for i in range(1, m):
        if u[i] < u[i - 1]:
            mu.append(i - descent)
            descent = i
    mu.append(i + 1 - descent)
    return mu
    
def is_mu_sorted(m, mu, u):
    """
    Returns True if the permutation u is mu-sorted and False otherwise.
    """

    descent, j = 0, 0
    for i in range(1, m):
        if u[i] < u[i - 1] and i != mu[j] + descent:
            return False
        elif i == mu[j] + descent:
            descent = i
            j += 1
    return True

def mu_sort(mu, u):
    """
    Returns the mu-sorted permutation of u.
    """
    
    l, i = [], 0
    for n in mu:
        l += sorted(u[i:i + n])
        i += n
    return l

def eligible_neighbor_rec(m, mu, u, v, i, j, is_mu_sorted, visited):
    """
    Returns a neighboring (i.e. differing by a single transposition) 
    permutation w of u that is mu-sorted and not in visited such that 
    w is in the segment [u, v]; candidate permutations are enumerated 
    recursively, starting with the transposition immediately following (i, j). 
    Returns None if no such permutation exists. 
    """

    if j < m:
        j += 1
    elif i < m - 1:
        i += 1
        j = i + 1
    else:
        return None, None, None
    if v[i - 1] <= u[j - 1] and u[j - 1] < u[i - 1] and \
            u[i - 1] <= v[j - 1] or v[j - 1] <= u[i - 1] and \
            u[i - 1] < u[j - 1] and u[j - 1] <= v[i - 1]:
        w = rinversion(m, u, i, j)
        if is_mu_sorted(m, mu, w) and w not in visited:
            return w, i, j
    return eligible_neighbor_rec(m, mu, u, v, i, j, is_mu_sorted, visited)

def eligible_neighbor_iter(m, mu, u, v, i, j, is_mu_sorted, visited):
    """
    Returns a neighboring (i.e. differing by a single transposition) 
    permutation w of u that is mu-sorted and not in visited such that 
    w is in the segment [u, v]; candidate permutations are enumerated 
    iteratively, starting with the transposition immediately following (i, j). 
    Returns None if no such permutation exists. 
    """

    w = None
    while True:
        if j < m:
            j += 1
        elif i < m - 1:
            i += 1
            j = i + 1
        else:
            w, i, j = None, None, None
            break
        if v[i - 1] <= u[j - 1] and u[j - 1] < u[i - 1] and \
           u[i - 1] <= v[j - 1] or v[j - 1] <= u[i - 1] and \
           u[i - 1] < u[j - 1] and u[j - 1] <= v[i - 1]:
            w = rinversion(m, u, i, j)
            if is_mu_sorted(m, mu, w) and w not in visited:
                break
    return w, i, j

def find_one(m, mu, u, v, translate, is_mu_sorted, mu_sort):
    """
    Returns the number of edges visited, distinguishability and a center 
    of the segment [u, v]. If translate is True, then the center of the 
    segment [uv', id] is computed and translated to the center of [u, v] 
    via right composition by v. Otherwise, the center of [u, v] is computed 
    directly. The function is_mu_sorted specifies whether a permutation is 
    mu-sorted or not and the function mu_sort mu-sorts a permutation; these 
    functions are appropriately defined by the caller of find_one in order to 
    search for the center in S_m, S_mu, or S_m_mu_hybrid.
    """

    if translate:
        q, u, v = v, compose(u, inv(v)), id(m)
    duv = dist(u, v)
    ideal, d1, radius, center = ideal_radius(u, v), 0, duv, u
    visited, dfs = [], []
    i, j = 1, 1
    while True:
        w, i, j = eligible_neighbor_iter(m, mu, u, v, i, j, 
                                         is_mu_sorted, visited)
        if w == None:
            if len(dfs) == 0:
                break
            else:
                u, i, j, d1 = dfs.pop()
        else:
            visited.append(w)
            d2 = d1 + 2 * abs(u[i - 1] - w[i - 1])
            if d2 == ideal:
                radius = ideal
                center = w
                break
            elif d2 > ideal:
                if d2 < radius:
                    radius = d2
                    center = w
            else:
                if duv - d2 < radius:
                    radius = duv - d2
                    center = w
                if radius == ideal:
                    break
                dfs.append((u, i, j, d1))
                u, i, j, d1 = w, 1, 1, d2
    return len(visited), radius, mu_sort(mu, compose(center, q) \
                                         if translate else center)

def find_one_S_m(m, mu, u, v):
    """
    Returns the number of edges visited, distinguishability and a center 
    of the segment [u, v] by searching for it in S_m.
    """

    return find_one(m, mu, u, v, True, lambda m, mu, u: True, lambda m, u: u)

def find_one_S_mu(m, mu, u, v):
    """
    Returns the number of edges visited, distinguishability and a center 
    of the segment [u, v] by searching for it in S_mu.
    """

    return find_one(m, mu, u, v, False, is_mu_sorted, lambda m, u: u)

def find_one_S_m_mu_hybrid(m, mu, u, v):
    """
    Returns the number of edges visited, distinguishability and a center 
    of the segment [u, v] by searching for it in S_m and then mu-sorting 
    the result.
    """

    return find_one(m, mu, u, v, True, lambda m, mu, u: True, mu_sort)

def find_all(m, mu, u, v, radius, translate, is_mu_sorted, mu_sort):
    """
    Returns all the centers of the segment [u, v] having the 
    distinguishability specified by the radius parameter. If translate is 
    True, then the centers of the segment [uv', id] is computed and 
    translated to the centers of [u, v] via right composition   
    by v. Otherwise, the centers of [u, v] are computed directly. The 
    function is_mu_sorted specifies whether a permutation is mu-sorted or 
    not and the function mu_sort mu-sorts a permutation; these functions 
    are appropriately defined by the caller of find_all in order to 
    search for the centers in S_m, S_mu, or S_m_mu_hybrid.
    """

    if translate:
        q, u, v = v, compose(u, inv(v)), id(m)
    centers = []
    duv = dist(u, v)
    if duv == radius:
        centers.append(mu_sort(mu, compose(u, q) if translate else u))
    d1, center = 0, u
    visited, dfs = [], []
    i, j = 1, 1
    while True:
        w, i, j = eligible_neighbor_iter(m, mu, u, v, i, j, 
                                         is_mu_sorted, visited)
        if w == None:
            if len(dfs) == 0:
                break
            else:
                u, i, j, d1 = dfs.pop()
        else:
            visited.append(w)
            d2 = d1 + 2 * abs(u[i - 1] - w[i - 1])
            if d2 == radius:
                centers.append(mu_sort(mu, compose(w, q) if translate else w))
            elif d2 < radius:
                dfs.append((u, i, j, d1))
                u, i, j, d1 = w, 1, 1, d2
    return centers

def find_all_S_m(m, mu, u, v):
    """
    Returns all the centers of the segment [u, v] as a generator by searching 
    for them in S_m. 
    """

    edges, radius, center = find_one_S_m(m, mu, u, v)
    return find_all(m, mu, u, v, radius, True, lambda m, mu, u: True, 
                    lambda m, u: u)

def find_all_S_mu(m, mu, u, v):
    """
    Returns all the centers of the segment [u, v] as a generator by searching 
    for them in S_mu. 
    """

    edges, radius, center = find_one_S_m(m, mu, u, v)
    return find_all(m, mu, u, v, radius, False, is_mu_sorted, lambda m, u: u)

def find_all_S_m_mu_hybrid(m, mu, u, v):
    """
    Returns all the centers of the segment [u, v] as a generator by searching 
    for them in S_mu. 
    """

    edges, radius, center = find_one_S_m(m, mu, u, v)
    return find_all(m, mu, u, v, radius, False, is_mu_sorted, lambda m, u: u)

def main(args):
    """
    Entry point.
    """

    # Parse and validate command-line arguments.
    parser = argparse.ArgumentParser(description = """Computes the segment 
                 center(s) and distinguishability of two permutations 
                 u and v""")
    parser.add_argument("-c", dest = "count", type = str, 
                        required = False, default = "one", 
                        choices = ["one", "all"],
                        help = """one or all centers; default = one""")
    parser.add_argument("-a", dest = "algorithm", type = str, 
                        required = False, default = "S_m", 
                        choices = ["S_m", "S_mu", "S_m_mu_hybrid"],
                        help = """search algorithm; default = S_m""")
    parser.add_argument("-u", dest = "u", type = str, metavar = "...", 
                        required = True, 
                        help = """permutation u; e.g., \"2 3 4 5 1\"""")
    parser.add_argument("-v", dest = "v", type = str, metavar = "...", 
                        required = True, 
                        help = """permutation v; e.g., \"3 5 1 4 2\"""")
    parser.add_argument("-mu", dest = "mu", type = str, metavar = "...", 
                        required = False, 
                        help = """composition; e.g., \"2 2 1\"; 
                        default = \"1 1 ... 1\"""")
    args = parser.parse_args()
    u = map(int, args.u.split())
    v = map(int, args.v.split())
    if not is_valid(u):
        print "ERROR: u is not a valid permutation"
        sys.exit(1)
    if not is_valid(v):
        print "ERROR: v is not a valid permutation"
        sys.exit(1)
    if len(u) != len(v):
        print "ERROR: u and v have different lengths"
        sys.exit(1)
    m = len(u)
    if args.mu == None:
        mu = [1 for i in range(m)]
    else:
        mu = map(int, args.mu.split())
    if sum(mu) != m:
        print "ERROR: mu is not a valid composition"
        sys.exit(1)
    elif not is_mu_sorted(m, mu, u):
        print "ERROR: u is not mu-sorted"
        sys.exit(1)
    elif not is_mu_sorted(m, mu, v):
        print "ERROR: v is not mu-sorted"
        sys.exit(1)

    # Find and print the segment center(s) and distinguishability of 
    # permutations u and v.
    centers = getattr(rdd, "find_%s_%s" \
                          %(args.count, args.algorithm))(m, mu, u, v)
    print "u: %s" %(u)
    print "v: %s" %(v)
    print "Ideal radius: %s" %(ideal_radius(u, v))
    if args.count == "one":
        print "Radius: %s" %(centers[1])
        print "Center: %s" %(centers[2])
    else:
        first = centers[0]
        print "Radius: %s" %(radius(u, v, first))
        print "Centers:"
        print "  %s" %(first)
        for center in centers[1:]:
            print "  %s" %(center)

if __name__ == "__main__":
    main(sys.argv)
