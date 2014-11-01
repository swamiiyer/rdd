""" 
This script generates data for estimating the algorithmic complexity of the 
find_one_S_m(), find_one_S_mu(), and find_one_S_m_mu_hybrid() functions 
in rdd.py. The data consist of r (for replicates) files, each consisting of 
n permutations of order m. The permutations generated are such that the   
mean and standard deviation of their distances from identity are within a 
specified tolerance level of the analytic values. The data files are named 
<m>_<n>_<r>_<t>.txt, where m, n, r, and t denote the permutation order, number 
of permutations per replicate, replicate number, and tolerance respectively.
"""

import argparse, numpy, random, rdd, sys

def main(args):
    """
    Entry point.
    """

    parser = argparse.ArgumentParser(description = """Generates data for 
             estimating the algorithmic complexity of the find_one_S_m(), 
             find_one_S_mu(), and find_one_S_m_mu_hybrid() functions in 
             rdd.py""")
    parser.add_argument("-m", dest = "order", type = int, metavar = "...", 
                        required = True, 
                        help = """permutation order""")
    parser.add_argument("-r", dest = "replicates", type = int, 
                        required = False, metavar = "...", default = 1, 
                        help = """number of replicates to generate""")
    parser.add_argument("-n", dest = "samples", type = int, 
                        required = True, metavar = "...", 
                        help = """number of permutations per replicate""")
    parser.add_argument("-t", dest = "tolerance", type = float, 
                        required = False, metavar = "...", default = 0.01, 
                        help = """tolerance for randomness of samples; 
                        default = 0.01""")
    args = parser.parse_args()
    expected_mean = (args.order ** 2 - 1) / 3.0
    expected_std = ((args.order + 1) * (2 * args.order ** 2 + 7) / 45.0) ** 0.5
    replicate = 1
    v = rdd.id(args.order)
    while replicate <= args.replicates:
        print "Generating order %d, replicate %d..." %(args.order, replicate)
        samples = []
        u = range(1, args.order + 1)
        error_mean, error_std = float('inf'), float('inf')
        count = 1
        while count <= args.samples:
            random.shuffle(u)
            if u not in samples:
                samples.append(list(u))
                count += 1
        distances = [rdd.dist(u, v) for u in samples]
        mean = numpy.mean(distances)
        std = numpy.std(distances)
        error_mean = abs(expected_mean - mean) / expected_mean
        error_std = abs(expected_std - std) / expected_std
        if error_mean < args.tolerance and error_std < args.tolerance:
            fh = open("%d_%d_%d_%s.txt" %(args.order, args.samples, 
                                       replicate, args.tolerance), "w")
            for u in samples:
                fh.write("%s\n" %(" ".join(map(str, u))))
            fh.close()
            replicate += 1

if __name__ == "__main__":
    main(sys.argv)
