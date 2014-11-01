""" 
This script generates plots (as PDFs) for visualizing the algorithmic 
complexity of the find_one_S_m (algorithm A), find_one_S_mu (algorithm B), 
and find_one_S_m_mu_hybrid (algorithm C) functions in rdd.py. These 
plots include: average search time, average number of edges traversed, and 
excess of distinguishability over the ideal radius; the x-axis in each of 
the plots is the permutation order. The input to the script are the 
names of the result files (one per replicate) produced by complexity.py, 
specified via STDIN. Each data point is an average over all the replicates.
"""

import numpy, pandas, pylab, rdd, sympy, sympy.abc, sys

def main(args):
    """
    Entry point.
    """

    files = sys.stdin.readlines()
    replicates = len(files)
    orders = None
    times1, times2, times3 = [], [], []
    edges1, edges2, edges3 = [], [], []
    excess = []
    for file in files:
        file = file.strip()
        data = pandas.read_table(file, sep = " ", header = 0)
        orders = data["m"].tolist()
        times1.append(data["tA"])
        times2.append(data["tB"])
        times3.append(data["tC"])
        edges1.append(data["EA"])
        edges2.append(data["EB"])
        edges3.append(data["EC"])
        excess.append(data["excess"])
    x = numpy.linspace(orders[0], orders[-1], 50)
    times1 = (sum(times1) / replicates).tolist()
    times2 = (sum(times2) / replicates).tolist()
    times3 = (sum(times3) / replicates).tolist()
    edges1 = (sum(edges1) / replicates).tolist()
    edges2 = (sum(edges2) / replicates).tolist()
    edges3 = (sum(edges3) / replicates).tolist()
    excess = (sum(excess) / replicates).tolist()

    # Figure 1: average search time versus permutation order.
    f1 = numpy.poly1d(numpy.polyfit(orders, times1, 3))
    f1_str = sympy.latex(sympy.Poly([numpy.around(c, 4) for c in f1.coeffs], 
                                    sympy.abc.m).as_expr())
    y1 = f1(x)
    f2 = numpy.poly1d(numpy.polyfit(orders, times2, 3))
    f2_str = sympy.latex(sympy.Poly([numpy.around(c, 4) for c in f2.coeffs], 
                                    sympy.abc.m).as_expr())
    y2 = f2(x)
    f3 = numpy.poly1d(numpy.polyfit(orders, times3, 3))
    f3_str = sympy.latex(sympy.Poly([numpy.around(c, 4) for c in f3.coeffs], 
                                    sympy.abc.m).as_expr())
    y3 = f3(x)
    pylab.figure(1, figsize = (7, 4.5), dpi = 500)
    pylab.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    pylab.xlabel(r"permutation order, $m$", fontsize = "x-small")
    pylab.ylabel(r"average search time, $t$ (in $\mu s$)", 
                 fontsize = "x-small")
    pylab.plt.tick_params(labelsize = "x-small")
    pylab.plot(orders, times1, "k^", alpha = 0.6, 
               label = r"$\mathcal{A}: %s$" %(f1_str))
    pylab.plot(x, y1, "k-", linewidth = 1, alpha = 0.6)
    pylab.plot(orders, times2, "k.", alpha = 0.6, 
               label = r"$\mathcal{B}: %s$" %(f2_str))
    pylab.plot(x, y2, "k-", linewidth = 1, alpha = 0.6)
    pylab.plot(orders, times3, "kv", alpha = 0.6, 
               label = r"$\mathcal{C}: %s$" %(f3_str))
    pylab.plot(x, y3, "k-", linewidth = 1, alpha = 0.6)
    leg = pylab.legend(loc = "best", frameon = False)
    pylab.setp(leg.get_texts(), fontsize = "x-small")
    pylab.savefig("search_time.pdf", format = "pdf")
    pylab.close(1)

    # Figure 2: average number of edges traversed versus permutation order.
    f1 = numpy.poly1d(numpy.polyfit(orders, edges1, 2))
    f1_str = sympy.latex(sympy.Poly([numpy.around(c, 4) for c in f1.coeffs], 
                                    sympy.abc.m).as_expr())
    y1 = f1(x)
    f2 = numpy.poly1d(numpy.polyfit(orders, edges2, 2))
    f2_str = sympy.latex(sympy.Poly([numpy.around(c, 4) for c in f2.coeffs], 
                                    sympy.abc.m).as_expr())
    y2 = f2(x)
    f3 = numpy.poly1d(numpy.polyfit(orders, edges3, 2))
    f3_str = sympy.latex(sympy.Poly([numpy.around(c, 4) for c in f3.coeffs], 
                                    sympy.abc.m).as_expr())
    y3 = f3(x)
    pylab.figure(2, figsize = (7, 4.5), dpi = 500)
    pylab.xlabel(r"permutation order, $m$", fontsize = "x-small")
    pylab.ylabel(r"average number of edges traversed, $\mathcal{E}$", 
                 fontsize = "x-small")
    pylab.plt.tick_params(labelsize = "x-small")
    pylab.plot(orders, edges1, "k^", alpha = 0.6, 
               label = r"$\mathcal{A}: %s$" %(f1_str))
    pylab.plot(x, y1, "k-", linewidth = 1, alpha = 0.6)
    pylab.plot(orders, edges2, "k.", alpha = 0.6, 
               label = r"$\mathcal{B}: %s$" %(f2_str))
    pylab.plot(x, y2, "k-", linewidth = 1, alpha = 0.6)
    pylab.plot(orders, edges3, "kv", alpha = 0.6,
               label = r"$\mathcal{C}: %s$" %(f3_str))
    pylab.plot(x, y3, "k-", linewidth = 1, alpha = 0.6)
    leg = pylab.legend(loc = "best", frameon = False)
    pylab.setp(leg.get_texts(), fontsize = "x-small")
    pylab.savefig("edges_traversed.pdf", format = "pdf")
    pylab.close(2)

    # Figure 3: excess of distinghishability over idea radius versus 
    # permutation order.
    pylab.figure(3, figsize = (7, 4.5), dpi = 500)
    pylab.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    pylab.xlabel(r"permutation order, $m$", fontsize = "x-small")
    pylab.ylabel(r"excess of distinguishability over ideal radius, $\delta$", 
                 fontsize = "x-small")
    pylab.plt.tick_params(labelsize = "x-small")
    pylab.plot(orders, excess, "r-", linewidth = 1, alpha = 0.6)
    pylab.savefig("excess.pdf", format = "pdf")
    pylab.close(3)
        
if __name__ == "__main__":
    main(sys.argv[1:])
