# qqplot.py --- Plot a Q-Q plot for (single SNP) p-values, from PLINK output.

# contact: chloe-agathe.azencott@mines-paristech.fr


import sys
from math import log10

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc

def main(args):
    usage = """python %s <PLINK association file> <qqplot file>
    Plot a Q-Q plot for the given PLINK output file.
    """ % args[0]

    if len(args) != 3:
        print usage
        sys.exit(0)

    plinkOutF = args[1]
    qplotOutF = args[2]

    
    with open(plinkOutF) as f:
        f.readline()
        observed = [float(line.split()[8]) for line in f]
    observed.sort()
    observed = [-log10(x) for x in observed]

    L = len(observed)
    expected = range(1,L+1)
    L = float(L)
    expected = [-log10(float(x)/L) for x in expected]



    plt.ion()
    rc('text', usetex=True)
    rc('font', family='serif')

    plt.figure(figsize=(6,6))
    pltIdx = 1
    ax = plt.subplot(1, 1, pltIdx)

    plt.plot(expected, observed, linestyle='-', color='#cc0000')
    plt.plot(expected, expected, linestyle='-', color='#cccccc')

    plt.xlabel('Expected (-logP)')
    plt.ylabel('Observed (-logP)')
    plt.xlim([0, 6])
    plt.ylim([0, 6])

    plt.savefig(qplotOutF)
    sys.stdout.write('Q-Q plot saved to %s\n' % qplotOutF)
    plt.clf()
    plt.close()
    


if __name__ == "__main__":
    main(sys.argv)
