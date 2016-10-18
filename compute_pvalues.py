# compute_pvalues.py --- Compute p-values corresponding to
# the t-test scores output by GLIDE.

# contact: chloe-agathe.azencott@mines-paristech.fr

import gzip
import os
import string
import sys

from scipy import stats

GPULIST = [3]
GPUNUM = len(GPULIST)
BLOCKSIZE = 4096 # <-- this is the parameter passed to GLIDE with flag -p 
THREADSIZE = 16  # <-- this is Bx and By in GLIDE's Makefile

from CONST import blockToIndex


def pval(tTest, numObs):
    """
    Return the p-value corresponding to the given tTest score.
    Degree of freedom = numObs - numConstraints
    """
    return 2 * stats.t.cdf(-abs(tTest), (numObs - 4))


def main(args):
    usage = """python %s <GLIDE output> <blocks> <.bim> <number of samples> <p-values file>
    Compute the p-values corresponding to the t-test value of the four regression coefficients.
    Retrieve from the .bim file the SNP ID and position corresponding to each line.\n""" % args[0]
    if len(args) != 6:
        sys.stderr.write(usage)
        sys.exit(0)
    else:
        glideOutF = args[1]
        blocks    = args[2]
        bimF      = args[3]
        outputF   = args[5]
        try:
            numPatients = int(args[4])            
        except ValueError:
            sys.stderr.write("<number of samples> should be an integer!\n")
            sys.exit(0)


    [leftBlock, rightBlock] = blocks.split("_")
    blockIdxLeft = blockToIndex[leftBlock]
    blockIdxRight = blockToIndex[rightBlock]


    with open(bimF) as f:
        snpsList = ['%s\t%s\t%s' % (line.split()[1], line.split()[0], line.split()[3]) \
                    for line in f]


    with open(glideOutF) as f, open(outputF, 'w') as g:
        g.write("SNP1\tchr1\tpos1\tSNP2\tchr2\tpos2\t")
        g.write("t-test\tpvalue(intercept)\tpvalue(x1)\tpvalue(x2)\tpvalue(x1:x2)\n")
        for line in f:
            ls = line.split()
            if len(ls) < 9:
                continue
            tval = float(ls[9])
            if abs(tval) > 5.0:
                # print (blockIdxLeft + BLOCKSIZE * int(ls[0]) + \
                #        THREADSIZE * int(ls[2]) + \
                #        int(ls[4])), len(snpsList)
                snpId1 = snpsList[blockIdxLeft + BLOCKSIZE * int(ls[0]) + \
                                  THREADSIZE * int(ls[2]) + \
                                  int(ls[4])]
                snpId2 = snpsList[blockIdxRight + BLOCKSIZE * int(ls[1]) + \
                                  THREADSIZE * int(ls[3]) + \
                                  int(ls[5])]
                pvalStr = "\t".join(["%s" % pval(float(ls[i]), numPatients) for i in range(6, 10)])
                g.write("%s\t%s\t%s\t%s\n" % \
                        (snpId1, snpId2, tval, pvalStr))

    
if __name__ == "__main__":
    main(sys.argv)
