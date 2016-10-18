# transpose.py -- Transpose .raw file into .glideIn file.

# contact: chloe-agathe.azencott@mines-paristech.fr

import os
import sys

CHUNK = 25000 # number of elements that will be transposed at once.

def main(args):
    usage = """python %s <plink raw root>
    Transpose .raw file (PLINK format)) into .glideIn file
    .raw: numIndividuals lines x numSNPs columns
    .glideIn: numSNPs lines x numIndividuals columns
    """ % args[0]

    if len(args) != 2:
        print usage
        sys.exit(-1)
        
    proot = args[1]

    g = open("%s.glideIn" % proot, 'w')

    iChunk = 0
    keepGoing = True
    while keepGoing:
        print "chunk #%d..." % (iChunk+1)
        f = open("%s.raw" % proot)
        hdr = f.readline()
        
        numSnps = len(hdr.split())-6
        iMin = iChunk*CHUNK+6
        iMax = min((iChunk+1)*CHUNK+6, numSnps+6)
        print "  ", iMin, iMax
        if iMax == numSnps+6:
            keepGoing = False
        
        snps = [[] for i in range(iMax-iMin)]
        for line in f:
            for (i, snpVal) in enumerate(line.split()[iMin:iMax]):
                snps[i].append(snpVal)
        f.close()

        for snpVals in snps:
            g.write("%s\n" % " ".join(snpVals))

        iChunk += 1



if __name__ == "__main__":
    main(sys.argv)

    

