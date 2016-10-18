# naive_impute.py -- Impute missing values in .glideIn file,
# using the most frequent value for this phenotype.

# contact: chloe-agathe.azencott@mines-paristech.fr


import os
import sys

from collections import Counter

def main(args):
    usage = """python %s <glideIn root>
    .glideIn: numSNPs lines x numIndividuals columns 
    Produce .glideInImputed file, where missing values (NA)
    have been replaced by the most frequent SNP value for the phenotype.    
    """ % args[0]

    if len(args) != 2:
        print usage
        sys.exit(-1)
        
    proot = args[1]

    f = open("%s.pheno" % proot)
    phenos = [int(line.split()[0]) for line in f]
    f.close()
    pos = set([i for (i, p) in enumerate(phenos) if p == 1])
    neg = set([i for (i, p) in enumerate(phenos) if p == 0])
    print len(pos), "pos", len(neg), "neg"


    f = open("%s.glideIn" % proot)
    g = open("%s.glideInImputed" % proot, 'w')
    for line in f:
        vals = line.split()
        miss = [i for (i, val) in enumerate(vals) if val == 'NA']
        if len(miss):
            pval = Counter([vals[i] for i in pos if vals[i] != 'NA']).most_common(1)[0][0]
            nval = Counter([vals[i] for i in neg if vals[i] != 'NA']).most_common(1)[0][0]

            newvals = [val if val != 'NA' else (pval if i in pos else nval) \
                       for (i, val) in enumerate(vals)]
            g.write("%s\n" % (" ".join(newvals)))
        else:
            g.write(line)
    f.close()
    g.close()

    
if __name__ == "__main__":
    main(sys.argv)

    
