# rerun_confounder_h5.py --- Run linear regression again,
# on a small number of SNP pairs, including confounders.
# Read genotypes from .h5 file (created with plink2h5.py).

# contact: chloe-agathe.azencott@mines-paristech.fr


import numpy as np
import os
import sys
import tables


from utils import testPair_logistic as testPair


def main(args):
    usage = """python %s <glideIn> <snps list> <pheno> <confounders> <significant pairs of SNPs> <reevaluated pairs of SNPs>
    Rerun linear regression evaluation on the significant pairs of SNPs.
    E.g.: py rerun_confounder.py mydata_final_clean.h5 mydata_final_clean.bim mydata_samples.sex mydata_final_clean.phenoGlide mydata_final_clean.sortedpvals mydata_final_clean.sortedpvals.rerun_sex\n""" % args[0]
    if len(args) != 7:
        sys.stderr.write(usage)
        sys.exit(0)

    h5fname = args[1]
    bimF = args[2]
    cfdrsF   = args[3]
    phenoF   = args[4]
    sigPairF = args[5]
    outputF  = args[6]

    threshold = 1e-11

    pheno = np.loadtxt(phenoF)
    confounders = np.loadtxt(cfdrsF)


    with open(bimF) as f:
        snpsDict = {('%s\t%s\t%s' % \
                     (line.split()[1], line.split()[0], line.split()[3])):idx \
                    for (idx, line) in enumerate(f)}

    newPvalsDict = {} # pvalue, line to write
    with tables.openFile(h5fname) as h5f, open(sigPairF) as f:
        print "%d SNPs x %d samples" % (h5f.root.genotype.shape[0], h5f.root.genotype.shape[1])
        f.readline()        
        for i, line in enumerate(f):
            ls = line.split()
            try:
                snp1idx = snpsDict['%s\t%s\t%s' % (ls[0], ls[1], ls[2])]
                try:
                    snp2idx = snpsDict['%s\t%s\t%s' % (ls[3], ls[4], ls[5])]
                    snp1x = np.ma.masked_values(h5f.root.genotype[snp1idx], 3)
                    snp2x = np.ma.masked_values(h5f.root.genotype[snp2idx], 3)

                    pval, pvalStr = testPair(snp1x, snp2x, pheno, confounders)
                    line_to_write = "%s\t%s\n" % ("\t".join(ls), pvalStr)
                                                        
                    if not newPvalsDict.has_key(pval):
                        newPvalsDict[pval] = [line_to_write]
                    else:
                        newPvalsDict[pval].append(line_to_write)

                except KeyError:
                    print "Didn't find %s in .bim file" % ls[3]
                    sys.exit(-1)
            except KeyError:
                print "Didn't find %s in .bim file" % ls[0]
                sys.exit(-1)
    f.close()
    h5f.close()

    with open(outputF, 'w') as g:
        g.write("SNP1\tchr1\tpos1\tSNP2\tchr2\tpos2\t")
        g.write("t-testGLIDE\tpvalGLIDE(intercept)\tpvalGLIDE(x1)\tpvalGLIDE(x2)\tpvalGLIDE(x1:x2)")
        g.write("\tpval(intercept)\tpval(x1)\tpval(x2)\tpval(confounder(s))\tpval(x1:x2)\n")

        sortedPvals = newPvalsDict.keys()
        sortedPvals.sort()
        for pval in sortedPvals:
            for line_to_write in newPvalsDict[pval]:
                g.write(line_to_write)
    g.close()
        
        

if __name__ == "__main__":
    main(sys.argv)
