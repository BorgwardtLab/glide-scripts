# rerun_confounder.py --- Run linear regression again,
# on a small number of SNP pairs, including confounders.

# contact: chloe-agathe.azencott@mines-paristech.fr

import numpy as np
import os
import sys

from scipy import stats as st

from utils import readSnp, testPair_logistic as testPair



def main(args):
    usage = """python %s <glideIn> <snps list> <pheno> <confounders> <significant pairs of SNPs> <reevaluated pairs of SNPs>
    Rerun linear regression evaluation on the significant pairs of SNPs.
    E.g.: py rerun_confounder.py mydata_final_clean.glideIn mydata_final_clean.snpNames mydata_samples.sex mydata_final_clean.pheno nydata_clean.sortedpvals mydata_final_clean.sortedpvals.rerun_sex\n""" % args[0]
    if len(args) != 7:
        sys.stderr.write(usage)
        sys.exit(0)

    glideInF = args[1]
    snpsLstF = args[2]
    cfdrsF   = args[3]
    phenoF   = args[4]
    sigPairF = args[5]
    outputF  = args[6]

    threshold = 1e-11

    pheno = np.loadtxt(phenoF)
    confounders = np.loadtxt(cfdrsF)

    f = open(snpsLstF)
    snpsDict = {line.split()[0]:idx for (idx, line) in enumerate(f)}
    f.close()

    newPvalsDict = {} # pvalue, line to write
    with open(sigPairF) as f:
        f.readline()        
        for i, line in enumerate(f):
            ls = line.split()
            try:
                snp1idx = snpsDict[ls[0]]                           
                try:
                    snp2idx = snpsDict[ls[1]]
                    snp1x = readSnp(snp1idx, glideInF) 
                    snp2x = readSnp(snp2idx, glideInF) 

                    pval, pvalStr = testPair(snp1x, snp2x, pheno, confounders)
                    line_to_write = "%s\t%s\n" % ("\t".join(ls), pvalStr)
                                                        
                    if not newPvalsDict.has_key(pval):
                        newPvalsDict[pval] = [line_to_write]
                    else:
                        newPvalsDict[pval].append(line_to_write)

                except KeyError:
                    print "Didn't find %s in SNPs list" % ls[3]
                    sys.exit(-1)
            except KeyError:
                print "Didn't find %s in SNPs list" % ls[0]
                sys.exit(-1)
    f.close()

    with open(outputF, 'w') as g:
        g.write("SNP1\tchr1\tpos1\tSNP2\tchr2\tpos2\t")
        g.write("t-testGLIDE\tpvalGLIDE(intercept)\tpvalGLIDE(x1)\tpvalGLIDE(x2)\tpvalGLIDE(x1:x2)")
        g.write("\tpval(intercept)\tpval(x1)\tpval(x2)\tpval(confounder(s)\tpval(x1:x2)\n")

        sortedPvals = newPvalsDict.keys()
        sortedPvals.sort()
        for pval in sortedPvals:
            for line_to_write in newPvalsDict[pval]:
                g.write(line_to_write)
    g.close()


        
        

if __name__ == "__main__":
    main(sys.argv)
