### meta_analysis_intermediate.py --- Run a meta analysis of the pairs of SNPs selected on different migraine datasets,
### for a given intermediate phenotype

## Copyright 2013 Chloe-Agathe Azencott
## Machine Learning and Computational Biology Research Group
## MPI Tuebingen (Germany)

import sys

import numpy as np
import tables


from utils import testPair_logistic_intermediate as testPair 

dataDirRoot = '/agbs/agkbshare/data/Migraine_Verneri_20130205/Clinical_studies/'
rDir = '/is/ei/cazencott/Research/Migraine/Results'


def main(args):
    usage = """python %s <intermediate phenotye>
    Run a meta analysis of the top-ranked pairs for the different cohorts.
    Save results in %s/meta_<intermediate phenotype>.sortedpvals.
    """ % (args[0], rDir)

    if len(args) == 2:
        phenotypeName = args[1]
    else:
        print usage
        sys.exit(0)
    
    # Get the list of pairs to run
    threshold = 1e-12
    sigPairsList = set([])
    for cohort in ['Dutch_MA', 'Dutch_MO', 'Finnish_MA', 'German_MA']:
        with open("%s/%s_%s.sortedpvals" % (rDir, cohort, phenotypeName)) as f:
            f.readline()
            for line in f:
                ls = line.split()
                if float(ls[-1]) > threshold:
                    break
                pair = ['%s_%s_%s' % (ls[0], ls[1], ls[2]),
                        '%s_%s_%s' % (ls[3], ls[4], ls[5])]
                pair.sort()
                pair = "/".join(pair)
                sigPairsList.add(pair)
        f.close()
    sigPairsList = list(sigPairsList)

    # Only keep pairs of SNPs that appear in all studies
    for cohort in ['Dutch_MA', 'Dutch_MO', 'Finnish_MA', 'German_MA']:
        f = open('%s/%s/%s_clean_%s.bim' % (dataDirRoot, cohort, cohort, phenotypeName))
        if not allSnps:
            allSnps = set([line.split()[1] for line in f])
        else:
            allSnps = allSnps.intersection(set([line.split()[1] for line in f]))
        f.close()
    
    print len(allSnps), "SNPs appear in all studies"

    sigPairsList = [x for x in sigPairsList if \
                    ((x.split("/")[0].split("_")[0] in allSnps) \
                     and (x.split("/")[1].split("_")[0] in allSnps))]
    print len(sigPairsList), "pairs to run."


    # Test pairs for each study
    sigPairsDict = {pair:[] for pair in sigPairsList} 
    for cohort in ['Dutch_MA', 'Dutch_MO', 'Finnish_MA', 'German_MA']:
        print "Testing for %s" % cohort
        dataDir =  '%s/%s' % (dataDirRoot, cohort)
        phenoF  = '%s/%s_clean_%s.phenoGlide' % (dataDir, cohort, phenotypeName)
        bimF    = '%s/%s_clean_%s.bim' % (dataDir, cohort, phenotypeName)
        h5fname = '%s/%s_clean_%s.h5' % (dataDir, cohort, phenotypeName)
    
        pheno = np.loadtxt(phenoF)

        with open(bimF) as f:
            snpsDict = {('%s_%s_%s' % \
                         (line.split()[1], line.split()[0], line.split()[3])):idx \
                        for (idx, line) in enumerate(f)}
        f.close()

        with tables.openFile(h5fname) as h5f:
            for pair in sigPairsList:
                pairSplit = pair.split("/")
                snp1idx = snpsDict[pairSplit[0]]                           
                snp2idx = snpsDict[pairSplit[1]]
                snp1x = np.ma.masked_values(h5f.root.genotype[snp1idx], 3)
                snp2x = np.ma.masked_values(h5f.root.genotype[snp2idx], 3)
                sigPairsDict[pair].append(testPair(snp1x, snp2x, pheno))
        h5f.close()


    # Combine into a meta-analysis
    pvalsDict = {}
    for pair, outputs in sigPairsDict.iteritems():
        betas = np.array([x[0] for x in outputs])
        sesqs = np.array([x[1] for x in outputs])

        betas[np.where(betas==0)] = 1e-10
        sesqs[np.where(sesqs<=0)] = 1e-10
        sesqinvs = 1./sesqs
        
        betam = np.sum(betas*sesqinvs)/np.sum(sesqinvs)
        sesqm = np.sqrt(1./np.sum(sesqinvs))
        zmeta = betam/sesqm
        pval  = 2 * st.norm.sf(zmeta) # two-sided

        try:
            pvalsDict[pval].append([pair, zmeta])
        except KeyError:
            pvalsDict[pval] = [[pair, zmeta]]
            
    # sort p-values and save       
    pvals = pvalsDict.keys()
    pvals.sort()
    with open('%s/meta_logistic_%s.sortedpvals' % (rDir, phenotypeName), 'w') as f:
        f.write("SNP1 & chr & pos & SNP2 & chr & pos & meta z-score & meta pval \n")
        for pval in pvals:
            for [pair, zmeta] in pvalsDict[pval]:
                [snp1idx, snp2idx] = pair.split("/")
                snp1 = snp1idx
                snp2 = snp2idx
                f.write("%s & %s & %.2e & %.2e\n" % \
                        (" & ".join(snp1.split("_")),
                         " & ".join(snp2.split("_")),
                         zmeta, pval))
        f.close()
     

if __name__ == "__main__":
    main(sys.argv)
