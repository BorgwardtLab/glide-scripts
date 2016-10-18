### processIntermediate.py --- Pre-process the data for intermediate phenotypes analysis

## Copyright 2013 Chloe-Agathe Azencott
## Machine Learning and Computational Biology Research Group
## MPI Tuebingen (Germany)

import csv
import os
import pandas as pd
import subprocess 
import sys

CHUNK = 10000

plink = os.path.expanduser('~/plink/plink-1.07-mac-intel/plink')

# phenotypesList = ['length', 'unilaterality', 'pulsation', 'intensity',
#                   'aggravation', 'photo', 'phono', 'nausea', 'vomiting']
# phenotypesList = ['length', 'unilaterality', 'intensity',
#                   'aggravation', 'phono', 'nausea', 'vomiting']
phenotypesList = ['pulsation', 'photo']

phenotypeIndexDict = {'length': 9,
                      'unilaterality': 10,
                      'pulsation': 11,
                      'intensity': 12,
                      'aggravation': 13,
                      'photo': 14,
                      'phono': 15,
                      'nausea': 16,
                      'vomiting': 17}

def getRowsToSkip(oldBimF, newBimF):
    """
    List the indices of SNPs from the oldBimF not present in the newBimF
    """
    print oldBimF, newBimF
    with open(newBimF) as f:
        newSnpsList = set([line.split()[1] for line in f])
    with open(oldBimF) as f:
        rows = []
        for (ix, line) in enumerate(f):
            #print ix, line.split()
            if line.split()[1] not in newSnpsList:
                rows.append(ix)
        return rows
        #return [ix for (ix, line) in enumerate(oldBimF) if line.split()[1] not in newSnpsList]
    

def main(args):
    usage = """python %s <data folder> <csv pheno file> <data root> 
    E.g. py %s Dutch_MA All_clinics_phenotypes_Dutch_MA.csv Dutch_MA_final_maf1hwe
    Pre-process the data for the intermediate phenotypes.
    Create:
      .glideInImputed 
      .phenoGlide
      .individuals\n""" % (args[0], args[0])
    if len(args) != 4:
        sys.stderr.write(usage)
        sys.exit(0)
        
    dataF     = '/agbs/agkbshare/data/Migraine_Verneri_20130205/Clinical_studies/%s' % args[1]
    intPhenoF = '%s/%s' % (dataF, args[2])
    root      = args[3]

    glideInF  = '%s/%s.glideInImputed' % (dataF, root)

    
    # Get intermediate phenotypic data
    phenotypeDict = {} # phenotype:[phenotypeDict, phenotypeInds]
    for phenotype in phenotypesList:
        phenotypeDict[phenotype] = [{}, # fam\tind:phenotype
                                   []]
    with open(intPhenoF) as csvf:
        csvrdr = csv.reader(csvf, delimiter='\t', )
        next(csvrdr, None) 
        for row in csvrdr:
            fam = row[0]
            ind = row[1]
            indId = '%s\t%s' % (fam, ind)
            for phenotype in phenotypesList:
                pIdx = phenotypeIndexDict[phenotype]
                phenotypeValue = int(row[pIdx])
                if phenotypeValue and int(phenotypeValue) in [1, 2]:
                    phenotypeDict[phenotype][0][indId] = int(phenotypeValue)-1
                    phenotypeDict[phenotype][1].append(indId)

    for phenotype in phenotypesList:
        sys.stdout.write("%d individuals with %s data\n" % \
                         (len(phenotypeDict[phenotype][1]), phenotype))

    with open('%s/%s.fam' % (dataF, root)) as f:
        allInds = ['%s\t%s' % (line.split()[0], line.split()[1]) for line in f]
        f.close()


    # Process phenotypes
    for phenotype in phenotypesList:
        with open('%s/%s_%s.pheno' % (dataF, root, phenotype), 'w') as f:
            pDict = phenotypeDict[phenotype][0]
            pInds = phenotypeDict[phenotype][1]
            f.write("\n".join(['%s\t%s' % (indId, pDict[indId]) for indId in pInds]))
            f.write("\n")
            f.close()
        sys.stdout.write("Wrote %s/%s_%s.pheno\n" % (dataF, root, phenotype))

        with open('%s/%s_%s.phenoGlide' % (dataF, root, phenotype), 'w') as f:
            f.write("\n".join(['%s' % pDict[indId] for indId in pInds]))
            f.write("\n")
            f.close()
        sys.stdout.write("Wrote %s/%s_%s.phenoGlide\n" % (dataF, root, phenotype))
    
        individualsF = '%s/%s_%s.individuals' % (dataF, root, phenotype)
        with open(individualsF, 'w') as f:
            f.write("\n".join(pInds))
            f.write("\n")
            f.close()
        sys.stdout.write("Wrote %s\n" % individualsF)

        cmdList = [plink, '--noweb', '--bfile', '%s/%s' % (dataF, root),
                   '--maf', '0.01', '--hwe', '1e-6', '--geno', '0.1', '--make-bed',
                   '--keep', individualsF,
                   '--out', '%s/%s_%s' % (dataF, root, phenotype)]
        print " ".join(cmdList)
        subprocess.call(cmdList)        

        rowsToSkip = getRowsToSkip('%s/%s.bim' % (dataF, root),
                                   '%s/%s_%s.bim' % (dataF, root, phenotype))

        individualsIndices = [ix for (ix, ind) in enumerate(allInds) if ind in set(pInds)]
        glideOutF = '%s/%s_%s.glideInImputed' % (dataF, root, phenotype)
        reader = pd.read_table(glideInF, sep=' ', chunksize=CHUNK, header=None, skiprows=rowsToSkip)
        with open(glideOutF, 'w') as fh:
            for (i, chunk) in enumerate(reader):
                chunk[individualsIndices].to_csv(fh, sep=' ', header=False, index=False),
                sys.stdout.write("%d-th chunk of glideIn processed\n" % i)
            fh.close()
        sys.stdout.write("Wrote %s\n" % glideOutF)


if __name__ == "__main__":
    main(sys.argv)
