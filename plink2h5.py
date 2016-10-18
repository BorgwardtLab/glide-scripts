# plink2h5.py -- Convert binary PLINK files into h5 file.

# contact: chloe-agathe.azencott@mines-paristech.fr


import sys
import tables

from plinkio import plinkfile as pf


def main(args):
    usage = """python %s <plink root> <h5 file>
    Convert binary PLINK files into h5 file.
    E.g.: py plink2h5.py mydata_final_clean mydata_final_clean.h5\n""" % args[0]
    if len(args) != 3:
        sys.stderr.write(usage)
        sys.exit(0)

    plinkRoot = args[1]
    h5fname   = args[2]

    plinkTitle = plinkRoot.split("/")[-1]

    # Read binary PLINK files
    plinkF = pf.open(plinkRoot)    
    numSnps = len(plinkF.get_loci())
    numSamples = len(plinkF.get_samples())
    print "%d SNPs x %d samples" % (numSnps, numSamples) 
    
    # Create the empty array to store genotypes
    atom = tables.Int8Atom()
    h5F  = tables.openFile(h5fname, 'w', title=plinkTitle)
    genotype = h5F.createCArray(h5F.root, 'genotype', atom,
                                (numSnps, numSamples),
                                title='Genotype',
                                filters=tables.Filters(complevel=5,
                                                       complib='blosc'))

    # populate
    for counter, row in enumerate(plinkF):
        genotype[counter,:] = list(row)
        if counter % 10000 == 9999:
            print (counter + 1), 'SNPs read'
    plinkF.close()
    h5F.close()



if __name__ == "__main__":
    main(sys.argv)
