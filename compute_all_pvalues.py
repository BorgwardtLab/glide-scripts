# compute_all_pvalues.py --- Compute all p-values corresponding to
# the t-test scores output by GLIDE.
# ie. run compute_pvalues.py on all outputs.

# contact: chloe-agathe.azencott@mines-paristech.fr

import glob
import os
import sys
import subprocess

computePvalues = "compute_pvalues.py"
PYTHON = "/usr/local/bin/python2.7"

def main(args):
    usage = """python %s <GLIDE output root> <.bim file> <number of samples>
    For all GLIDE output starting with the given output root:
        Compute the p-values corresponding to the t-test value of the four regression coefficients.
        Also retrieves the snpId corresponding to each line.
        Write results in <GLIDE output root>_<leftBlock>_<rightBlock>.pvals\n""" % args[0]
    if len(args) != 4:
        print usage
        sys.exit(-1)
    glideOutRoot = args[1]
    bimF         = args[2]
    numSamples   = args[3]

    for glideOutF in glob.glob("%s_*.output" % glideOutRoot):
        glidePvalF = "%s.pvals" % glideOutF[:-7]

        gs = glideOutF.split(".")[-2].split("_")
        blocks = "%s_%s" % (gs[-2], gs[-1])
        cmdList = [PYTHON, computePvalues, glideOutF, blocks,
                   bimF, numSamples, glidePvalF]
        cmdList = ["%s" % x for x in cmdList]
        subprocess.call(cmdList)        

        
if __name__ == "__main__":
    main(sys.argv)
