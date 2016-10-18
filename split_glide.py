# split_glide.py --- Write the batch file that will run GLIDE sequentially
# on the "tiles" created beforehand using 
# > split -l <tile size> <GLIDE input file> <GLIDE input root>.glide_

# contact: chloe-agathe.azencott@mines-paristech.fr


import glob
import os
import string
import subprocess
import sys

GPULIST = [0]
GPUNUM = len(GPULIST)
BLOCKSIZE = 4096 # <-- Will be passed to the -p flag of GLIDE; number of SNPs in each GPU block

def main(args):
    usage = """python %s <GLIDE input root> <phenotype file> <number of subjects> <t-test threshold> <bash commands file root>
    Create the list of cudalin commands necessary to compute all interactions.
    Write one file per GPU (between 0 and GPUNUM)\n""" % args[0]
    if len(args) != 6:
        print args
        sys.stderr.write(usage)
        sys.exit(0)
    else:
        glideInputRoot = args[1]
        phenoF = args[2]
        try:
            numSubjects = int(args[3])
        except ValueError:
            sys.stderr.write("<number of subjects> should be an integer.\n")
            sys.exit(-1)           
        try:
            tThresh = float(args[4])
        except ValueError:
            sys.stderr.write("<t-test threshold> should be a float.\n")
            sys.exit(-1)
        outputBashF = args[5]

    # Create list of input files 
    inputFilesList = []
    cudaRoot = "%s.glide_*" % glideInputRoot
    for inputF in glob.glob(cudaRoot):
        inputFilesList.append(inputF)
    print len(inputFilesList), "input files detected"

    # Pair up input files
    gList = [open("%s_%s.sh" % (outputBashF, gpu), 'w') for gpu in GPULIST]
    for (i, file1) in enumerate(inputFilesList):
        split1 = file1.split("_")[-1]
        # Number of lines
        proc = subprocess.Popen("wc -l %s" % file1,
                                shell=True,
                                stdout=subprocess.PIPE)
        numSnps1 = int(proc.communicate()[0].split()[0])

        for (j, file2) in enumerate(inputFilesList[i:]):
            split2 = file2.split("_")[-1]
            # Name of the output file                
            outputF = "%s_%s_%s.output" % \
                      (glideInputRoot, split1, split2)
            # Keep going only if output file does not exist
            if not os.path.exists(outputF):
                # Number of lines
                proc = subprocess.Popen("wc -l %s" % file2,
                                        shell=True,
                                        stdout=subprocess.PIPE)
                numSnps2 = int(proc.communicate()[0].split()[0])

                # Command
                gpu = 0 #(i+j) % GPUNUM
                cmd = "./GLIDE -f1 %s -f2 %s -fp %s -n %s -m %s -m2 %s -p %s -t %s -o %s -g %s\n" % \
                      (file1, file2, phenoF, numSubjects, numSnps1, numSnps2,
                       BLOCKSIZE, tThresh, outputF, gpu)
                #gList[gpu].write(cmd)
                gList[0].write(cmd)

    [g.close() for g in gList]
    print gList

if __name__ == "__main__":
    main(sys.argv)
