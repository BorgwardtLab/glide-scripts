# compute_block_to_index.py -- Compute the dictionary (stored in CONST.py)
# used to convert between tile IDs and SNP indices.

# contact: chloe-agathe.azencott@mines-paristech.fr

import string
import sys

TILESIZE = 24576 # <-- How many SNPs each sequentially issued GLIDE command will deal with
#TILESIZE = 10000

def main():
    # Convert block names to offset indices in the list of SNPs
    allTheLetters = string.lowercase
    blockToIndex = {}
    for leftLetterIdx in range(2): # compute just for a and b
        leftLetter = allTheLetters[leftLetterIdx]
        for rightLetterIdx in range(len(allTheLetters)):
            rightLetter = allTheLetters[rightLetterIdx]
            blockIdx = (len(allTheLetters) * leftLetterIdx + rightLetterIdx) * TILESIZE
            blockToIndex['%s%s' % (leftLetter, rightLetter)] = blockIdx
    print blockToIndex


if __name__ == "__main__":
    main()
