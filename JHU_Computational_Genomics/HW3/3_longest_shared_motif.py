import sys
import os
import re

if __name__ == '__main__':
    #print(sys.argv)
    inFile = sys.stdin
else:
    workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW3'
    os.chdir(workDir)
    sys.path.append(os.path.abspath(os.getcwd()))
    inputFilePath = 'data/3_longest_shared_motif_1.txt'
    inFile = open(inputFilePath, 'r')

import Fasta

# read inputs
dnas = Fasta.parse_fasta(inFile)
inFile.close()

# find the smallest string
dnaLen = [len(d) for d in dnas]
smallestLen = min(dnaLen)
smallestIndex = dnaLen.index(smallestLen)
smallestDna = dnas[smallestIndex]

# splice dna
exists = [False]
motif = ''
for curLen in range(smallestLen,0,-1):
    for start in range(smallestLen-curLen+1):
        substr = smallestDna[start:start+curLen]
        exists = [d.find(substr)>=0 for d in dnas]
        if all(exists):
            motif = substr
            break
    if len(motif) > 0:
        break

print(motif)
