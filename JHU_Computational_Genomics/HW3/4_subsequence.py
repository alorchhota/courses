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
    inputFilePath = 'data/4_subsequence_1.txt'
    inFile = open(inputFilePath, 'r')

import Fasta

# read inputs
dnas = Fasta.parse_fasta(inFile)
inFile.close()

s = dnas[0]
t = dnas[1]

subseqIndexes = []

si = -1
for ti in range(len(t)):
    tchar = t[ti]
    while 1:
        si += 1
        if s[si] == tchar:
             subseqIndexes.append(si)
             break

print(' '.join([str(i+1) for i in subseqIndexes]))