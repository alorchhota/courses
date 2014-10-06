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
    inputFilePath = 'data/5_edit_0.txt'
    inFile = open(inputFilePath, 'r')

import Fasta

# read inputs
proteins = Fasta.parse_fasta(inFile)
inFile.close()

protein1 = proteins[0]
protein2 = proteins[1]

sigma = 0

## create scores and backtrack arrays
len1 = len(protein1)
len2 = len(protein2)
scores = [[0]*(len2+1) for x in range(len1+1)]
backtracks = [[-1]*(len2+1) for x in range(len1+1)]

## put first row and first column of scores and backtrack arrays
for i in range(1,len1+1):
    scores[i][0] = scores[i-1][0] + sigma
    backtracks[i][0] = 1
    
for j in range(1,len2+1):
    scores[0][j] = scores[0][j-1] + sigma
    backtracks[0][j] = 2

## update scores in greedy approach
for i in range(1,len1+1):
    for j in range(1,len2+1):
        topScore = scores[i-1][j] + sigma
        leftScore = scores[i][j-1] + sigma
        diagScore = scores[i-1][j-1] + int(protein1[i-1]==protein2[j-1])
        candidateScores = [diagScore, topScore, leftScore]
        scores[i][j] = max(candidateScores)
        backtracks[i][j] = candidateScores.index(scores[i][j])

## max score
maxScore = scores[len1][len2]

# backtrack
lcs = []
i = len1
j = len2
while i>0 and j>0:
    if backtracks[i][j] == 0:
        if protein1[i-1] == protein2[j-1]:
            lcs.append(protein1[i-1])
        i -= 1
        j -= 1
    elif backtracks[i][j] == 1:
        i -= 1
    elif backtracks[i][j] == 2:
        j -= 1

lcs1 = ''.join([lcs[i] for i in range(len(lcs)-1,-1,-1)])
print(lcs1)


#print('\n'.join( '\t'.join([str(backtracks[i][j]) for j in range(len2+1)]) for i in range(len1+1)))

#print('\n'.join( '\t'.join([str(scores[i][j]) for j in range(len2+1)]) for i in range(len1+1)))