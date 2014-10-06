import sys
import os
import itertools

if __name__ == '__main__':
    inFile = sys.stdin
else:
    workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW3'
    os.chdir(workDir)
    sys.path.append(os.path.abspath(os.getcwd()))
    inputFilePath = 'data/2_permutations-1.txt'
    inFile = open(inputFilePath, 'r')

# read input
lines = inFile.readlines()
inFile.close()
n = int(lines[0].strip())

# number of permutations
nperm = 1
for i in range(1,n+1):
    nperm *= i
print(nperm)

# generate permutations and print
items = [i for i in range(1,n+1)]
perm = itertools.permutations(items)
for p in perm:
    print(' '.join([str(item) for item in p]))
