import sys
import os

if __name__ == '__main__':
    inputPath = sys.argv[1]
    P = sys.argv[2]
else:
    workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW2'
    os.chdir(workDir)
    inputPath = 'data/complete_works.txt'
    P = 'achievements'
    #P = 'acquaintance'
    #P = 'remembrances'

if(len(P)!=12):
    raise RuntimeError('Length of pattern must be 12.')

class SeedIndexHash(object):
    
    def __init__(self, t, seed):
        """ Create index, extracting subsequences based on seed """
        self.t = t
        self.seed = seed
        self.index = {}
        ln = len(seed)
        one_indexes = [si for si in range(ln) if seed[si]=='1']
        for i in range(0, len(t)-ln+1):
            subseq = ''.join([t[i+si] for si in one_indexes])
            if subseq in self.index:
                self.index[subseq].append(i) # subseq already in dictionary
            else:
                self.index[subseq] = [i] # add to dictionary
    def query(self, p):
        """ Return candidate alignments for p """
        return self.index.get(p, [])

def hamm(s1, s2):
    ''' returns hamming distance between two strings'''
    l1 = len(s1)
    l2 = len(s2)
    #
    if l1 != l2:
        raise RuntimeError('unequal length!')
    #
    mismatch = [s1[i]!=s2[i] for i in range(l1)]
    d = sum(mismatch)
    return d


# inputs
with open(inputPath, 'r') as fh:
    T = fh.read()
lenT = len(T)

# create inverted index hash
seed = '10101010101'
ih = SeedIndexHash(T, seed)


n = len(P)

# split p
p1 = ''.join([P[i] for i in range(0,len(P),2)])
p2 = ''.join([P[i] for i in range(1,len(P),2)])

# find splitted-p matches from index hash
p1_idx = ih.query(p1)
p2_idx = ih.query(p2)

# calculate distance between T and P
hamm1 = [(idx, hamm(T[idx:idx+n], P))
         for idx in p1_idx 
         if idx <= lenT-n]
hamm2 = [(idx-1, hamm(T[idx-1:idx-1+n], P)) 
         for idx in p2_idx
         if idx >= 1]

# distance 0 means, exact matches
exactMatchIdx1 = [idx for (idx, d) in hamm1 if d==0]
exactMatchIdx2 = [idx for (idx, d) in hamm2 if d==0]
distinctMatchIndexes = set(exactMatchIdx1 + exactMatchIdx2)

# distance <= 1 means, all approximate matches
approxMatchIdx1 = [idx for (idx, d) in hamm1 if d<=1]
approxMatchIdx2 = [idx for (idx, d) in hamm2 if d<=1]
distinctApproxMatchIndexes = set(approxMatchIdx1 + approxMatchIdx2)

numExactMatch = len(distinctMatchIndexes)
numApproxMatch = len(distinctApproxMatchIndexes)
num1Mismatch = numApproxMatch - numExactMatch

numComparisons = len(p1_idx) + len(p2_idx)
specificity = (numApproxMatch + 0.0) / numComparisons

# print output
sys.stdout.write("Exact match: " + str(numExactMatch) + '\n')
sys.stdout.write("1 mismatch: " + str(num1Mismatch) + '\n')
sys.stdout.write("Specificity: " + str(specificity) + '\n')
