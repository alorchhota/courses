class SpacedSeedIndexHash(object):
    
    def __init__(self, t, seed):
        """ Create index, extracting subsequences based on seed """
        self.t = t
        self.seed = seed
        self.index = {}
        ln = len(seed)
        for i in range(0, len(t)-ln+1):
            subseq = ''.join([t[i+si] for si in range(ln) if seed[si]=='1'])
            if subseq in self.index:
                self.index[subseq].append(i) # subseq already in dictionary
            else:
                self.index[subseq] = [i] # add to dictionary
    def query(self, p):
        """ Return candidate alignments for p """
        return self.index.get(p, [])

class Aligner(object):
    def __init__(self, t,seed):
        self.seed = seed
        self.indexHash = SpacedSeedIndexHash(t,seed)
        self.indexLen = seed.count('1')
    def find(self, P, maxMismatch = 0):
        # divide p
        return 0 
        
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


# Set working directory to work in eclipse
workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW2'
import os
os.chdir(workDir)

# inputs
inputPath = 'data/complete_works.txt'
with open(inputPath, 'r') as fh:
    T = fh.read()
lenT = len(T)

# create inverted index hash
seed = '10101010101'
ih = SeedIndexHash(T, seed)

#P = 'achievements'
#P = 'acquaintance'
P = 'remembrances'

n = len(P)

p1 = ''.join([P[i] for i in range(0,len(P),2)])
p2 = ''.join([P[i] for i in range(1,len(P),2)])

p1_idx = ih.query(p1)
p2_idx = ih.query(p2)


hamm1 = [(idx, hamm(T[idx:idx+n], P))
         for idx in p1_idx 
         if idx <= lenT-n]
hamm2 = [(idx-1, hamm(T[idx-1:idx-1+n], P)) 
         for idx in p2_idx
         if idx >= 1]

exactMatchIdx1 = [idx for (idx, d) in hamm1 if d==0]
exactMatchIdx2 = [idx for (idx, d) in hamm2 if d==0]
distinctMatchIndexes = set(exactMatchIdx1 + exactMatchIdx2)

approxMatchIdx1 = [idx for (idx, d) in hamm1 if d<=1]
approxMatchIdx2 = [idx for (idx, d) in hamm2 if d<=1]
distinctApproxMatchIndexes = set(approxMatchIdx1 + approxMatchIdx2)

numExactMatch = len(distinctMatchIndexes)
numApproxMatch = len(distinctApproxMatchIndexes)
num1Mismatch = numApproxMatch - numExactMatch

numComparisons = len(p1_idx) + len(p2_idx)
specificity = (numApproxMatch + 0.0) / numComparisons

print(numExactMatch)
print(num1Mismatch)
print(specificity)