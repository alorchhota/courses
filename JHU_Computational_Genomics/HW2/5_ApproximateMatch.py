class IndexHash(object):
    
    def __init__(self, t, ln, ival=1):
        """ Create index, extracting substrings of length 'ln' """
        self.t = t
        self.ln = ln
        self.ival = ival
        self.index = {}
        for i in range(0, len(t)-ln+1):
            substr = t[i:i+ln]
            if substr in self.index:
                self.index[substr].append(i) # substring already in dictionary
            else:
                self.index[substr] = [i] # add to dictionary
    def query(self, p):
        """ Return candidate alignments for p """
        return self.index.get(p[:self.ln], [])

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
#workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW2'
workDir = '/Users/ashis/Desktop/work/github/courses/JHU_Computational_Genomics/HW2'
import os
os.chdir(workDir)

# inputs
inputPath = 'data/complete_works.txt'
with open(inputPath, 'r') as fh:
    T = fh.read()
lenT = len(T)

# create inverted index hash
L = 6   # length of substr
ih = IndexHash(T, L)

#P = 'achievements'
#P = 'acquaintance'
P = 'remembrances'

p1 = P[0:L]
p2 = P[L:2*L]

p1_idx = ih.query(p1)
p2_idx = ih.query(p2)


hamm1 = [(idx, hamm(T[idx+L:idx+2*L], p2))
         for idx in p1_idx 
         if idx <= lenT-2*L]
hamm2 = [(idx-L, hamm(T[idx-L:idx], p1)) 
         for idx in p2_idx
         if idx >= L]

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