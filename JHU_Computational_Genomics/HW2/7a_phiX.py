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
    def find(self, P, maxMismatch = 0):
        # split p
        n = len(self.indexHash.seed)
        n_split = maxMismatch+1
        splits = []
        for i in range(0,n_split):
            splits.append(P[i*n:(i+1)*n])
        #
        lenT = len(self.indexHash.t)
        lenP = len(P)
        #
        distinctMatchIndexes = set()
        for si in range(0,n_split):
            p = splits[si]
            p_idx = self.indexHash.query(p)
            #
            hamm_p = [(idx-si*n, hamm(self.indexHash.t[idx-si*n:idx-si*n+lenP], P))
                     for idx in p_idx 
                     if idx-si*n >= 0 and idx-si*n+lenP <= lenT]
            #
            matchIdx = [idx for (idx, d) in hamm_p if d<=maxMismatch]
            for idx in matchIdx:
                distinctMatchIndexes.add(idx)
        #
        return distinctMatchIndexes 


# Set working directory to work in eclipse
workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW2'
import os
os.chdir(workDir)

# # inputs
# inputPath = 'data/complete_works.txt'
# with open(inputPath, 'r') as fh:
#     T = fh.read()
# lenT = len(T)
# 
# # create inverted index hash
# #seed = '111111111111111111111111111111'
# seed = '111111'
# al = Aligner(T, seed)
# #P = 'achievements'
# #P = 'acquaintance'
# P = 'Semembrances'
# maxMismatch = 1
# matchIdx = al.find(P, maxMismatch)
# lenP = len(P)
# hamm_dist = [hamm(T[idx:idx+lenP], P) for idx in matchIdx]
# print('Exact Match: ' + str(hamm_dist.count(0)))
# print('1 Mismatch: ' + str(hamm_dist.count(1)))

# copied from http://nbviewer.ipython.org/gist/BenLangmead/8376306
def parse_fastq(fh):
    """ Parse reads from a FASTQ filehandle.  For each read, we
        return a name, nucleotide-string, quality-string triple. """
    reads = []
    while True:
        first_line = fh.readline()
        if len(first_line) == 0:
            break  # end of file
        name = first_line[1:].rstrip()
        seq = fh.readline().rstrip()
        fh.readline()  # ignore line starting with +
        qual = fh.readline().rstrip()
        reads.append((name, seq, qual))
    return reads

def parse_fasta(fh):
    ''' 
    parse a fasta file.
    return the dna string.
    '''
    lines = fh.readlines()
    t = ''.join([l.strip() for l in lines if not l.startswith('>')])
    return t

## read inputs
fastaInputPath = 'data/phix.fa.txt'
fastqInputPath = 'data/phix_reads.fastq.txt'

with open(fastaInputPath, 'r') as fh:
    dna = parse_fasta(fh)

with open(fastqInputPath, 'r') as fh:
    reads = parse_fastq(fh)

seed = '1' * 30
al = Aligner(dna, seed)
maxMismatch = 4

best_distances = []
for r in reads:
    read = r[1]
    matchIdx = al.find(read, maxMismatch)
    lenR = len(read)
    hamm_dist = [hamm(dna[idx:idx+lenR], read) for idx in matchIdx]
    best_dist = min(hamm_dist) if len(hamm_dist) != 0 else maxMismatch+1
    best_distances.append(best_dist)
    
counts = [best_distances.count(c) for c in range(5)]
print(counts)