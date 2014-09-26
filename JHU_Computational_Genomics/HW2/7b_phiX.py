import sys
import os

if __name__ == '__main__':
    print(sys.argv)
    fastaInputPath = sys.argv[1]
    fastqInputPath = sys.argv[2]
else:
    workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW2'
    os.chdir(workDir)
    fastaInputPath = 'data/phix.fa.txt'
    fastqInputPath = 'data/phix_reads.fastq.txt'


def phred33_to_q(qual):
    """ Turn Phred+33 ASCII-encoded quality into Phred-scaled integer """
    return ord(qual)-33
    
def q_to_phred33(Q):
    """ Turn Phred-scaled integer into Phred+33 ASCII-encoded quality """
    return chr(Q + 33)
    
def q_to_p(Q):
    """ Turn Phred-scaled integer into error probability """
    return 10.0 ** (-0.1 * Q)
    
def p_to_q(p):
    """ Turn error probability into Phred-scaled integer """
    import math
    return int(round(-10.0 * math.log10(p)))

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
    # this aligner class works for seeds with 1's only (no 0s).
    def __init__(self, t,seed):
        self.seed = seed
        self.indexHash = SpacedSeedIndexHash(t,seed)
    def find(self, P, maxMismatch = 0):
        # split P
        n = len(self.indexHash.seed)
        n_split = maxMismatch+1
        splits = []
        for i in range(0,n_split):
            splits.append(P[i*n:(i+1)*n])
        # variables needed in the program
        lenT = len(self.indexHash.t)
        lenP = len(P)
        # find matches
        distinctMatchIndexes = set()
        for si in range(0,n_split):
            p = splits[si]
            p_idx = self.indexHash.query(p)
            # calculate hamming distance of possible strings
            hamm_p = [(idx-si*n, hamm(self.indexHash.t[idx-si*n:idx-si*n+lenP], P))
                     for idx in p_idx 
                     if idx-si*n >= 0 and idx-si*n+lenP <= lenT]
            # take allowed hamming distance
            matchIdx = [idx for (idx, d) in hamm_p if d<=maxMismatch]
            for idx in matchIdx:
                distinctMatchIndexes.add(idx)
        #
        return distinctMatchIndexes 



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

with open(fastaInputPath, 'r') as fh:
    dna = parse_fasta(fh)

with open(fastqInputPath, 'r') as fh:
    reads = parse_fastq(fh)

seed = '1' * 30
al = Aligner(dna, seed)
maxMismatch = 4

matchIndexes = []
disagreedIndexes = set()
for r in reads:
    read = r[1]
    matchIdx = al.find(read, maxMismatch)
    matchIndexes.append(matchIdx)
    lenR = len(read)
    for idx in matchIdx:
        refSubstr = dna[idx:idx+lenR]
        disagreed = [idx+ti for ti in range(lenR) if refSubstr[ti] != read[ti]]
        for disIdx in disagreed:
            disagreedIndexes.add(disIdx)

readLength = 150


def findOverlapAt(offset):
    '''
    find reads which have overlap with reference at offset
    return readIndex and matchIndex
    '''
    minMatchIndex = offset-readLength+1
    overlaps = [(ri,mi) for ri in range(len(reads)) for mi in matchIndexes[ri] if minMatchIndex <= mi <= offset]  
    return overlaps 



# check every disagreed index and generate a report whenever necessary
disagreedIndexes = sorted(disagreedIndexes)
for di in disagreedIndexes:
    refNt = dna[di]
    weight = {} # weights of nucleotides
    overlappedReadIndexes = findOverlapAt(di)
    for ov in overlappedReadIndexes:
        ov_read = reads[ov[0]]      # overlapped read 
        nt = ov_read[1][di-ov[1]]   # mismatched nt
        phq = ov_read[2][di-ov[1]]  # phred quality
        q = phred33_to_q(phq)       # quality value
        # update weight of nt
        if(nt in weight):
            weight[nt] = weight[nt] + q 
        else:
            weight[nt] = q
    # skip, if no read overlapped or all match with reference
    if(len(weight)==0 or (len(weight)==1 and weight.keys()[0]==refNt)):
        continue
    # sort nucleotides according to weights descending
    sortedNt = sorted(weight.items(), key=lambda x: -x[1])
    if(sortedNt[0][0] != refNt):
        if len(sortedNt) > 1:
            printTuple = [di, sortedNt[0][0], sortedNt[0][1], sortedNt[1][0], sortedNt[0][1], refNt ]
        else:
            # if only one type of nucleotides found, report refNt as second one with 0 weight
            printTuple = [di, sortedNt[0][0], sortedNt[0][1], refNt, 0, refNt ]
        sys.stdout.write(' '.join([str(item) for item in printTuple]) + '\n')
