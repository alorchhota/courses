import sys
import os
import re
import urllib
import itertools


########## settings to work in personal computer ########
if __name__ != '__main__':
    workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW4'
    os.chdir(workDir)
    sys.path.append(os.path.abspath(os.getcwd()))

################################################################
############################ part-1#### ######################## 
################################################################

# Download the file containing the reads to "reads.fa" in current directory
inputFilePath = 'reads.fa'
urllib.urlretrieve("http://www.cs.jhu.edu/~langmea/resources/f2014_hw4_reads.fa", inputFilePath)
inFile = open(inputFilePath, 'r')

outFile = open('overlaps.txt', 'w')


def parse_fasta(fh):
    ''' Parse FASTA into a dictionary '''
    fa = {}
    name = None
    # Part 1: compile list of lines for each sequence
    for ln in fh:
        if ln[0] == '>':  # new sequence
            name = ln[1:].split()[0]
            fa[name] = []
        else:
            # append nucleotides to current sequence
            fa[name].append(ln.rstrip())
    # Part 2: join lists into strings
    for name, nuc_list in fa.iteritems():
        fa[name] = ''.join(nuc_list)  # join into one long string
    return fa


def make_kmer_table(seqs, k):
    ''' Given dictionary (e.g. output of parse_fasta) and integer k,
        return a dictionary that maps each k-mer to the set of names
        of reads containing the k-mer. '''
    table = {}  # maps k-mer to set of names of reads containing k-mer
    for name, seq in seqs.iteritems():
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer not in table:
                table[kmer] = set()
            table[kmer].add(name)
    return table

def suffixPrefixMatch(str1, str2, min_overlap):
    ''' Returns length of longest suffix of str1 that is prefix of
        str2, as long as that suffix is at least as long as min_overlap. '''
    if len(str2) < min_overlap: return 0
    str2_prefix = str2[:min_overlap]
    str1_pos = -1
    while True:
        str1_pos = str1.find(str2_prefix, str1_pos + 1)
        if str1_pos == -1: return 0
        str1_suffix = str1[str1_pos:]
        if str2.startswith(str1_suffix): return len(str1_suffix)



kmer_tables = {}

def getKmerTable(k):
    ''' Lazy initialization of k-mer table construction.
        If the k-mer table is present, returns it.
        Ohterwise, build it and then returns it.'''
    if k not in kmer_tables:
        kmer_tables[k] = make_kmer_table(reads, k)
    return kmer_tables[k]
    

## read fasta file
reads = parse_fasta(inFile)


################################################################
############################ part-2 ############################ 
################################################################



############## build an overlap graph #######
reads_to_process = set([rn for rn in reads.keys()])
overlap_graph = []
for k in range(100,39,-1):
    kmer_table = getKmerTable(k)
    processed_read = set()
    for r1 in reads_to_process:
        suffix = reads[r1][-k:]
        # candidate read must have this suffix as a kmer
        candidate_reads = [r2 for r2 in kmer_table[suffix]
                              if reads[r2][:k]==suffix and r1 != r2]
        # if exactly one match found, report it 
        if len(candidate_reads) == 1:
            overlap_graph.append((r1, candidate_reads[0], k))
        # any match found means the read is processed.
        if len(candidate_reads) >= 1:
            processed_read.add(r1)
    # removed processed reads for next round
    reads_to_process = reads_to_process.difference(processed_read)

# save overlap graph
outFile.write('\n'.join([' '.join([str(item) for item in edge])  for edge in overlap_graph]) )

inFile.close()
outFile.close()


#############################################################
############### part-3 ######################################
#############################################################

inFile = open('overlaps.txt','r')
outFile = open('unitigs.txt', 'w')

############# read right buddy inputs from file #########
def parse_overlap_grap(fh):
    ''' Parse overlap file '''
    rb = {}
    # for each line, add an entry in rightBuddy
    for ln in fh:
        #print(ln)
        if len(ln.strip())==0:
            continue
        items = re.split(' ', ln)
        #print(items)
        rb[items[0]] = {'right':items[1], 'score':int(items[2])}
    return rb

rightBuddy = parse_overlap_grap(inFile)


############ construct left buddy from right buddy #########
leftBuddy = {}
for lnode, rb in rightBuddy.iteritems():
    rnode = rb['right']
    score = rb['score']
    if rnode not in leftBuddy:
        leftBuddy[rnode] = {'left':lnode, 'score':score}
    elif score > leftBuddy[rnode]['score']:
        # better buddy found
        leftBuddy[rnode] = {'left':lnode, 'score':score}
    elif score == leftBuddy[rnode]['score']:
        # multiple best buddy, mark for deletion
        leftBuddy[rnode] = {'left':None, 'score':score}


########## construct unitigs by joining mutual best buddies ##########
unitigs = []

def getUnitigIndexToJoin(node):
    '''Given node name, returns the index of the contig node belongs to.'''
    # nodes can be joined either at start or end
    for idx in range(len(unitigs)):
        if(unitigs[idx][0]==node or unitigs[idx][-1]==node):
            return idx
    return None

def join(node1, node2):
    '''joins two nodes to build unitig.
       the joined contigs are saved in unitigs variable'''
    global unitigs
    idx1 = getUnitigIndexToJoin(node1)
    idx2 = getUnitigIndexToJoin(node2)
    if idx1 is not None and idx2 is not None:
        # merge two different unitigs into one
        unitigs[idx1] = unitigs[idx1] + unitigs[idx2]
        unitigs = unitigs[:idx2] + unitigs[idx2+1:]
        return 
    if idx1 is not None:
        # node1 must be at the end of the contig
        # because, node2 will be right to node1
        unitigs[idx1] = unitigs[idx1] + [node2]
        return
    if idx2 is not None:
        # node2 must be at the start of the contig
        # because, node1 will be left to node2
        unitigs[idx2] = [node1] + unitigs[idx2]
        return
    # new contig found
    unitigs.append([node1, node2])
    return

# join possible nodes
for lnode, rb in rightBuddy.iteritems():
    rnode = rb['right']
    if rnode in leftBuddy and leftBuddy[rnode]['left']==lnode:
        join(lnode, rnode)


############ save unitigs in a file ############ 
#print(len(unitigs))
for idx in range(len(unitigs)):
    utig = unitigs[idx]
    outFile.write('START UNITIG ' + str(idx+1) + ' ' + utig[0] + '\n')
    for ui in range(1,len(utig)):
        outFile.write(' ' + str(utig[ui]) + ' ' + str(rightBuddy[utig[ui-1]]['score']) + '\n')
    outFile.write('END UNITIG ' + str(idx+1) + '\n')

inFile.close()
outFile.close()


################################################################
############################ part-4 ############################ 
################################################################

genomeOutFile = open('solution.fa', 'w')

def genomeOfUnitig(utig):
    ''' Given unitigs, returns genome'''
    gen = reads[utig[0]]
    for i in range(0, len(utig)-1):
        #read1 = reads[utig[i-1]]
        read2 = reads[utig[i+1]]
        matchScore = rightBuddy[utig[i]]['score']
        gen = gen + read2[matchScore:]
    return gen

utigGenomes = [genomeOfUnitig(utig) for utig in unitigs]

### explore pairwise suffix-prefix match length
### and use pairs with long matches.
perms = itertools.permutations(range(len(unitigs)), 2)
pairwiseSuffixPrefixMatch = [(p[0], p[1], suffixPrefixMatch(utigGenomes[p[0]], utigGenomes[p[1]], min_overlap=1)) for p in perms]
pairwiseSuffixPrefixMatch.sort(key=lambda x: x[2], reverse=True) 

# the longest matches are as follows:
# utig1 - utig2 - overlap
#   1   -   2   - 99
#   3   -   2   - 99
#   2   -   0   - 98
#   2   -   1   - 98
# So, the unitigis are ordered as (3-2-1-2-0) in the final genome.

genome = utigGenomes[3] + utigGenomes[2][99:] + utigGenomes[1][98:] + utigGenomes[2][99:] + utigGenomes[0][98:] 

#print(len(genome))

############# write genome in fasta format ######
def write_solution(genome, per_line=60, out=sys.stdout):
    offset = 0
    out.write('>solution\n')
    while offset < len(genome):
        nchars = min(len(genome) - offset, per_line)
        line = genome[offset:offset+nchars]
        offset += nchars
        out.write(line + '\n')

write_solution(genome, per_line=60, out=genomeOutFile)

genomeOutFile.close()

