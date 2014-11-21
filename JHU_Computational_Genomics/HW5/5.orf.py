import sys
import os
import re
import csv


########## settings to work in personal computer ########
inFile = sys.stdin
#inFile = open('data/rosalind_orf.txt', 'r')
if len(sys.argv) < 2:
    codonMapFile = 'inputs/rna-codon.txt'
else:
    codonMapFile = sys.argv[1]


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


def dna2rna(dna):
    """ transcription from dna to rna """
    return re.sub('T', 'U', dna)


def reverseComplement(dna):
    complementMap = {'A':'T', 'T':'A', 'G':'C','C':'G'}
    complement = "".join([complementMap[nucleotide] for nucleotide in dna])
    reverseComplement = complement[::-1]
    return reverseComplement

def rna2peptide(rna, codonMap=None):
    if codonMap==None:
        raise ValueError('codonMap must be provided')
    codons = [rna[start:(start+3)] for start in range(0, len(rna)-2, 3)]
    peptide = "".join(["-" if codonMap[c]=="" else codonMap[c] for c in codons])
    return peptide

# read dna
dnaStrs = parse_fasta(inFile)
dna = dnaStrs[dnaStrs.keys()[0]]

# convert dna to rna
rna = dna2rna(dna)

# read codon map
codonMap = None
with open(codonMapFile) as f:
    freader = csv.reader(f, delimiter='\t')
    codonMap = {row[0]:row[1] for row in freader if len(row)>0}

# list of start and stop codons
startCodons = ['AUG']
stopCodons = ['UAG', 'UGA', 'UAA']

# store all distinct peptides
peptides = set()

## search open read frames in forward direction
for start in range(0,3):                # loop for reading frame start
    startCodonPos = [i for i in range(start, len(rna), 3) if rna[i:i+3] in startCodons]
    stopCodonPos = [i for i in range(start, len(rna), 3) if rna[i:i+3] in stopCodons]
    for pos1 in startCodonPos:
        for pos2 in stopCodonPos:
            if pos2 > pos1:
                curRna = rna[pos1:pos2]
                curPeptide = rna2peptide(curRna, codonMap)
                peptides.add(curPeptide)
                break

## search open read frames in reverse direction
dna = reverseComplement(dna)
rna = dna2rna(dna)
for start in range(0,3):                # loop for reading frame start
    startCodonPos = [i for i in range(start, len(rna), 3) if rna[i:i+3] in startCodons]
    stopCodonPos = [i for i in range(start, len(rna), 3) if rna[i:i+3] in stopCodons]
    for pos1 in startCodonPos:
        for pos2 in stopCodonPos:
            if pos2 > pos1:
                curRna = rna[pos1:pos2]
                curPeptide = rna2peptide(curRna, codonMap)
                peptides.add(curPeptide)
                break


## print output
for p in peptides:
    print(p)


inFile.close()
