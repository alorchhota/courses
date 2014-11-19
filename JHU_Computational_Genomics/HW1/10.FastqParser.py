#!/usr/bin/python

import sys
import collections

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

def countNucleotides(dna):
    """ count neuclides in a dna string"""
    counter = collections.Counter(dna)
    return {c:counter[c] for c in 'ACGT'}

    
def summarize_reads(reads):
    if len(reads) == 0:
        raise Exception('Note enough reads')
    nr = len(reads)
    nnt = len(reads[0][1]) 
    summary = []
    for pos in range(nnt):
        dna = ''.join([r[1][pos] for r in reads])
        counts = countNucleotides(dna)
        counts['X'] = nr - sum(counts.values())
        # quality
        phred = [r[2][pos] for r in reads]
        quality = [phred33_to_q(ph) for ph in phred]
        counts['q'] = sum([1 for q in quality if q <20])
        counts['Q'] = nr - counts['q']
        # summary
        summary.append([counts[c] for c in 'ACGTXqQ'])
    # return
    return(summary)
 
#/home/ashis/work/github/genomics-class/data/hq1_reads.fastq.txt
     
# formated I/O
if __name__ == "__main__":
    with sys.stdin as inFile:
        reads = parse_fastq(inFile)
        summary = summarize_reads(reads)

        with sys.stdout as outFile:
            outStr = '\n'.join([' '.join([str(c) for c in s]) for s in summary])
            outFile.write(outStr)
        
