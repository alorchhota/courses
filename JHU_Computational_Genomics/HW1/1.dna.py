#!/usr/bin/python

import sys
import collections

def countNucleotides(dna):
    """ count neuclides in a dna string"""
    counter = collections.Counter(dna)
    return {c:counter[c] for c in 'ACGT'}

# formated I/O
if __name__ == "__main__":
    with sys.stdin as inFile:
        lines = inFile.readlines()
        dna = lines[0]
        freq = countNucleotides(dna)

        with sys.stdout as outFile:
            outStr = ' '.join([str(freq[c]) for c in 'ACGT'])
            outFile.write(outStr)
        

