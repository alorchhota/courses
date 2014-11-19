#!/usr/bin/python

import sys
from Rna import Rna 


# format I/O
if __name__ == '__main__':
    codonMapFile = 'inputs/rna-codon.txt'
    with sys.stdin as inFile:
        lines = inFile.readlines()
        rna = Rna(lines[0].strip(), codonMapFile)
        protein = rna.translate()

    with sys.stdout as outFile:
        outFile.write(protein)


