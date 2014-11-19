#!/usr/bin/python

import sys
import re
from Rna import Rna 

def locations(dna, sub):
    return [m.start()+1 for m in re.finditer('(?=' + sub + ')', dna)]


# format I/O
if __name__ == '__main__':
    with sys.stdin as inFile:
        lines = inFile.readlines()
        dna = lines[0].strip()
        sub = lines[1].strip()
        loc = locations(dna, sub)

    with sys.stdout as outFile:
        outFile.write(' '.join([str(l) for l in loc]))


