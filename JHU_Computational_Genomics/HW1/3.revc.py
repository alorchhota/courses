#!/usr/bin/python

import sys
from Dna import *

# format I/O
if __name__ == '__main__':
    with sys.stdin as inFile:
        lines = inFile.readlines()
        dna = Dna(lines[0].strip())
        rc = dna.revc()

    with sys.stdout as outFile:
        outFile.write(rc)




