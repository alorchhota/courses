#!/usr/bin/python

import sys
import re

def transcribe(dna):
    """ transcription from dna to rna """
    return re.sub('T', 'U', dna)

if __name__ == '__main__':
    with sys.stdin as inFile:
        lines = inFile.readlines()
        dna = lines[0]
        rna = transcribe(dna)

        with sys.stdout as outFile:
            outFile.write(rna)
