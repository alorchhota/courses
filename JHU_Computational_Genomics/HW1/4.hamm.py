#!/usr/bin/python

import sys

def hamm(str1, str2):
    l1,l2 = len(str1), len(str2)
    if l1 != l2:
        raise ValueError('length must be same.')

    dist = sum([1 for i in range(0,l1) if str1[i]!=str2[i]])
    return dist


# format I/O
if __name__ == '__main__':
    with sys.stdin as inFile:
        lines = inFile.readlines()
        dna1 = lines[0].strip()
        dna2 = lines[1].strip()
        d = hamm(dna1, dna2)

    with sys.stdout as outFile:
        outFile.write(str(d))




