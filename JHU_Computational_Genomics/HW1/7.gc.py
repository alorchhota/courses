#!/usr/bin/python

import sys
from Dna import *

# format I/O
if __name__ == '__main__':
    with sys.stdin as inFile:
        lines = inFile.readlines()
        lines = [l.strip() for l in lines]

    # create a dictionary of dnas
    starts = [i for i in range(0,len(lines)) if lines[i].startswith(">")]
    keys = [lines[i][1:].strip() for i in starts]
    dnas = [''.join(lines[(starts[i]+1):starts[i+1]])  for i in range(len(starts)-1)]
    dnas = dnas + [''.join(lines[(starts[-1]+1):])]
    seqs = dict(zip(keys, dnas))
    
    ## calculate gc-content and find max
    gcs = [(k,Dna(v).gc()) for k,v in seqs.items()]
    maxGC = max(gcs, key=lambda x:x[1])

    with sys.stdout as outFile:
        outFile.write('\n'.join([str(item) for item in maxGC]))




