import sys
import os
import re

if __name__ == '__main__':
    #print(sys.argv)
    inFile = sys.stdin
else:
    workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW3'
    os.chdir(workDir)
    sys.path.append(os.path.abspath(os.getcwd()))
    inputFilePath = 'data/3_longest_shared_motif_1.txt'
    inFile = open(inputFilePath, 'r')

import Fasta
import SuffixTree

# read inputs
dnas = Fasta.parse_fasta(inFile)
inFile.close()


# find the smallest string
def lcp (s1, s2):
    p = ""
    for i in range(min([len(s1), len(s2)])):
        if s1[i] == s2[i]:
            p += s1[i]
        else:
            break
    return p

def traverseAndPrune(n1, n2, l1i, l2i):
    '''
    traverse node n1 & n2
    label matching starts from index l1i & l2i in n1.lab & n2.lab
    '''
    if n1.lab is  not None and n2.lab is not None:
        lab1 = n1.lab[l1i:]
        lab2 = n2.lab[l2i:]
    else:
        lab1 = ""
        lab2 = ""
    if lab1 == lab2:
        '''traverse next common nodes'''
        n1keys = n1.out.keys()
        n2keys = n2.out.keys()
        for n1k in n1keys:
            if n1k not in n2keys:
                del n1.out[n1k]
        for n2k in n2keys:
            if n2k not in n1keys:
                del n2.out[n2k]
        for k in n1.out.keys():
            traverseAndPrune(n1.out[k], n2.out[k], 0, 0)
    else:
        p = lcp(lab1, lab2)
        if len(p) < len(lab1) and len(p) <len(lab2):
            '''mismatch found, prune all children of both nodes'''
            n1.lab = n1.lab[:l1i+len(p)]
            n2.lab = n2.lab[:l2i+len(p)]
            for k in n1.out.keys():
                del n1.out[k]
            for k in n2.out.keys():
                del n2.out[k]
        elif len(p) < len(lab1):
            '''try to move to a child of n2'''
            l1i += len(p)
            n2keys = n2.out.keys()
            for n2k in n2keys:
                if n2k != n1.lab[l1i]:
                    del n2.out[n2k]
            if n1.lab[l1i] in n2.out.keys():
                traverseAndPrune(n1, n2.out[n1.lab[l1i]], l1i, 0)
        else:
            '''try to move to a child of n1'''
            l2i += len(p)
            n1keys = n1.out.keys()
            for n1k in n1keys:
                if n1k != n2.lab[l2i]:
                    del n1.out[n1k]
            if n2.lab[l2i] in n1.out.keys():
                traverseAndPrune(n1.out[n2.lab[l2i]], n2, 0, l2i)
        
            
## create the first suffix tree
st1 = SuffixTree.SuffixTree(dnas[0])
## create suffix tree for other dnas 
## and prune st1 if some branch is not present in the new suffix tree
for i in range(1, len(dnas)):
    st2 = SuffixTree.SuffixTree(dnas[i])
    traverseAndPrune(st1.root, st2.root, 0, 0)

## find the longest substring from remaining common substrings
def longestSubstr(node):
    lab = re.sub('\$', '', node.lab, 0) if node.lab is not None else ""
    if len(node.out) == 0:
        return lab
    # longest substr is the label + longest child subtree 
    childSubs = [re.sub('\$', '', longestSubstr(node.out[key]), 0) for key in node.out]
    childSubs.sort(key=lambda s: len(s), reverse=True)
    lsub = lab + childSubs[0]
    return lsub

lss = longestSubstr(st1.root)
print(lss)