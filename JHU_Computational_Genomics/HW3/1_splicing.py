import sys
import os
import re

if __name__ == '__main__':
    #print(sys.argv)
    inFile = sys.stdin
    codonMapFile = sys.argv[1]
else:
    workDir = '/home/ashis/work/github/courses/JHU_Computational_Genomics/HW3'
    os.chdir(workDir)
    sys.path.append(os.path.abspath(os.getcwd()))
    inputFilePath = 'data/rosalind_splc.txt'
    inFile = open(inputFilePath, 'r')
    codonMapFile = 'inputs/rna-codon.txt'

import Dna;
import Rna;
import Fasta; 
# read inputs

lines = Fasta.parse_fasta(inFile)
inFile.close()

dna = lines[0]
introns = lines[1:]

# splice dna
re_pattern = '|'.join(introns)
splicedDna = re.sub(re_pattern, '', dna, 0)
sdna = Dna.Dna(splicedDna)

# transcribe and translate
rna = Rna.Rna(sdna.transcribe(), codonMapFile=codonMapFile);
protein = rna.translate()
print(protein)
