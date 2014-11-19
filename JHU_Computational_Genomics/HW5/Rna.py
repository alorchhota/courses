import csv

class Rna:
    codonMap = None

    def constructCodonMap(self, mapFilePath):
        """ maps codons from a text file
        map file should be a tab delimited file."""

        with open(mapFilePath, 'r') as codonFile:
            codonReader = csv.reader(codonFile, delimiter="\t")
            Rna.codonMap = {row[0]:row[1] for row in codonReader if len(row)>0}

    def __init__(self, rna, codonMapFile=None):
        self.rna = rna
        if codonMapFile is not None:
            self.constructCodonMap(codonMapFile)

    def translate(self):
        if Rna.codonMap is None:
            raise RuntimeError('codon map must be generated first.')
        protein =''.join([Rna.codonMap[self.rna[i:i+3]] for i in range(0, len(self.rna),3) if Rna.codonMap[self.rna[i:i+3]]!='stop'])
        return protein


