import collections
import re

class Dna:
    def __init__(self, dna):
        self.dna = dna

    def transcribe(self):
        '''transcription from dna to rna'''
        return re.sub('T', 'U', self.dna)

    def countNucleotides(self):
        '''Count nucleotides in dna'''
        counter = collections.Counter(self.dna)
        return {c:counter[c] for c in 'ACGT'}

    def revc(self):
        ''' get reverse complement of dna'''
        cmap = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        rc = [cmap[nt] for nt in self.dna]
        rc.reverse()
        return ''.join(rc)

    def gc(self):
        '''return GC content (percentage of C or G)'''
        freq = self.countNucleotides()
        gccontent = 100.0 * sum([float(freq[nt]) for nt in 'GC']) / sum([freq[nt] for nt in 'ACGT'])
        return gccontent






