def parse_fasta(fh):
    ''' 
    parse a fasta file.
    return dna strings.
    '''
    lines = fh.readlines()
    dnaStrings = []
    curDna = None
    for l in lines:
        if l.startswith('>'):
            # save previous dna
            if curDna is not None:
                dnaStrings.append(curDna)
            # start a new dna 
            curDna = ''
        else:
            curDna += l.rstrip();
    #
    if curDna is not None:
        dnaStrings.append(curDna)
    return dnaStrings
