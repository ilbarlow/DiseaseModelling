#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Function to convert the wormmine protein sequence to  FASTA seqeunces suitable
for input into online alignment tools"""


import sys
import pandas as pd
import os

def wormMinetoFASTA(inputfile):
    """ function to generate FASTA.txt file that can be used for protein alignments"""

    #read the csv
    inDF = pd.read_csv(inputfile, index_col = False)
    inDF = inDF.drop_duplicates('Gene.CDSs.protein.symbol')

    #generate output file
    outputFile = os.path.join(os.path.dirname(inputfile), 'FASTA_sequences.txt')

    with open(outputFile, 'w') as fid:
        for i,entry in inDF.iterrows():
            fid.write('>')
            fid.write(entry['Gene.CDSs.protein.symbol'])
            fid.write('\n')
            fid.write(entry['Gene.CDSs.protein.sequence.residues'])
            fid.write('\n')


if __name__ == '__main__':
    inputfile = sys.argv[1]
    wormMinetoFASTA(inputfile)