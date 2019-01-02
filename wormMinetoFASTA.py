#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Function to convert the wormmine protein sequence to  FASTA seqeunces suitable
for input into online alignment tools"""


import sys
import pandas as pd

def wormMinetoFASTA(inputfile, outputfile):
    inDF = pd.read_csv(inputfile, index_col = False)
    with open(outputfile, 'w') as fid:
        for entry in inDF.iterrows():
            fid.write('>')
            fid.write(entry[1]['Gene.CDSs.protein.primaryAccession'])
            fid.write('\n')
            fid.write(entry[1]['Gene.CDSs.protein.sequence.residues'])
            fid.write('\n')                

if __name__ == '__main__':
    inputfile = sys.argv[1]
    outputfile= sys.argv[2]
    wormMinetoFASTA(inputfile, outputfile)