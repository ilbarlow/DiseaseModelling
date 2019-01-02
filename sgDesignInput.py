#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 16:17:50 2018

@author: ibarlow
"""

""" function to convert WormMine output to a csv with columns Gene name, 
chromosome, start, and end of targeted region from left to right"""

import pandas as pd
import sys

def sgDesignInput(input_file, output_file):
    """Inputs:
    input_file - csv of from WormMine output with columns called Gene.symbol, Gene.chromosome.primaryIdentifier,
    Gene.locations.start, Gene.locations.end 

    Output:
    output_file - csv of genes with their genomic locations. Headers are name, chrom, start, and end ready for design_guides.py"""
    df = pd.read_csv(input_file, index_col=0)
    
    df2 = pd.DataFrame(df['Gene.symbol'])
    df2['chrom'] = df['Gene.chromosome.primaryIdentifier'].apply(lambda x: 'chr'+x)
    
    df2 = pd.concat([df2, df[['Gene.locations.start', 'Gene.locations.end']].astype(int)],axis=1)
    
    df2 = df2.rename(columns = {'Gene.symbol': 'name', 'Gene.locations.start':'start',\
                          'Gene.locations.end':'end'})
    
    df2.to_csv(output_file, index=False)
    
if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    sgDesignInput(input_file, output_file)
    
    