#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" This function is to find the existing strains or sequence information for
the selected wormbase ID"""

import sys
from intermine.webservice import Service
import pandas as pd
import csv
import os
import io
import numpy as np

service = Service("http://intermine.wormbase.org/tools/wormmine/service")

def find_strain(inputID, outputDF):
    """ Accesses the wormmine database and compiles the query and outpu
    depending on whether strains or sequences are wanted
    Inputs:
        inputID = wormbaseID
        outputDF = dataframe to which the results will be attached
    Output:
        outputDF = dataframe containing the wormbase information"""    
    query = service.new_query('Gene')
    query.add_view("strains.*")
    _q = query.where('Gene.primaryIdentifier','=', inputID).results('dict')
    print ('{} strains found for {}'.format(len(_q), inputID))
    for row in _q:
	    outputDF = outputDF.append(pd.Series(row), ignore_index=True)
    return outputDF

def find_sequence(inputID, outputDF):
    query = service.new_query('Gene')
    query.add_view('Gene.CDSs.*',
                   'transcripts.length')
    _q =query.where('Gene.primaryIdentifier', '=', inputID).results('dict')
    for row in _q:
	    outputDF = outputDF.append(pd.Series(row), ignore_index=True)
    return outputDF
           
    
def findStrains(inputIDs, outputFile):
    """ Function imports the wormbase IDs, runs find_strains to run wormmine and
    then output an output file containing the scraped data
    Input:
        inputIDs = .csv containing of the requested WormBaseIDs
        outputFile = filename for the outputfile (will be saved in the same directory
        as input file)
        strains= Bool
        sequence = Bool
    Output:
        if strains = True, CGCstrains.csv and noCGCstrains.csv f
        output file with scraped data"""
        
    #open csv
    with io.open(inputIDs, 'r',  encoding='utf-8-sig') as fid:
        IDs = fid.read().split(',')
        #remove dupliates
        IDs = np.unique(IDs)
        print (IDs)
		    
    CGCstrainsDF = pd.DataFrame()
    for i in IDs:
        print ('Finding {} strains'.format(i))    				    
        try:
            CGCstrainsDF = find_strain(i,
                                       CGCstrainsDF,
                                       )
        except StopIteration:
            ('strain error for {}'.format(i))
    CGCstrainsDF = CGCstrainsDF[CGCstrainsDF['Gene.strains.genotype'].notnull()]
    CGCstrainsDF = CGCstrainsDF.drop_duplicates(subset = 'Gene.strains.primaryIdentifier')
    CGCstrainsDF.to_csv(os.path.join(os.path.dirname(inputIDs),outputFile),
                        index=False)
    CGCstrainsDF['gene_summary'] = list(zip(CGCstrainsDF['Gene.primaryIdentifier'],CGCstrainsDF['Gene.symbol']))
	
    #make a dictionary of WBIDs and gene names
    WBID_dict = {}
    [WBID_dict.update({k:v}) for k,v in np.unique(CGCstrainsDF['gene_summary'])]
    noCGC_strains = list(set(WBID_dict.keys()).symmetric_difference(IDs))
        
    with open(os.path.join(os.path.dirname(inputIDs), 'CGCstrains.csv'),'w') as f:
        w = csv.writer(f)
        w.writerows(WBID_dict.items())		    
    with open(os.path.join(os.path.dirname(inputIDs), 'noCGCstrains.csv'), 'w') as f:
        w = csv.writer(f)
        for i in noCGC_strains:
            w.writerow([i])
    return

def findSequences(inputIDs, outputFile):
    """ Function imports the wormbase IDs, runs find_strains to run wormmine and
    then output an output file containing the scraped data
    Input:
        inputIDs = .csv containing of the requested WormBaseIDs
        outputFile = filename for the outputfile (will be saved in the same directory
        as input file)
        strains= Bool
        sequence = Bool
    Output:
        if strains = True, CGCstrains.csv and noCGCstrains.csv f
        output file with scraped data"""
        
    #open csv
    with io.open(inputIDs, 'r',  encoding='utf-8-sig') as fid:
        IDs = fid.read().split(',')
        #remove dupliates
        IDs = np.unique(IDs)
        print (IDs)
		    
    sequenceDF = pd.DataFrame()
    for i in IDs:
        print ('Finding {} sequences'.format(i))
        try:
            sequenceDF = find_sequence(i,
                                     sequenceDF,
                                     )
        except Exception as error:
            print (error)            
    sequenceDF.to_csv(os.path.join(os.path.dirname(inputIDs),outputFile), index=False)

if __name__ == '__main__':
    inputIDs = sys.argv[1]
    outputFile = sys.argv[2]
    strains = bool(sys.argv[3])
    sequences = bool(sys.argv[4])
 
    if sequences:
        findSequences(inputIDs, outputFile)

    elif strains:
        findStrains(inputIDs, outputFile)
    

