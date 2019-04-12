#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from intermine.webservice import Service
import pandas as pd
import csv
import os
import io
import numpy as np
import itertools

service = Service("http://intermine.wormbase.org/tools/wormmine/service")

def find_strain(inputID, strainDF):
	query = service.new_query('Gene')
	query.add_view("strains.*")
	_q = query.where('Gene.primaryIdentifier', '=', inputID).results('dict')
	print ('{} strains found for {}'.format(len(_q), inputID))
	for row in _q:
	    strainDF = strainDF.append(pd.Series(row), ignore_index=True)
	return strainDF
    

if __name__ == '__main__':

	inputIDs = sys.argv[1]
	outputFile = sys.argv[2]

	def findStrains(inputIDs, outputFile):

		#open csv
		with io.open(inputIDs, 'r',  encoding='utf-8-sig') as fid:
			IDs = fid.read().split(',')
			#remove dupliates
		IDs = np.unique(IDs)
		print (IDs)
		CGCstrainsDF = pd.DataFrame()
		for i in IDs:
		    print ('Finding {}'.format(i))		    
		    try:
		        CGCstrainsDF = find_strain(i, CGCstrainsDF)

		    except StopIteration:
		        print ('strain error for {}'.format(i))
		        # query = service.new_query('Gene')
		        # query.add_view(
		        # 		"strains.*"
		        # 		)
		        # CGCstrainsDF = find_strain(i, CGCstrainsDF)
	        
	        
		# strainsDF[strainsDF['Gene.strains.genotype']!='']
		CGCstrainsDF = CGCstrainsDF.drop_duplicates(subset = 'Gene.strains.primaryIdentifier')
		CGCstrainsDF.to_csv(outputFile, index=False)

		CGCstrainsDF['gene_summary'] = list(zip(CGCstrainsDF['Gene.primaryIdentifier'], CGCstrainsDF['Gene.symbol']))

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

findStrains(inputIDs, outputFile)