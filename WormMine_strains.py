#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from intermine.webservice import Service
import pandas as pd
import os
import io
import numpy as np

service = Service("http://intermine.wormbase.org/tools/wormmine/service")

if __name__ == '__main__':

	query = service.new_query('Gene')
	query.add_view(
		"strains.CGCReceived",\
		"strains.genotype",\
		"strains.laboratory",\
		"strains.madeBy",\
		"strains.mutagen",\
		"strains.ncbiTaxonomyID",\
		"strains.otherName",\
		"strains.outcrossed",\
		"strains.primaryIdentifier",\
		"strains.species",\
		"strains.remark",\
		"strains.strainHistory"\
		)

	inputIDs = sys.argv[1]
	outputFile = sys.argv[2]

	def findStrains(inputIDS, outputFile):

		#open csv
	    with io.open(inputIDs, 'r',  encoding='utf-8-sig') as fid:
	    	 IDs = fid.read().split(',')
	    	#remove dupliates
	    IDs = np.unique(IDs)
	    print (IDs)

	    strainsDF = pd.DataFrame()
	    for i in IDs:
	    	print ('Finding {}'.format(i))
	    	total_strains = 0
	    	for r in query.rows():
	    		if i in r and r['Gene.strains.genotype']!=None:
	    			total_strains+=1
	    			print ('{} strains found'.format(i))
	    			strainsDF = strainsDF.append(
	    							pd.Series(r.to_d()),\
	    							ignore_index=True
	    							)
	    	print ('{} strains found for {}'.format(total_strains, i))
	    # strainsDF[strainsDF['Gene.strains.genotype']!='']
	    strainsDF = strainsDF.drop_duplicates(subset = 'Gene.strains.primaryIdentifier')
	    strainsDF.to_csv(outputFile, index=False)
	    return

	findStrains(inputIDs, outputFile)
	    
