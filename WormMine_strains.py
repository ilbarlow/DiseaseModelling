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

	def findStrains(inputIDs, outputFile):

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
            try:
                _q = query.where('Gene.primaryIdentifier', '=', i).results('dict')
                print ('{} strains found for {}'.format(len(_q), i))
                for row in _q:
                    strainsDF = strainsDF.append(pd.Series(row), ignore_index=True)
            except StopIteration:
                print ('0 strains found for {}'.format(i))
          
	    # strainsDF[strainsDF['Gene.strains.genotype']!='']
	    strainsDF = strainsDF.drop_duplicates(subset = 'Gene.strains.primaryIdentifier')
	    strainsDF.to_csv(outputFile, index=False)
        
        strain_summary = {
                'CGCstrains': list(set(np.unique(strainsDF['Gene.primaryIdentifier'])).intersection(IDs)),
                'noCGCstrains': list(set(np.unique(strainsDF['Gene.primaryIdentifier'])).symmetric_difference(IDs))
                }
	    with open(os.path.join(os.path.dirname(inputIDs), 'strainSummary.csv'),'w') as f:
            w = csv.writer(f)
            w.writerow(strain_summary.keys())
            w.writerows(itertools.zip_longest(*strain_summary.values()))

        return

	findStrains(inputIDs, outputFile)
	    
