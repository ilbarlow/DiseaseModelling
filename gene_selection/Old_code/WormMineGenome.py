#!/usr/bin/env python3

"""To mine the WormBase and get all the strains for a WormBaseID"""
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/


# The following four lines will be needed in every python script:
import sys
from intermine.webservice import Service
import pandas as pd
service = Service("http://intermine.wormbase.org/tools/wormmine/service")
#WormID = '/Users/ibarlow/Documents/MATLAB/WormDiseaseModelling/Selected_GeneLists/WormIDs.csv'
# '~/Documents/DiseaseModel_wormStrains.csv'

def WormMineGenome(file_in, output_file):
    # Get a new query on the class (table) you will be querying:
    query = service.new_query("Gene")

    # The view specifies the output columns
    query.add_view(
        "name", "primaryIdentifier", "locations.start", "locations.end", \
        "locations.strand", "chromosome.primaryIdentifier"
        )

    with open(file_in, 'r', encoding='utf-8-sig') as fid:
    	worms = fid.read().split(',')

    print (worms)

    #find the associated genomic locations of the wormbaseIDs
    expStrains = []
    for item in worms:
        if item =='':
            continue
        else:
            for row in query.rows():
                if item in row:
                    print(row['name'], row['primaryIdentifier'], \
                        row['locations.start'], row['chromosome.primaryIdentifier'])
                    expStrains.append(row) 

    df = pd.DataFrame()
    for i in expStrains:
    	df = df.append(pd.Series(i.to_d()), ignore_index=True)

    df.to_csv(output_file)

if __name__ == '__main__':
    file_in = sys.argv[1]
    output_file = sys.argv[2]
    WormMineGenome(file_in, output_file)