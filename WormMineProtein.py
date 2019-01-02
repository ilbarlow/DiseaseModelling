import sys
from intermine.webservice import Service
import pandas as pd
service = Service("http://intermine.wormbase.org/tools/wormmine/service")
#WormID = '/Users/ibarlow/Documents/MATLAB/WormDiseaseModelling/Selected_GeneLists/WormIDs.csv'
# '~/Documents/DiseaseModel_wormStrains.csv'

def WormMineProtein(file_in, output_file):
    # Get a new query on the class (table) you will be querying:
    query = service.new_query("Gene")

    # The view specifies the output columns
    query.add_view(
         "CDSs.protein.primaryAccession", "CDSs.protein.sequence.residues",\
         "CDSs.protein.sequence.md5checksum", "CDSs.protein.sequence.length"
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
                    print(row["CDSs.protein.primaryAccession"])
                    expStrains.append(row) 

    df = pd.DataFrame()
    for i in expStrains:
    	df = df.append(pd.Series(i.to_d()), ignore_index=True)

    df = df.drop_duplicates(subset = 'Gene.CDSs.protein.primaryAccession', \
        keep = 'first')
    df.to_csv(output_file)

if __name__ == '__main__':
    file_in = sys.argv[1]
    output_file = sys.argv[2]
    WormMineProtein(file_in, output_file)