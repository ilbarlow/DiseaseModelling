import sys
from intermine.webservice import Service
import pandas as pd
import numpy as np
service = Service("http://intermine.wormbase.org/tools/wormmine/service")

def WormMineCDNAsearch(inputProteins, outputfile):
	# Get a new query on the class (table) you will be querying:
    query = service.new_query("Protein")
    # The view specifies the output columns
    query.add_view(
    	"primaryAccession", "CDSs.transcripts.sequence.residues",
    	"CDSs.transcripts.locations.start", "CDSs.transcripts.locations.end",
    	"CDSs.transcripts.locations.strand", "CDSs.gene.name", "CDSs.gene.primaryIdentifier",
    	"CDSs.gene.locations.start", "CDSs.gene.locations.end", \
    	"CDSs.chromosomeLocation.locatedOn.primaryIdentifier", "CDSs.protein.sequence.residues",\
    	"CDSs.protein.sequence.length", "CDSs.transcripts.exons.locations.end", "CDSs.transcripts.exons.locations.start",\
        "CDSs.transcripts.exons.locations.strand", "CDSs.transcripts.exons.sequence.residues"
    	)

    #open csv
    with open(inputProteins, 'r', encoding='utf-8-sig') as fid:
    	 homologues = fid.read().split(',')

    print (homologues)

    #look for CEID in query and add to a dataframe
    hom_DF = pd.DataFrame()
    for item in homologues:
        for row in query.rows():
            if item in row:
                print (row["primaryAccession"] + '+' + row["CDSs.chromosomeLocation.locatedOn.primaryIdentifier"])
                if pd.Series(row.to_d()).isna().sum()>0:
                    temp = pd.Series(row.to_d()).fillna(np.nan)
                    hom_DF = hom_DF.append(temp, ignore_index=True)
                else:
                    hom_DF = hom_DF.append(pd.Series(row.to_d()), ignore_index=True)

    #save output to csv
    hom_DF.to_csv(outputfile, index=False)

if __name__ == '__main__':
	inputProteins = sys.argv[1]
	outputfile = sys.argv[2]
	WormMineCDNAsearch(inputProteins, outputfile)
    

