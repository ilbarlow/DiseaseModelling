import sys
from intermine.webservice import Service
import pandas as pd
import numpy as np
service = Service("http://intermine.wormbase.org/tools/wormmine/service")

def WormMineCDNAsearch(inputIDs, margin, outputDir):
	""" Function to find cDNAs from list of WormBaseIDs or protein IDs queries (as csv) and outputs
    cDNA sequences, gDNA sequence, sgRNAsearch.csv (name, chrom, start, end) for chopchop query, and .csv
    with all the WormMine output information
    Input:
    inputIDs - .csv containing list of wormbase or protein ids

    margin - two-element vector of bp around the gene start to look for sgRNAs - for use in chopchop
    eg. [50 250] will put 'start' locations 50bp before 1st exon start and 'end' 250bp after 1st exon start

    outputDir - name of directory into which the output .csvs and sequences should be saved

    Output:
    .ape cDNA sequences
    
    sgInputs.csv - for use in chopchop

    WormMineOutput.csv - query results
    """

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
    with open(inputIDs, 'r', encoding='utf-8-sig') as fid:
    	 IDs = fid.read().split(',')
    #remove dupliates
    IDs = np.unique(IDs)
    
    print (IDs)

    #look for CEID or WormBaseID in query and add to a dataframe
    input_DF = pd.DataFrame()
    for item in homologues:
        for row in query.rows():
            if item in row:
                print (row["primaryAccession"] + '+' + row["CDSs.chromosomeLocation.locatedOn.primaryIdentifier"])
                if pd.Series(row.to_d()).isna().sum()>0:
                    temp = pd.Series(row.to_d()).fillna(np.nan)
                    input_DF = input_DF.append(temp, ignore_index=True)
                else:
                    input_DF = input_DF.append(pd.Series(row.to_d()), ignore_index=True)

    #to save outputs need to first check if outputdirectory exits
    if not os.path.exist(outputDir):
        os.makedirs(outputDir)
        print (outputDir ' created')
    
    #make a subdirectory for cDNA sequences
    cDNAdir = os.path.join(outputDir, 'cDNAsequences')
    os.makedirs(cDNAdir)

    #group dataframe by wormbaseID and to identify first and last exons and export cDNA sequence as .ape
    startExons = pd.DataFrame()
    grouped=  input_DF.groupby('Protein.CDSs.gene.primaryIdentifier')
    for i in IDs:
        t = grouped.get_group(i)
        temp = pd.Series()
        temp['name'] = t.iloc[0]['Protein.CDSs.gene.name']
        temp['chrom'] = 'chr'+t.iloc[0]['Protein.CDSs.chromosomeLocation.locatedOn.primaryIdentifier']
        
        #find if first exon is at beginning or end of locations - this depends on the strand
        if sum(t['Protein.CDSs.transcripts.exons.locations.strand'] == 1) == t.shape[0]:
            temp['start'] = int(t['Protein.CDSs.transcripts.exons.locations.start'].min()-margin[0])
            temp['end'] = int(temp['start']+margin[1])
        elif sum(t['Protein.CDSs.transcripts.exons.locations.strand'] == -1) == t.shape[0]:
            temp['start'] = int(t['Protein.CDSs.transcripts.exons.locations.end'].max()+margin[0])
            temp['end'] = int(temp['start']-margin[1])
        startExons = startExons.append(temp.to_frame().transpose())
        
        #save cDNA sequence
        with open(os.path.join(cDNAdir, temp['name'] + '_cDNA.ape'), 'w+') as fopen:
            fopen.write(t.iloc[0]['Protein.CDSs.transcripts.sequence.residues'])

    #save sgRNA input csv from startExons
    startExons.to_csv(os.path.join(outputDir, 'sgInputs.csv'), index=False)
    
    #save WormMine output to csv
    input_DF.to_csv(os.path.join(outputDir, 'WormMineoutput.csv'), index=False)

if __name__ == '__main__':
	inputIDs = sys.argv[1]
	margin = sys.argv[2]
    outputDir = sys.argv[3]
	WormMineCDNAsearch(inputIDs, margin, outputDir)
    

