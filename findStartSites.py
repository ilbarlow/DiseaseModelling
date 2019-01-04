import sys
import pandas as pd
import numpy as np
import os

def findStartSites(CDNAcsv, output_file):
    """function to identify the first exon of the transcript based on strand and also outputs the cDNA sequences into a new folder
    Input - 
    cDNAcsv - csv file output by WormMineCDNAsearch.py
    
    output_file - csv file with name, chrom, start-50, end (for looking for guides) to be saved

    Ouput - 

    output_file - with details added

    output_cDNAs - cDNA sequences saved into a folder called 'cDNA_sequences'
    """

    #check if already outputDir
    outputDir = os.path.join(os.path.dirname(CDNAcsv), 'cDNA_sequences')
    #if not make a new directory
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    #import CDNAseqeunces from CDNAsearch
    inputDF = pd.read_csv(CDNAcsv)

    worms = list(np.unique(inputDF['Protein.CDSs.gene.primaryIdentifier']))

	#groupy the dataframe by wormID so that can loop through to to identify the first and last exons
	#go through finding the first exon and load into a new dataframe
    startExons = pd.DataFrame()
    grouped = inputDF.groupby('Protein.CDSs.gene.primaryIdentifier')
    for i in worms:
        t = grouped.get_group(i)
        temp = pd.Series()
        temp['name'] = t.iloc[0]['Protein.CDSs.gene.name']
        temp['chrom'] = 'chr'+t.iloc[0]['Protein.CDSs.chromosomeLocation.locatedOn.primaryIdentifier']
        #find if first exon is at beginning or end of locations - this depends on the strand
        if sum(t['Protein.CDSs.transcripts.exons.locations.strand'] == 1) == t.shape[0]:
            temp['start'] = int(t['Protein.CDSs.transcripts.exons.locations.start'].min()-50)
            temp['end'] = int(temp['start']+250)
        elif sum(t['Protein.CDSs.transcripts.exons.locations.strand'] == -1) == t.shape[0]:
            temp['start'] = int(t['Protein.CDSs.transcripts.exons.locations.end'].max()+50)
            temp['end'] = int(temp['start']-250)
        startExons = startExons.append(temp.to_frame().transpose())
        
        with open(os.path.join(outputDir, temp['name'] + '_cDNA.ape'), 'w+') as fopen:
            fopen.write(t.iloc[0]['Protein.CDSs.transcripts.sequence.residues'])

        del temp	
        
    startExons.to_csv(output_file, index=False)
        
if __name__ == '__main__':
    CDNAcsv = sys.argv[1]
    output_file = sys.argv[2]
    findStartSites(CDNAcsv, output_file)
    
    