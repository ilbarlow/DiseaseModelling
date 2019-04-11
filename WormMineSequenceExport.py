#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from intermine.webservice import Service
import pandas as pd
import numpy as np
import os
import io
service = Service("http://intermine.wormbase.org/tools/wormmine/service")

#sys.path.insert(0, '/Users/ibarlow/Documents/GitHub/pythonScripts/WormMine')

def wormMinetoFASTA(inputfile):
    """ function to generate FASTA.txt file that can be used for protein alignments"""

    #read the csv
    inDF = pd.read_csv(inputfile, index_col = False)
    inDF = inDF.drop_duplicates('Gene.CDSs.protein.symbol')

    #generate output file
    outputFile = os.path.join(os.path.dirname(inputfile), 'FASTA_sequences.txt')

    with open(outputFile, 'w') as fid:
        for i,entry in inDF.iterrows():
            fid.write('>')
            fid.write(entry['Gene.CDSs.protein.symbol'])
            fid.write('\n')
            fid.write(entry['Gene.CDSs.protein.sequence.residues'])
            fid.write('\n')


def WormMineSequenceExport(inputIDs, margin, outputDir, createFASTA):
    """ Function to find cDNAs from list of WormBaseIDs or protein IDs queries (as csv) and outputs
    cDNA sequences, gDNA sequence, sgRNAsearch.csv (name, chrom, start, end) for chopchop query, and .csv
    with all the WormMine output information
    Input:
    inputIDs - .csv containing list of wormbase or protein ids

    margin - two-element tuple of bp around the gene start to look for sgRNAs - for use in chopchop
    eg. [50 250] will put 'start' locations 50bp before 1st exon start and 'end' 250bp after 1st exon start

    outputDir - name of directory into which the output .csvs and sequences should be saved

    createFASTA - boolean. If True then a FASTA.txt file is generated for protein alignments

    Output:
    .ape cDNA sequences
    
    sgInputs.csv - for use in chopchop

    WormMineOutput.csv - query results

    FASTA.txt - text file of protein sequences
    """

    # Get a new query on the class (table) you will be querying:
    query = service.new_query('Gene')

    # The view specifies the output columns
    query.add_view(
        "name", "primaryIdentifier", "locations.start", "locations.end", \
        "locations.strand", "chromosome.primaryIdentifier", "CDSs.transcripts.sequence.residues",\
        "CDSs.transcripts.locations.start", "CDSs.transcripts.locations.end",\
        "CDSs.transcripts.locations.strand","CDSs.transcripts.exons.locations.end",\
        "CDSs.transcripts.exons.locations.start", "CDSs.transcripts.exons.locations.strand",\
        "CDSs.protein.sequence.residues",\
        "CDSs.protein.sequence.length", "CDSs.protein.symbol")

    #open csv
    with io.open(inputIDs, 'r',  encoding='utf-8-sig') as fid:
    	 IDs = fid.read().split(',')
    #remove dupliates
    IDs = np.unique(IDs)
    
    print (IDs)

    #look for CEID or WormBaseID in query and add to a dataframe
    # loop over rows instead first or find more efficient way of finding wormbase
    # IDs or extract positions of IDs and load positions
    # is there a way of storing the database locally
    input_DF = pd.DataFrame()
    for item in IDs:
        print ('Finding ' + item)
        for row in query.rows():
            if item in row:
                print (row['primaryIdentifier'] + '+' + row["chromosome.primaryIdentifier"])
                if pd.Series(row.to_d()).isna().sum()>0:
                    temp = pd.Series(row.to_d()).fillna(np.nan)
                    input_DF = input_DF.append(temp, ignore_index=True)
                else:
                    input_DF = input_DF.append(pd.Series(row.to_d()), ignore_index=True)
            else:
                continue
 
    #to save outputs need to first check if outputdirectory exits
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        print (outputDir + ' created')
    
    #make a subdirectory for cDNA sequences
    cDNAdir = os.path.join(outputDir, 'cDNAsequences')
    if not os.path.exists(cDNAdir):
        os.makedirs(cDNAdir)

    #group dataframe by wormbaseID and to identify first and last exons and export cDNA sequence as .ape
    startExons = pd.DataFrame()
    grouped=  input_DF.groupby('Gene.primaryIdentifier')
    for i in IDs:
        t = grouped.get_group(i)
        temp = pd.Series()
        temp['name'] = t.iloc[0]['Gene.name']
        temp['chrom'] = 'chr'+t.iloc[0]['Gene.chromosome.primaryIdentifier']     
        #find if first exon is at beginning or end of locations - this depends on the strand
        if (t['Gene.CDSs.transcripts.exons.locations.strand'] == '1').sum() == t.shape[0]:
            temp['start'] = int(t['Gene.CDSs.transcripts.exons.locations.start'].min()-margin[0])
            temp['end'] = int(temp['start']+margin[1])
        elif (t['Gene.CDSs.transcripts.exons.locations.strand'] == '-1').sum() == t.shape[0]:
            temp['start'] = int(t['Gene.CDSs.transcripts.exons.locations.end'].max()+margin[0])
            temp['end'] = int(temp['start']-margin[1])
        startExons = startExons.append(temp.to_frame().transpose())
        
        #save cDNA sequence
        with open(os.path.join(cDNAdir, temp['name'] + '_cDNA.ape'), 'w+') as fopen:
            fopen.write(t.iloc[0]['Gene.CDSs.transcripts.sequence.residues'])
            print (temp['name'] + ' cDNA sequence saved')

    #save sgRNA input csv from startExons
    startExons.to_csv(os.path.join(outputDir, 'sgInputs.csv'), index=False)

    #now run twoBittoFa to download the gDNAsequences
    #create subdirectory
    gDNAdir = os.path.join(outputDir, 'gDNAsequences')
    if not os.path.exists(gDNAdir):
        os.makedirs(gDNAdir)
 
    # save WormMine output to csv
    input_DF.to_csv(os.path.join(outputDir, 'WormMineoutput.csv'), index=False)

    #drop duplicates for saving gDNA
    input_DF = input_DF.drop_duplicates('Gene.symbol')

    for i,r in input_DF.iterrows():
        start = r['Gene.locations.start']
        end = r['Gene.locations.end']
        cmd = './twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/ce11/ce11.2bit ' + \
        os.path.join(gDNAdir, r['Gene.name'] + '_gDNA.fa')
        cmd += ' -seq=chr'+ r['Gene.chromosome.primaryIdentifier']
        cmd += ' -start=' + str(int(start))
        cmd += ' -end=' + str(int(end))
        print (cmd)

        os.system (cmd)

    if createFASTA == True:
        wormMinetoFASTA(os.path.join(outputDir, 'WormMineoutput.csv'))
    else:
        print('no FASTA.txt file generated')


if __name__ == '__main__':
    inputIDs = sys.argv[1]
    margin = (sys.argv[2], sys.argv[3])
    outputDir = sys.argv[4]
    WormMineSequenceExport(inputIDs, margin, outputDir)
    

