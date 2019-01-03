"""download gDNA sequences for selected genes"""

import sys
import os
import pandas as pd

def gDNA_download(GenomeLocs, output_dir):
	"""function to download the Genome Locations for each of the target genes
	Input -
	GenomeLocs - .csv file from WormMine containing the fields Gene.locations.start 
	and Gene.locations.end and Gene.symbol and Gene.chromosome.primaryIdentifier

	output_dir - directory into which all the gDNA sequences will be saved

	Output -
	gDNA sequences in fasta format"""

	inputDF = pd.read_csv(GenomeLocs)

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	for i,r in inputDF.iterrows():
		start = r['Gene.locations.start']
		end = r['Gene.locations.end']

		cmd = './twoBitToFa http://hgdownload.cse.ucsc.edu/gbdb/ce11/ce11.2bit ' + \
		os.path.join(output_dir, r['Gene.symbol'] + '.fa')
		cmd += ' -seq=chr'+ r['Gene.chromosome.primaryIdentifier']
		cmd += ' -start=' + str(r['Gene.locations.start'])
		cmd += ' -end=' + str(r['Gene.locations.end'])
		print (cmd)

		os.system (cmd)

if __name__ == '__main__':
	GenomeLocs = sys.argv[1]
	output_dir = sys.argv[2]
	gDNA_download(GenomeLocs, output_dir)

