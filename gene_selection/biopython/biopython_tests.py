#!/usr/bin/env/biopython python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 10:31:13 2020

@author: ibarlow
"""

import os
import sys

#mysql attempts
os.chdir('/Users/ibarlow/repositories/chopchop/')

name = 'FAM78B'

# mysql -h genome-mysql.soe.ucsc.edu -ugenome -A -e "select \
#    e.chrom, e.txStart, e.txEnd, e.strand, e.name, j.name as geneSymbol from ncbiRefSeqCurated e,\
#    ncbiRefSeqLink j where e.name = j.id AND e.chrom='${chrom}' AND \
#       ((e.txStart >= ${chromStart} - 10000 AND e.txStart <= ${chromEnd} + 10000) OR \ (e.txEnd >= ${chromStart} - 10000 AND e.txEnd <= ${chromEnd} + 10000)) \


cmd = 'mysql --host=genome-mysql.soe.ucsc.edu --user=genomep --password=password -A -P 3306 -e \"select l.name, l.name2, l.txStart, l.txEnd from ncbiRefSeq l where l.name2=\'{}\'" hg38'.format(name)
os.system(cmd)

# coudln't figure out how to save sql output to a file
# sytax :  INTO OUTFILE '/tmp/orders.txt'


#biopython version
SeqIO.parse("FAM78B", "genbank")


request = Entrez.epost("gene",id='FAM89B')
result = Entrez.read(request)

from Bio import SeqIO
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "i.barlow@lms.mrc.ac.uk"  # Always tell NCBI who you are
handle = Entrez.efetch(db="nucleotide", id=IDArray[0], rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
>>> print(record.id)
EU490707.1
>>> print(record.name)
EU490707
>>> print(record.description)
Selenipedium aequinoctiale maturase K (matK) gene, partial cds; chloroplast
>>> print(len(record.features))
3
>>> print(repr(record.seq))
Seq('ATTTTTTACGAACCTGTGGAAATTTTTGGTTATGACAATAAATCTAGTTTAGTA...GAA', IUPACAmbiguousDNA())


#find gene entrez id
sterm = 'FAM78B'+ '[sym] \"Homo Sapiens\"[orgn]'
handle = Entrez.esearch(db="gene", retmode = "xml", term = sterm )
record = Entrez.read(handle)
IDArray = record["IdList"]
