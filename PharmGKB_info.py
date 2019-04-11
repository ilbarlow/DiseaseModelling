#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:56:33 2019

@author: ibarlow
"""

""" Pairing up the Neuro Disease model selected genes with associated known 
pharmacogonomics data from PharmGKB data"""

import pandas as pd
import sys
import os
import glob

selected_genes_in = "/Volumes/behavgenom$/Ida/CRISPRs/DiseaseModels/NeuroDiseaseModelling/v2/RefiningGenes/TopGenes_FINAL.xlsx"
GKBin = "/Users/ibarlow/OneDrive - Imperial College London/Documents/DrugScreening/PharmGKB"
GKB_links = ['genes', 'phenotypes', 'drugs', 'chemicals', 'relationships', 'variants']

#load the PharmGKBdata
GKBfiles = glob.glob(os.path.join(GKBin, '**/*.tsv'))
GKBdata = {}
for file in GKBfiles:
    for link in GKB_links:
        if link in file:
            GKBdata[link] = pd.read_csv(file, sep='\t')
            
#load the selected genes and make a list of them
selected_genesDF = pd.read_excel(selected_genes_in, 'Top32')
EnsemblIDs = []
GeneSymbols = []
for i,r in selected_genesDF.iterrows():
    for j in r.EnsemblID.split(','):
        EnsemblIDs.append(j.lstrip())
    for j in r.HGNCSymbol.split(','):
        GeneSymbols.append(j.lstrip())
        
#find the pharmGKB ids for the selected genes
GKB_geneDF = pd.DataFrame()
for gene in EnsemblIDs:
    GKB_geneDF = GKB_geneDF.append(GKBdata['genes'][GKBdata['genes']['Ensembl Id']==gene], sort=True)
for gene in GeneSymbols:
    GKB_geneDF = GKB_geneDF.append(GKBdata['variants'][GKBdata['variants']['Gene Symbols']==gene], sort=True)

#tie up the relationships with the PharmGKB relationships
GKB_relationships = pd.DataFrame()
for i,r in GKB_geneDF.iterrows():
    GKB_relationships = GKB_relationships.append(
            GKBdata['relationships'][GKBdata['relationships'].apply(lambda row: row.astype(str).str.contains(r['PharmGKB Accession Id']).any(), axis=1)],
            ignore_index=True)
    print ('{} relationships for {} found'.format(
            (GKBdata['relationships'].apply(lambda row: row.astype(str).str.contains(r['PharmGKB Accession Id']).any(), axis=1)).sum(),
            r['PharmGKB Accession Id']))
    
GKB_relationships.to_csv(os.path.join(os.path.dirname(selected_genes_in), 'PharmGKB_Top32_associations.csv'))
