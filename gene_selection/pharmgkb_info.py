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
import numpy as np

selected_genes_in = "/Volumes/behavgenom$/Ida/CRISPRs/DiseaseModels/NeuroDiseaseModelling/v2/RefiningGenes/Round2/TopGenes_FINAL.xlsx"
GKBin = "/Users/ibarlow/OneDrive - Imperial College London/Documents/DrugScreening/PharmGKB"
GKB_links = ['genes', 'phenotypes', 'drugs', 'chemicals', 'relationships', 'variants']
sunyStrains = "/Volumes/behavgenom$/Ida/CRISPRs/DiseaseModels/NeuroDiseaseModelling/v2/RefiningGenes/SunyStrains_WBIDs.csv"


#load the PharmGKBdata
GKBfiles = glob.glob(os.path.join(GKBin, '**/*.tsv'))
GKBdata = {}
for file in GKBfiles:
    for link in GKB_links:
        if link in file:
            GKBdata[link] = pd.read_csv(file, sep='\t')
            
#load the selected genes and make a list of them and also compile into a dictionary
selected_genesDF = pd.read_excel(selected_genes_in, 'Top32').reset_index(drop=True)
EnsemblIDs = []
GeneSymbols = []
wormOrthologs = {}
for i,r in selected_genesDF.iterrows():
    wormOrthologs[r.CommonName]=[]
    for j in r.EnsemblID.split(','):
        EnsemblIDs.append(j.lstrip())
    for j in r.HGNCSymbol.split(','):
        GeneSymbols.append(j.lstrip())
        wormOrthologs[r.CommonName].append(j.lstrip())
        

#find the pharmGKB ids for the selected genes
GKB_geneDF = pd.DataFrame()
for gene in EnsemblIDs:
    GKB_geneDF = GKB_geneDF.append(GKBdata['genes'][GKBdata['genes']['Ensembl Id']==gene], sort=True)
for gene in GeneSymbols:
    GKB_geneDF = GKB_geneDF.append(GKBdata['variants'][GKBdata['variants']['Gene Symbols']==gene], sort=True)
GKB_geneDF.drop_duplicates(inplace=True).reset_index(drop=True)
GKB_geneDF.reset_index(drop=True, inplace=True)

#tie up the relationships with the PharmGKB relationships
GKB_relationships = pd.DataFrame()
no_relationships = []
for i,r in GKB_geneDF.iterrows():
    _temp = (GKBdata['relationships'][GKBdata['relationships'].apply(lambda row: (row == r['PharmGKB Accession Id']).any(), axis=1)]).copy()
    if type(r['Gene Symbols']) == str:
        _temp['Gene'] = r['Gene Symbols']
    elif type(r['Symbol']) == str:
        _temp['Gene'] = r['Symbol']
#    else:
#        _temp['Gene'] = np.nan
    GKB_relationships = GKB_relationships.append(
            _temp,
            ignore_index=True)
    try:
        print ('{} relationships for {} found'.format(_temp.shape[0], _temp.Gene.unique()[0]))
    except Exception:
        if type(r['Gene Symbols']) == str:
            print ('No relationships found for {}'.format(r['Gene Symbols']))
        if type(r['Symbol']) == str:
            print ('No relationships found for {}'.format(r.Symbol))
        no_relationships.append(r.Symbol)
    del _temp

print ('{} relationships found in total'.format(GKB_relationships.shape[0]))
    
GKB_relationships.to_csv(os.path.join(os.path.dirname(selected_genes_in), 'PharmGKB_Top32_associations.csv'))

#find all the chemical associations
GKB_drug_relationships = GKB_relationships[GKB_relationships.Entity1_type.str.contains('Chemical') | GKB_relationships.Entity2_type.str.contains('Chemical')]

GKB_drug_relationships.to_csv(os.path.join(os.path.dirname(selected_genes_in), 'PharmGKB_Top32_associations_Drugs.csv'), index_label = False)

#list the genes with associations
associatedGenes = list(GKB_drug_relationships.Gene.unique())
GKB_geneGroups = GKB_drug_relationships.groupby('Gene')
