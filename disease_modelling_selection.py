#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 15:32:35 2018

@author: ibarlow
"""

""" Refine gene selection for disease modelling - import csvs from matlab and 
combine indications """

import pandas as pd
import os
import glob
import numpy as np

#import excel file
fid = '/Volumes/behavgenom$/Ida/CRISPRs/DiseaseModels/NeuroDiseaseModelling/v2/GeneSelection/'

import_list = glob.glob(os.path.join(fid, '*.csv'))
dataRaw = pd.DataFrame()
for file in import_list:
    dataRaw = dataRaw.append(pd.read_csv(file, index_col =False), sort=True)

#find the unique human gene identifiers - for a single worm gene there are several human orthologs
deduplicatedData = dataRaw.drop_duplicates(subset = ['WormBaseID', 'EnsemblID'])
uniqueWormIDs = list(deduplicatedData.WormBaseID.unique())
uniqueHumanIDs = list(deduplicatedData.EnsemblID.unique())
#go through the wormIDs and concatenate OMIM phenotypes and EnsemblIDs
wormCompiled = pd.DataFrame()
for ID in uniqueWormIDs:
    _temp = deduplicatedData.groupby(by='WormBaseID').get_group(ID).drop_duplicates('WormBaseID')
    _temp.OMIMPhenotypes = np.unique(', '.join(deduplicatedData.groupby(by='WormBaseID').get_group(ID).OMIMPhenotypes))
    _temp.EnsemblID = ', '.join(deduplicatedData.groupby(by='WormBaseID').get_group(ID).EnsemblID)
    _temp.HGNCSymbol = ', '.join(deduplicatedData.groupby(by='WormBaseID').get_group(ID).HGNCSymbol)
    _temp.LocusID = ', '.join(deduplicatedData.groupby(by='WormBaseID').get_group(ID).LocusID)
    wormCompiled = wormCompiled.append(_temp)
    del _temp
wormCompiled =wormCompiled.reset_index(drop=True)
wormCompiled.to_csv(os.path.join(os.path.dirname(fid), 'NeuroGeneSelection.csv'))

#%% Same for the myosin models

#import excel file
fid = '/Users/ibarlow/Documents/MATLAB/WormDiseaseModelling/Selected_GeneLists/wormMyosinHomologue.xls'
dataRaw = pd.read_excel(fid, sheet_name=None)
dataConcat = pd.concat(dataRaw.values(), ignore_index=True, sort=True)
dataConcat = dataConcat.fillna('')

#find the unique human gene identifiers - for a single worm gene there are several human orthologs
#wormData = dataConcat.drop_duplicates(subset = 'WormBaseID')
uniqueWIDs = list(dataConcat.WormBaseID.unique())
#combine OMIM phenotypes
wormCompiled = pd.DataFrame()
for ID in uniqueWIDs:
    temp = dataConcat.groupby(by='WormBaseID').get_group(ID).drop_duplicates(subset = 'WormBaseID')
    temp.OMIMPhenotypes = ','.join(dataConcat.groupby(by='WormBaseID').get_group(ID).OMIMPhenotypes)
    temp.EnsemblID = ', '.join(dataConcat.groupby(by='WormBaseID').get_group(ID).EnsemblID)
    temp.HGNCSymbol = ', '.join(dataConcat.groupby(by='WormBaseID').get_group(ID).HGNCSymbol)
    wormCompiled = wormCompiled.append(temp)
    del temp
wormCompiled = wormCompiled.reset_index(drop=True)

#save
wormCompiled.to_csv(os.path.join(os.path.dirname(fid), 'wormMyosinHomologues_updated.csv'))
