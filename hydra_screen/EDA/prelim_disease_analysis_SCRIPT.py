#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:55:30 2020

@author: ibarlow

Script for having a first look at the disease models
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from tierpsytools.read_data.hydra_metadata import read_hydra_metadata, align_bluelight_conditions
from sklearn.decomposition import PCA

# features set types
feat256_fname = Path('/Users/ibarlow/tierpsy-tools-python/tierpsytools/extras/feat_sets/tierpsy_256.csv')
featsets = {'prestim': [],
            'bluelight': [],
            'poststim': []}

CONTROL_STRAIN = 'N2'

BAD_FEAT_THRESH = 3 # 3 standard deviations away from the mean
BAD_FEAT_FILTER = 0.1 # 10% threshold for removing bad features
BAD_WELL_FILTER = 0.3 # 30% threshold for bad well
# NAN_CUTOFF = 7000

FEAT_FILE = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/Results/20200626/features_summary_tierpsy_plate_20200720_113846.csv')
FNAME_FILE = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/Results/20200626/filenames_summary_tierpsy_plate_20200720_113846.csv')

METADATA_FILE = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/AuxiliaryFiles/metadata.csv')

root = Path('.').resolve() / 'Desktop/Disease_prelim'
root.mkdir(exist_ok=True)

# %%
if __name__ == '__main__':

    feat = pd.read_csv(FEAT_FILE,
                       comment='#')
    fname = pd.read_csv(FNAME_FILE,
                        comment='#')
    
    meta = pd.read_csv(METADATA_FILE, index_col=None)
    meta = meta[meta.worm_strain.notna()]
    meta.loc[meta.worm_strain == 'N2', 'worm_gene'] = 'N2'
    
    feat256 = []
    with open(feat256_fname, 'r') as fid:
        for r in fid.readlines():
            feat256.append(r.rstrip())
            
    # %% Match metaddata and features
    feat, meta = read_hydra_metadata(feat,
                                     fname,
                                     meta)
    feat, meta = align_bluelight_conditions(feat,
                                            meta)
    
    # Remove wells missing bluelight conditions
    imgst_cols = [col for col in meta.columns if 'imgstore_name' in col]
    miss = meta[imgst_cols].isna().any(axis=1)
    
    feat = feat.loc[~miss,:]
    meta = meta.loc[~miss,:]
    
    # remove features with too many nans
    feat = feat.loc[:, feat.isna().sum(axis=0)/feat.shape[0]<BAD_FEAT_FILTER]
    
    # remove features with std=0
    feat = feat.loc[:, feat.std(axis=0) != 0]
    
    # remove wells with too many nans
    feat = feat.loc[feat.isna().sum(axis=1)/feat.shape[0]<BAD_WELL_FILTER, :]
    
    # fillnans
    feat.fillna(feat.mean(axis=0), inplace=True)
    
    featZ = pd.DataFrame(data=stats.zscore(feat, axis=0),
                         columns=feat.columns,
                         index=feat.index)
    
    # feature sets
    pathcurvature_feats = [x for x in feat.columns if 'path_curvature' in x]
    for f in featsets.keys():
        featsets[f] = [x for x in feat.columns if f in x]
        featsets[f] = list(set(featsets[f]) - set(pathcurvature_feats))
       
    # %% PCA
    pca = PCA()
    
    for f in featsets.keys():
        X2 = pca.fit_transform(featZ[featsets[f]].values)
        cumvar = np.cumsum(pca.explained_variance_ratio_)
        thresh = cumvar <= 0.95 #set 95% variance threshold
        cut_off = int(np.argwhere(thresh)[-1])
    
        #make a plot
        sns.set_style('whitegrid')
        plt.figure()
        plt.plot(range(0, len(cumvar)), cumvar*100)
        plt.plot([cut_off, cut_off], [0, 100], 'k')
        plt.xlabel('Number of Principal Components', fontsize =16)
        plt.ylabel('variance explained', fontsize =16)
        #plt.savefig(os.path.join(directoryA[:-7], 'Figures', 'agarPCvar.png'), dpi =150)
        #plt.savefig(os.path.join(directoryA[:-7], 'Figures', 'agarPCvar.svg'),dpi = 150)
        
        #now put the 1:cut_off PCs into a dataframe
        PCname = ['PC_%d' %(p+1) for p in range (0,cut_off+1)]
        PC_df = pd.DataFrame(data = X2[:,:cut_off+1],
                             columns = PCname,
                             index=featZ.index)
        
        PC_plotting = pd.concat([PC_df,
                                 meta],
                                 axis=1)
        
        plt.figure()
        sns.scatterplot(x='PC_1',
                        y='PC_2',
                        data=PC_plotting,
                        hue='worm_gene',
                        palette='Paired',
                        linewidth=0,
                        alpha=0.8)
        plt.axis('equal')
        plt.xlim([-100, 100])
        plt.ylim([-100, 100])
        plt.legend(loc='center left', 
                        bbox_to_anchor=(1.0,0.5),
                        ncol = 1,
                        frameon= True)
        plt.tight_layout()
        plt.xlabel('PC_1 ({}%)'.format(np.round(cumvar[0]*100,2)))
        plt.ylabel('PC_2 ({})%'.format(np.round((cumvar[1]-cumvar[0])*100, 2)))
        plt.title (f)
        plt.savefig(root / 'PC1_PC2_{}.png'.format(f), dpi=300)
        
    # %% entire feature space
    pca=PCA()
    
    X2 = pca.fit_transform(featZ.values)
    
    
    
    
        
            
            
            
