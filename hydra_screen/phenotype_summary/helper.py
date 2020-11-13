#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 15:31:14 2020

@author: ibarlow

Helper functions for reading the disease data and making strain specific 
plots
"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from tierpsytools.read_data.hydra_metadata import read_hydra_metadata, align_bluelight_conditions

CONTROL_STRAIN = 'N2'
DATES_TO_DROP = '20200626'
BAD_FEAT_THRESH = 3 # 3 standard deviations away from the mean
BAD_FEAT_FILTER = 0.1 # 10% threshold for removing bad features
BAD_WELL_FILTER = 0.3 # 30% threshold for bad well


STIMULI_ORDER = {'prestim':1,
                 'bluelight':2,
                 'poststim':3}

def drop_nan_worms(feat, meta, saveto, export_nan_worms=False):

    # remove (and check) nan worms
    nan_worms = meta[meta.worm_gene.isna()][['featuresN_filename',
                                             'well_name',
                                             'imaging_plate_id',
                                             'instrument_name',
                                             'imaging_date_yyyymmdd']]
    if export_nan_worms:
        nan_worms.to_csv(saveto / 'nan_worms.csv',
                          index=False)
    feat = feat.drop(index=nan_worms.index)
    meta = meta.drop(index=nan_worms.index)

    return feat, meta

def read_disease_data(feat_file, fname_file, metadata_file, drop_nans=True, export_nan_worms=False):
    """
    
    Parameters
    ----------
    feat_file : TYPE
        DESCRIPTION.
    fname_file : TYPE
        DESCRIPTION.
    metadata_file : TYPE
        DESCRIPTION.
    export_nan_worms : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    feat : TYPE
        DESCRIPTION.
    meta : TYPE
        DESCRIPTION.

    """
    feat = pd.read_csv(feat_file,
                       comment='#')
    fname = pd.read_csv(fname_file,
                        comment='#')
    meta = pd.read_csv(metadata_file, index_col=None)
    assert meta.worm_strain.unique().shape[0] == meta.worm_gene.unique().shape[0]  
    
    feat, meta = read_hydra_metadata(feat,
                                     fname,
                                     meta)
    feat, meta = align_bluelight_conditions(feat,
                                            meta,
                                            how='inner') #removes wells that don't have all 3 conditions
    if drop_nans:
        drop_nan_worms(feat, meta, saveto=feat_file.parent)
    
    return feat, meta


def select_strains(candidate_gene, control_strain, feat_df, meta_df):
    """

    Parameters
    ----------
    strains : TYPE
        DESCRIPTION.
    feat_df : TYPE
        DESCRIPTION.
    meta_df : TYPE
        DESCRIPTION.
    control_strain : TYPE
        DESCRIPTION.

    Returns
    -------
    feat_df : TYPE
        DESCRIPTION.
    meta_df : TYPE
        DESCRIPTION.

    """
    
    gene_list = [g for g in meta_df.worm_gene.unique() if g != control_strain]
    gene_list = [g for g in gene_list if 'myo' not in g and 'unc-54' not in g]
       
    idx = [c for c,g in list(enumerate(gene_list)) if  g==candidate_gene]
    locs = list(meta_df.query('@candidate_gene in worm_gene').index)
    N2_locs = meta_df.query('@control_strain in worm_gene').sample(len(locs)).index
    locs.extend(N2_locs)

     #Only do analysis on the disease strains
    meta_df = meta_df.loc[locs,:]
    feat_df = feat_df.loc[locs,:]
    
    return feat_df, meta_df, idx, gene_list


def filter_features(feat_df, meta_df):
    """

    Parameters
    ----------
    feat_df : TYPE
        DESCRIPTION.
    meta_df : TYPE
        DESCRIPTION.

    Returns
    -------
    feat_df : TYPE
        DESCRIPTION.
    meta_df : TYPE
        DESCRIPTION.
    featlist : TYPE
        DESCRIPTION.
    featsets : TYPE
        DESCRIPTION.

    """
    imgst_cols = [col for col in meta_df.columns if 'imgstore_name' in col]
    miss = meta_df[imgst_cols].isna().any(axis=1)

    # remove data from dates to exclude
    bad_date = meta_df.date_yyyymmdd == float(DATES_TO_DROP)
    # bad wells
    good_wells_from_gui = meta_df.is_bad_well == False
    feat_df = feat_df.loc[good_wells_from_gui & ~bad_date & ~miss,:]
    meta_df = meta_df.loc[good_wells_from_gui & ~bad_date & ~miss,:]
    # remove features with too many nans
    feat_df = feat_df.loc[:,
                          feat_df.isna().sum(axis=0)/feat_df.shape[0]<BAD_FEAT_FILTER]
    # remove wells with too many nans
    feat_df = feat_df.loc[feat_df.isna().sum(axis=1)/feat_df.shape[0]<BAD_WELL_FILTER, :]
    meta_df = meta_df.loc[feat_df.index,:]
    # remove features with std=0
    feat_df = feat_df.loc[:,
                          feat_df.columns[feat_df.std(axis=0)!=0]]

     # feature sets
     # abs features no longer in tierpsy
    pathcurvature_feats = [x for x in feat_df.columns if 'path_curvature' in x]
    
    #remove these features
    feat_df = feat_df.drop(columns=pathcurvature_feats)
    
    featlist = list(feat_df.columns)
    # for f in featsets.keys():
    #     featsets[f] = [x for x in feat.columns if f in x]
    #     featsets[f] = list(set(featsets[f]) - set(pathcurvature_feats))
    
    featsets={}
    for stim in STIMULI_ORDER.keys():
        featsets[stim] = [f for f in featlist if stim in f]
    featsets['all'] = featlist
        
    return feat_df, meta_df, featsets
    

def make_colormaps(gene_list, idx, candidate_gene, CONTROL_STRAIN, featlist):
    """

    Parameters
    ----------
    gene_list : TYPE
        DESCRIPTION.
    idx : TYPE
        DESCRIPTION.
    candidate_gene : TYPE
        DESCRIPTION.
    CONTROL_STRAIN : TYPE
        DESCRIPTION.
    STIMULI_ORDER : TYPE
        DESCRIPTION.
    featlist : TYPE
        DESCRIPTION.

    Returns
    -------
    strain_cmap : TYPE
        DESCRIPTION.
    strain_lut : TYPE
        DESCRIPTION.
    stim_cmap : TYPE
        DESCRIPTION.
    feat_lut : TYPE
        DESCRIPTION.

    """
    
    cmap = list(np.flip((sns.color_palette('cubehelix',
                                           len(gene_list)*2+6))[3:-4:2]))
    N2_cmap = (0.6, 0.6, 0.6) 
    strain_cmap = [cmap[idx[0]], N2_cmap]
    strain_lut = dict(zip([candidate_gene,
                           CONTROL_STRAIN],
                          strain_cmap))
    
    stim_cmap = sns.color_palette('Pastel1',3)
    stim_lut = dict(zip(STIMULI_ORDER.keys(), stim_cmap))
    feat_lut = {f:v for f in featlist for k,v in stim_lut.items() if k in f}

    return strain_lut, stim_lut, feat_lut

def plot_colormaps(strain_lut, stim_lut, saveto):
    """
    

    Parameters
    ----------
    strain_lut : TYPE
        DESCRIPTION.
    stim_lut : TYPE
        DESCRIPTION.
    saveto : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    sns.set_style('dark')
    fig, ax = plt.subplots(1,2, figsize = [10,2])
    # ax=ax.flatten()
    for c, (axis, lut) in enumerate([(ax[0], strain_lut), (ax[1], stim_lut)]):
        axis.imshow([[v for v in lut.values()]])
        axis.axes.set_xticks(range(0, len(lut), 1))
        axis.axes.set_xticklabels(lut.keys(), rotation=90, fontsize=12)
        axis.axes.set_yticklabels([])
    fig.tight_layout()
    if saveto==None:
        return
    else:
        fig.savefig(saveto / 'colormaps.png')
        
    return
    
    
# def write_ordered_features(clusterfeats, saveto):
#     Path(saveto).touch(exist_ok=False)
#     with open(saveto, 'w') as fid:
#         for l in clusterfeats:
#             fid.writelines(l + ',\n') 
            
#     return

def make_clustermaps(featZ, meta, featsets, strain_lut, feat_lut, saveto, group_vars=['worm_gene','imaging_date_yyyymmdd']):
    """
    

    Parameters
    ----------
    featZ : TYPE
        DESCRIPTION.
    meta : TYPE
        DESCRIPTION.
    featsets : TYPE
        DESCRIPTION.
    strain_lut : TYPE
        DESCRIPTION.
    feat_lut : TYPE
        DESCRIPTION.
    saveto : TYPE
        DESCRIPTION.
    group_vars : TYPE, optional
        DESCRIPTION. The default is ['worm_gene','imaging_date_yyyymmdd'].

    Returns
    -------
    clustered_features : TYPE
        DESCRIPTION.

    """
    featZ_grouped = pd.concat([featZ,
                               meta],
                              axis=1
                              ).groupby(group_vars).mean()
    featZ_grouped.reset_index(inplace=True)
    
    row_colors = featZ_grouped['worm_gene'].map(strain_lut)
    col_colors = featZ_grouped[featsets['all']].columns.map(feat_lut)
    
    # make clustermaps
    clustered_features = {}
    
    for stim, fset in featsets.items():
        cg = sns.clustermap(featZ_grouped[fset],
                        row_colors=row_colors,
                        col_colors=col_colors,
                        vmin=-2,
                        vmax=2
                        )
        cg.ax_heatmap.axes.set_xticklabels([])
        cg.ax_heatmap.axes.set_yticklabels([])
        if saveto!=None:          
            cg.savefig(Path(saveto) / '{}_clustermap.png'.format(stim))
        
        clustered_features[stim] = np.array(featsets[stim])[cg.dendrogram_col.reordered_ind]
        plt.close('all')
        
    return clustered_features
        