#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:47:23 2020

@author: ibarlow

script to investigate cat-4 phenotypes

Import just cat-4 and N2 data to compare the features and find 
independent, simple and concise features

Look at prestim and poststim independentally

"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from tierpsytools.read_data.hydra_metadata import read_hydra_metadata, align_bluelight_conditions
from sklearn.decomposition import PCA


FEAT_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')
SAVETO = FEAT_FILE.parent.parent.parent / 'Figures' / 'paper_figures'
SAVETO.mkdir(exist_ok=True)
feat256_fname = Path('/Users/ibarlow/tierpsy-tools-python/tierpsytools/extras/feat_sets/tierpsy_256.csv')

# setting up features set types
feat256 = []
with open(feat256_fname, 'r', encoding='utf-8-sig') as fid:
    for l in fid.readlines():
        feat256.append(l.rstrip().lstrip())

stimuli_order = {'prestim':1,
                 'bluelight':2,
                 'poststim':3}
bluelight_feat256 = ['{}_{}'.format(f, b) for f in feat256 for b in stimuli_order.keys()]
featset256 ={}
for b in stimuli_order.keys():
    featset256[b] = ['{}_{}'.format(f,b) for f in feat256]
    
# other vars
CONTROL_STRAIN = 'N2'
DATES_TO_DROP = '20200626'
BAD_FEAT_THRESH = 3 # 3 standard deviations away from the mean
BAD_FEAT_FILTER = 0.1 # 10% threshold for removing bad features
BAD_WELL_FILTER = 0.3 # 30% threshold for bad well
# NAN_CUTOFF = 7000

candidate_gene = 'cat-4'
SAVETO = SAVETO / candidate_gene
SAVETO.mkdir(exist_ok=True)

# %% helper functions
def write_ordered_features(clusterfeats, saveto):
    Path(saveto).touch(exist_ok=False)
    with open(saveto, 'w') as fid:
        for l in clusterfeats:
            fid.writelines(l + ',\n') 
            
    return

#%%

if __name__ == '__main__':

    feat = pd.read_csv(FEAT_FILE,
                       comment='#')
    fname = pd.read_csv(FNAME_FILE,
                        comment='#')

    meta = pd.read_csv(METADATA_FILE, index_col=None)
    # meta = meta[meta.worm_strain.notna()]

    assert meta.worm_strain.unique().shape[0] == meta.worm_gene.unique().shape[0]
    
    #check files in summary vs files in metadata
    # assert len(set(fname.filename.unique()) - set(meta.imgstore_name.unique())) == 0

    imgstore_fname = ['/'.join(r.filename.split('/')[-3:-1]) for i,r in fname.iterrows()]
    missing_files = set(imgstore_fname) - set(meta.imgstore_name.unique())
    
    # missing files are from 20200626
    
    #%%
    feat, meta = read_hydra_metadata(feat, fname, meta)

    feat, meta = align_bluelight_conditions(feat,
                                            meta,
                                            how='inner') #removes wells that don't have all 3 conditions

    #%% wells to check
    nan_worms = meta[meta.worm_gene.isna()][['featuresN_filename',
                                             'well_name',
                                             'imaging_plate_id',
                                             'instrument_name',
                                             'imaging_date_yyyymmdd']]
    # nan_worms.to_csv(FEAT_FILE.parent / 'nan_worms.csv',
    #                   index=False)

    feat = feat.drop(index=nan_worms.index)
    meta = meta.drop(index=nan_worms.index)

    #%% select only N2s and strain of interest

    genes = [g for g in meta.worm_gene.unique() if g != CONTROL_STRAIN]
    genes = [g for g in genes if 'myo' not in g and 'unc-54' not in g]
       
    idx = [c for c,g in list(enumerate(genes)) if  g==candidate_gene]
    locs = list(meta.query('@candidate_gene in worm_gene').index)
    N2_locs = meta.query('@CONTROL_STRAIN in worm_gene').sample(len(locs)).index
    locs.extend(N2_locs)

     #Only do analysis on the disease strains
    meta = meta.loc[locs,:]
    feat = feat.loc[locs,:]
    
    # %% filtering
    imgst_cols = [col for col in meta.columns if 'imgstore_name' in col]
    miss = meta[imgst_cols].isna().any(axis=1)

    # remove data from dates to exclude
    bad_date = meta.date_yyyymmdd == float(DATES_TO_DROP)

    # bad wells
    good_wells_from_gui = meta.is_bad_well == False

    feat = feat.loc[good_wells_from_gui & ~bad_date & ~miss,:]
    meta = meta.loc[good_wells_from_gui & ~bad_date & ~miss,:]

    # remove features with too many nans
    feat = feat.loc[:, feat.isna().sum(axis=0)/feat.shape[0]<BAD_FEAT_FILTER]

    # remove features with std=0
    feat = feat.loc[:, feat.std(axis=0) > 1e-6 ]

    # remove wells with too many nans
    feat = feat.loc[feat.isna().sum(axis=1)/feat.shape[0]<BAD_WELL_FILTER, :]
    meta = meta.loc[feat.index,:]

     # feature sets
     # abs features no longer in tierpsy
    pathcurvature_feats = [x for x in feat.columns if 'path_curvature' in x]
    
    #remove these features
    feat = feat.drop(columns=pathcurvature_feats)
    
    featlist = list(feat.columns)
    # for f in featsets.keys():
    #     featsets[f] = [x for x in feat.columns if f in x]
    #     featsets[f] = list(set(featsets[f]) - set(pathcurvature_feats))
    
    featsets={}
    for stim in stimuli_order.keys():
        featsets[stim] = [f for f in featlist if stim in f]
    
    #%% make colormaps
    cmap = list(np.flip((sns.color_palette('cubehelix',
                                           len(genes)*2+6))[3:-4:2]))
    N2_cmap = (0.6, 0.6, 0.6) 
    strain_cmap = [cmap[idx[0]], N2_cmap]
    
    stim_cmap = sns.color_palette('Pastel1',3)
    
    strain_lut = dict(zip([candidate_gene,
                           CONTROL_STRAIN],
                          strain_cmap))
    
    stim_lut = dict(zip(stimuli_order.keys(), stim_cmap))
    feat_lut = {f:v for f in featlist for k,v in stim_lut.items() if k in f}

    # colorbar to map colors to strains
    sns.set_style('dark')
    # plt.figure(figsize=[10,4])
    fig, ax = plt.subplots(1,2, figsize = [10,2])
    # ax=ax.flatten()
    ax[0].imshow([strain_cmap])
    ax[0].axes.set_xticks(range(0, len(strain_lut), 1))
    ax[0].axes.set_xticklabels(strain_lut.keys(), rotation=90, fontsize=12)
    ax[0].axes.set_yticklabels([])
    # plt.savefig(SAVETO / 'strain_colormap.png')
    
    #colorbar to map colors to stimuli
    ax[1].imshow([stim_cmap])
    ax[1].axes.set_xticks(range(0, len(stim_lut), 1))
    ax[1].axes.set_xticklabels(stim_lut.keys(), rotation=90, fontsize=12)
    ax[1].axes.set_yticklabels([])
    fig.tight_layout()
    fig.savefig(SAVETO / 'colormaps.png')
        
    #%%  fillnans and normalise - used later for clustering and PCA/LDA
    feat_nonan = feat.copy()
    feat_nonan.fillna(feat_nonan.mean(axis=0),
                      inplace=True)

    featZ = pd.DataFrame(data=stats.zscore(feat_nonan, axis=0),
                         columns=featlist,
                         index=feat_nonan.index)
    
    #%% clustermaps for first  inspection

    groupby_vars = ['worm_gene',
                    # 'imaging_plate_id',
                    'imaging_date_yyyymmdd']
    
    featZ_grouped = pd.concat([featZ,
                               meta],
                              axis=1
                              ).groupby(groupby_vars).mean()
    featZ_grouped.reset_index(inplace=True)
    
    row_colors = featZ_grouped['worm_gene'].map(strain_lut)
    col_colors = featZ_grouped[featlist].columns.map(feat_lut)
    
    cg = sns.clustermap(featZ_grouped[featlist],
                        row_colors=row_colors,
                        col_colors=col_colors,
                        vmin=-2,
                        vmax=2
                        # dendrogram_ratio=(0.2,0.2)
                        )
    cg.ax_heatmap.axes.set_xticklabels([])
    cg.ax_heatmap.axes.set_yticklabels([])
    cg.savefig(SAVETO / 'allstim_clustermap.png')
    
    _reordered_feats = np.array(featlist)[cg.dendrogram_col.reordered_ind]
    write_ordered_features(_reordered_feats,
                           SAVETO / 'all_clusterfeats_ordered.txt')
    
    for fset in featsets.keys():
        cg = sns.clustermap(featZ_grouped[featsets[fset]],
                        row_colors=row_colors,
                        # col_colors=col_colors,
                        vmin=-2,
                        vmax=2
                        # dendrogram_ratio=(0.2,0.2)
                        )
        cg.ax_heatmap.axes.set_xticklabels([])
        cg.ax_heatmap.axes.set_yticklabels([])
        cg.savefig(SAVETO / '{}_clustermap.png'.format(fset))
        
        _reordered_feats = np.array(featsets[fset])[cg.dendrogram_col.reordered_ind]
        write_ordered_features(_reordered_feats,
                               SAVETO / '{}_clusterfeats_ordered.txt'.format(fset))               
        plt.close('all')
    
    #%% k significant features for each prestim, bluelight and poststim
    from tierpsytools.analysis.significant_features import k_significant_feat      
    sns.set_style('white')
    label_format = '{:.4f}'
    kfeats = {}
    for stim, fset in featsets.items():
        kfeats[stim], scores, support = k_significant_feat(
            feat_nonan[fset],
            meta.worm_gene,
            k=100,
            plot=False)
        
        for i in range(0,5):
            fig, ax = plt.subplots(4, 5, sharex=True, figsize = [20,20])
            for c, axis in enumerate(ax.flatten()):
                counter=(20*i)-(20-c)
                sns.boxplot(x=meta['worm_gene'],
                            y=feat[kfeats[stim][counter]],
                            palette=strain_cmap,
                            ax=axis)
                axis.set_ylabel(fontsize=8, ylabel=kfeats[stim][counter])
                axis.set_yticklabels(labels = [label_format.format(x) for x in axis.get_yticks()], fontsize=6)
                axis.set_xlabel('')
            plt.tight_layout()
            fig.fontsize=11
            plt.savefig(SAVETO / '{}_ksig_feats.png'.format(i*20),
                        dpi=400)
    plt.close('all')
    # now for paired stats tests
    from tierpsytools.analysis.paired_stats_tests import paired_stats_tests
    
    pVals, bhP_values, group_classes = paired_stats_tests(feat,
                                                           meta.loc[:, 'worm_gene'],
                                                           control_group='N2')
    bhP_values['worm_gene'] = group_classes
    bhP_values.set_index('worm_gene', inplace=True)

    bhP_values.notna().sum().sum()
    
    bhP_values[['length_50th_poststim',
               'curvature_std_midbody_abs_50th_poststim',
               'd_curvature_std_midbody_abs_50th_poststim']]
    # # %% clustermaps
    

    
    # row_colors = meta['worm_gene'].map(lut)
   
    # cg = sns.clustermap(feat[featsets['prestim']],
    #                    z_score=1,
    #                    row_colors=row_colors,
    #                    vmin=-2,
    #                    vmax=2,
    #                    method='average')
    
    # cg.dendrogram_col.reordered_ind
    # reordered_feats = np.array(featsets['prestim'])[cg.dendrogram_col.reordered_ind]


    # groupby_vars = ['worm_gene',
    #                 # 'imaging_plate_id',
    #                 'imaging_date_yyyymmdd']
    
    # featZ_grouped = pd.concat([featZ,
    #                                 meta],
    #                                axis=1
    #                                ).groupby(groupby_vars).mean()
    # featZ_grouped.reset_index(inplace=True)
    
    # row_colors = featZ_grouped['worm_gene'].map(lut)
    # col_colors = featZ_grouped[featlist].columns.map(feat_lut)
    
    # cg = sns.clustermap(featZ_grouped[featlist],
    #                     row_colors=row_colors,
    #                     col_colors=col_colors,
    #                     vmin=-2,
    #                     vmax=2
    #                     # dendrogram_ratio=(0.2,0.2)
    #                     )
    # cg.ax_heatmap.axes.set_xticklabels([])
    # cg.ax_heatmap.axes.set_yticklabels([])
    
    # cg.dendrogram_col.reordered_ind
    # reordered_feats = np.array(featsets['prestim'])[cg.dendrogram_col.reordered_ind]

    # sns.heatmap(featZ_grouped[reordered_feats],
    #             vmin=-2,
    #             vmax=2, yticklabels=featZ_grouped['worm_gene'])
    
    # label_format = '{:.4f}'
    # cor={}
    # cor_feats={}
    # for fset in featsets.keys():
    #     cor[fset] = feat[featsets[fset]].corr()
    # # want to keep features that are not correlated and simple
    #     cor_feats[fset]= cor[fset].columns[(cor[fset]>0.7).sum()<(cor[fset]>0.7).sum().mean()]
    
    # cor = feat.corr()
    # sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
    # plt.show()
    
