#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 16:23:42 2020

@author: ibarlow

script to investigate cat-4 phenotypes

Import just cat-4 and N2 data to compare the features and find 
independent, simple and concise features

Look at prestim and poststim independentally
"""


import pandas as pd
import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from matplotlib.gridspec import GridSpec
from tierpsytools.analysis.significant_features import k_significant_feat      

sys.path.insert(0, '/Users/ibarlow/Documents/GitHub/pythonScripts/DiseaseModelling/hydra_screen/phenotype_summary')

from helper import read_disease_data, select_strains, filter_features, make_colormaps, write_ordered_features, STIMULI_ORDER

FEAT_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')

CONTROL_STRAIN = 'N2'
CANDIDATE_GENE='cat-4'

SAVETO = FEAT_FILE.parent.parent.parent / 'Figures' / 'paper_figures' / CANDIDATE_GENE
SAVETO.mkdir(exist_ok=True)
feat256_fname = Path('/Users/ibarlow/tierpsy-tools-python/tierpsytools/extras/feat_sets/tierpsy_256.csv')

# # setting up features set types
# feat256 = []
# with open(feat256_fname, 'r', encoding='utf-8-sig') as fid:
#     for l in fid.readlines():
#         feat256.append(l.rstrip().lstrip())

#%%
if __name__ == '__main__':
    
    feat, meta = read_disease_data(FEAT_FILE,
                                   FNAME_FILE,
                                   METADATA_FILE,
                                   export_nan_worms=False)
    
    feat, meta, idx, gene_list = select_strains(CANDIDATE_GENE,
                                                CONTROL_STRAIN,
                                                feat,
                                                meta)
    
    feat, meta, featsets = filter_features(feat,
                                           meta)
    
    strain_lut, stim_lut, feat_lut = make_colormaps(gene_list,
                                                    idx,
                                                    CANDIDATE_GENE,
                                                    CONTROL_STRAIN,
                                                    featsets['all'])
    
    #%% colorbars to map colors to strains
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
    # stim_cmap = 
    ax[1].imshow([[v for v in stim_lut.values()]])
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
    
    assert featZ.isna().sum().sum() == 0
    #%% make a nice clustermap / heatmap
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
    
    # make clustermaps
    clustered_features = {}
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
    clustered_features['all'] = np.array(featlist)[cg.dendrogram_col.reordered_ind]
    # write_ordered_features(_reordered_feats,
    #                        SAVETO / 'all_clusterfeats_ordered.txt')
    
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
        # cg.savefig(SAVETO / '{}_clustermap.png'.format(fset))
        
        clustered_features[fset] = np.array(featsets[fset])[cg.dendrogram_col.reordered_ind]
        # write_ordered_features(_reordered_feats,
        #                        SAVETO / '{}_clusterfeats_ordered.txt'.format(fset))               
        plt.close('all')
        
    #%% k significant features for each prestim, bluelight and poststim
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
            plt.savefig(SAVETO / '{}_{}_ksig_feats.png'.format(i*20, stim),
                        dpi=400)
    plt.close('all')
    
    #%% pairwise statistics to find features that are different from N2
    
    from tierpsytools.analysis.paired_stats_tests import paired_stats_tests
    
    pVals, bhP_values, group_classes = paired_stats_tests(feat,
                                                           meta.loc[:, 'worm_gene'],
                                                           control_group='N2')
    bhP_values['worm_gene'] = group_classes
    bhP_values.rename(mapper={0:'p<0.05'},
                      inplace=True)
    # bhP_values.set_index('worm_gene', inplace=True)
    
    #%%
    # nice figure with single barcode for each strains, and asterisks of signficicantly different features
    selected_feats = ['length_50th_poststim',
                  'curvature_std_midbody_abs_50th_poststim',
                  'd_curvature_std_midbody_abs_50th_poststim']
    
    font_settings = {'fontsize':14}
    
    for fset, values in clustered_features.items():
        
        heatmap_df = [pd.concat([featZ,
                                meta],
                               axis=1
                               ).groupby('worm_gene').mean()[values]]
                                                             
        heatmap_df.append(np.log10(bhP_values[values]))
        
        _stim = pd.DataFrame(data=[i.split('_')[-1] for i in values],
                           columns=['stim_type'])
        _stim['stim_type'] = _stim['stim_type'].map(STIMULI_ORDER)
        _stim = _stim.transpose()
        _stim.rename(columns={c:v for c,v in enumerate(values)}, inplace=True)
        heatmap_df.append(_stim)
        
        heatmap_df = pd.concat(heatmap_df)
        # heatmap_df.rename(mapper={0:'p<0.05'},
        #                   inplace=True)
       
        cm = ['inferno', 'inferno', 'gray', 'Pastel1'] #'Reds', 'Greens', 'Oranges', 'Purples', 'bone', 'winter']
        vmin_max = [(-2,2), (-2,2), (-20, 0), (1,3)]
        f = plt.figure(figsize= (20,5))
        gs = GridSpec(4, 1, wspace=0, hspace=0, height_ratios=[3,3,1, 1])
        cbar_ax = f.add_axes([.91, .3, .03, .4])
    
        # f, axs = plt.subplots(heatmap_df.shape[0],1, gridspec_kw={'wspace': 0})
        for n, ((ix,r), c, v) in enumerate(zip(heatmap_df.iterrows(), cm, vmin_max)):
    
            axis = f.add_subplot(gs[n])
            sns.heatmap(r.to_frame().transpose().astype(float),
                        yticklabels=[ix],
                        xticklabels=[],
                        ax=axis,
                        cmap=c,
                        cbar=n==0, #only plots colorbar for first plot
                        cbar_ax=None if n else cbar_ax,
                        vmin=v[0],
                        vmax=v[1])
            axis.set_yticklabels(labels=[ix], rotation=0, fontsize=14)
            
            if n>2:
                c = sns.color_palette('Pastel1',3)
                sns.heatmap(r.to_frame().transpose(),
                        yticklabels=[ix],
                        xticklabels=[],
                        ax=axis,
                        cmap=c,
                        cbar=n==0, #only plots colorbar for first plot
                        cbar_ax=None if n else cbar_ax,
                        vmin=v[0],
                        vmax=v[1])
                axis.set_yticklabels(labels=[ix], rotation=0, fontsize=14)
        f.tight_layout(rect=[0, 0, .9, 1])

        for sf in selected_feats:
            try:
                axis.text(heatmap_df.columns.get_loc(sf), 1, '*', fontdict=font_settings)
            except KeyError:
                print('{} not in featureset'.format(sf))
        f.savefig(SAVETO / '{}_heatmap.png'.format(fset))
