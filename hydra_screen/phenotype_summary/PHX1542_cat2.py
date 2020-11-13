#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 12:55:01 2020

@author: ibarlow
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

from helper import read_disease_data, select_strains, filter_features, make_colormaps, plot_colormaps, make_clustermaps

FEAT_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')

CONTROL_STRAIN = 'N2'
CANDIDATE_GENE='cat-2'

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
                                                    featlist=featsets['all'])
    
    # colorbars to map colors to strains
    plot_colormaps(strain_lut, stim_lut, SAVETO)
    
    #%%  fillnans and normalise - used later for clustering and PCA/LDA
    feat_nonan = feat.copy()
    feat_nonan.fillna(feat_nonan.mean(axis=0),
                      inplace=True)

    featZ = pd.DataFrame(data=stats.zscore(feat_nonan, axis=0),
                         columns=featsets['all'],
                         index=feat_nonan.index)
    
    assert featZ.isna().sum().sum() == 0
    #%% make a nice clustermap / heatmap
    
    clustered_features = make_clustermaps(featZ, meta, featsets, strain_lut, feat_lut, saveto=None)
    
        
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
                            palette=strain_lut.values(),
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

    print('{} out of {} significant'.format(bhP_values.notna().sum().sum(), bhP_values.shape[1]))
    #%%
    # nice figure with single barcode for each strains, and asterisks of signficicantly different features
    selected_feats = ['relative_to_body_radial_velocity_head_tip_w_forward_10th_bluelight',
                  'speed_head_tip_w_forward_50th_bluelight',
                  'speed_head_tip_w_forward_IQR_prestim']
    
    font_settings = {'fontsize':14}
    
    for stim, fset in clustered_features.items():
        
        heatmap_df = [pd.concat([featZ,
                                meta],
                               axis=1
                               ).groupby('worm_gene').mean()[fset]]
                                                             
        heatmap_df.append(np.log10(bhP_values[fset]))
        
        _stim = pd.DataFrame(data=[i.split('_')[-1] for i in fset],
                           columns=['stim_type'])
        _stim['stim_type'] = _stim['stim_type'].map(STIMULI_ORDER)
        _stim = _stim.transpose()
        _stim.rename(columns={c:v for c,v in enumerate(fset)}, inplace=True)
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
        f.savefig(SAVETO / '{}_heatmap.png'.format(stim))
        
    #TODO make nice plots of the bluelight responese - either ts or windows
    #%%
    window_files = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/Results/window_summaries')
    bluelight_window_dict = {0: 'prelight',
                         1: 'bluelight',
                         2: 'postlight',
                         3: 'prelight',
                         4: 'bluelight',
                         5: 'postlight',
                         6: 'prelight',
                         7: 'bluelight',
                         8: 'postlight'}
    
    def find_window(fname):
        import re
        window_regex = r"(?<=_window_)\d{0,9}"
        window = int(re.search(window_regex, str(fname))[0])
        return window
    
    summary_files = list(window_files.rglob('*_window_*'))   
    feat_files = [f for f in summary_files if 'features' in str(f)]
    feat_files.sort(key=find_window)
    fname_files = [f for f in summary_files if 'filenames' in str(f)]
    fname_files.sort(key=find_window)
    
    assert (find_window(f[0]) == find_window(f[1]) for f in list(zip(feat_files, fname_files)))
    
    feat_df = []
    meta_df = []
    for c,f in enumerate(list(zip(feat_files, fname_files))):
        _feat, _meta = read_disease_data(f[0], f[1], METADATA_FILE, drop_nans=False)
        _meta['window'] = find_window(f[0])
        
        meta_df.append(_meta)
        feat_df.append(_feat)
     
    meta = pd.concat(meta_df)
    meta.reset_index(drop=True, inplace=True)
    feat = pd.concat(feat_df)
    feat.reset_index(drop=True, inplace=True)
    
    
    feat, meta = drop_nan_worms(feat, meta, saveto=None)
    feat_df, meta_df, idx, gene_list = select_strains(candidate_gene, control_strain, feat_df, meta_df)
        
    
        