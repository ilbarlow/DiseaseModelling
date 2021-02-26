#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 13:51:04 2021

@author: ibarlow

Script for making introductory/figure 1 of paper:
    - Overview of phenotypic space
        hierachical clustergram
        PCA
    - "robustness" of phenotype by looking at day to day variation of N2 and 
    compare to strains: if N2 > strains then hard to detect phenotype. Calculate the 
    coefficient of variation (standard deviation / mean) for each strain
    
"""

import pandas as pd
# import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
from sklearn.decomposition import PCA
from tierpsytools.preprocessing.preprocess_features import impute_nan_inf

sys.path.insert(0, '/Users/ibarlow/Documents/GitHub/pythonScripts/DiseaseModelling/hydra_screen/phenotype_summary')
from helper import (read_disease_data,
                    filter_features,
                    make_colormaps,
                    strain_gene_dict,
                    select_strains,
                    BLUELIGHT_WINDOW_DICT,
                    DATES_TO_DROP,
                    STIMULI_ORDER)
from plotting_helper import (CUSTOM_STYLE,
                             make_clustermaps,
                             plot_cmap_text,
                             feature_box_plots)

ROOT_DIR = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen')
FEAT_FILE =  Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/eleni_filters/features_summary_tierpsy_plate_filtered_traj_compiled.csv') #list(ROOT_DIR.rglob('*filtered/features_summary_tierpsy_plate_20200930_125752.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/eleni_filtered/features_summary_tierpsy_plate_filtered_traj_compiled.csv')# #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/eleni_filters/filenames_summary_tierpsy_plate_filtered_traj_compiled.csv')#list(ROOT_DIR.rglob('*filtered/filenames_summary_tierpsy_plate_20200930_125752.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/eleni_filtered/filenames_summary_tierpsy_plate_filtered_traj_compiled.csv')#  #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = list(ROOT_DIR.rglob('*wells_annotated_metadata.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')

RAW_DATA_DIR = Path('/Volumes/Ashur Pro2/DiseaseScreen')
WINDOW_FILES = RAW_DATA_DIR / 'Results' / 'window_summaries'

CONTROL_STRAIN = 'N2'

EXAMPLES = {'curvature_head_norm_abs_50th_prestim': ['bbs-1',
                                                      'avr-14'],
             'speed_midbody_norm_50th_prestim': ['glr-1',
                                                 'unc-80'],
             'relative_to_body_radial_velocity_head_tip_w_forward_10th_poststim': ['cat-2',
                                                                                  'glc-2']
             }

saveto = ROOT_DIR / 'Figures' / 'summary_figures'
saveto.mkdir(exist_ok=True)

#%%
if __name__ == '__main__':

    #set style for all figures
    plt.style.use(CUSTOM_STYLE)
    sns.set_style('ticks')

    feat, meta = read_disease_data(FEAT_FILE,
                                   FNAME_FILE,
                                   METADATA_FILE,
                                   export_nan_worms=False)
    
    meta.worm_gene.replace({'C43B7.2':'figo-1'}, inplace=True)
    
    #make colormaps
    genes = [g for g in meta.worm_gene.unique() if g != CONTROL_STRAIN]
    genes = [g for g in genes if 'myo' not in g and 'unc-54' not in g]
    genes.sort()
    # strains = genes.copy()
    # strains.extend([CONTROL_STRAIN])
    
    # select only neuro disease models
    meta = meta.query('@genes in worm_gene or @CONTROL_STRAIN in worm_gene')
    feat = feat.loc[meta.index,:]
     
    # filter features
    feat_df, meta_df, featsets = filter_features(feat,
                                                 meta)
    
    strain_lut, stim_lut, feat_lut = make_colormaps(genes,
                                                    featlist=featsets['all']
                                                    )
    
    plot_cmap_text(strain_lut, 70)
    plt.savefig(saveto / 'strain_cmap.png', bbox_inches="tight", dpi=200)
    plt.close('all')
    
    # impute nans and inf and z score  
    feat_nonan = impute_nan_inf(feat_df)

    featZ = pd.DataFrame(data=stats.zscore(feat_nonan[featsets['all']], axis=0),
                         columns=featsets['all'],
                         index=feat_nonan.index)

    assert featZ.isna().sum().sum() == 0
    
    #%% make clustermaps
    (saveto / 'clustermaps').mkdir(exist_ok=True)
    (saveto / 'clustermaps' / 'day_grouped').mkdir(exist_ok=True)
    
    clustered_features = make_clustermaps(featZ,
                                          meta_df,
                                          featsets,
                                          strain_lut,
                                          feat_lut,
                                          group_vars=['worm_gene', 'imaging_date_yyyymmdd'],
                                          saveto= saveto / 'clustermaps')
    plt.close('all')
   
    #%% now get an idea of N2 robustness - measure CoV
    
    # genes.append(CONTROL_STRAIN)
    
    cov_by_date = []
    for g in genes:
        _meta = meta_df.query('@g in worm_gene')
        
        _feat_grouped = pd.concat([feat_nonan.loc[_meta.index],
                             _meta], axis=1).groupby('date_yyyymmdd').mean()
        
        _feat_cov = _feat_grouped[featsets['all']].apply(stats.variation).apply(np.abs)
        # pd.Series(np.abs(stats.variation(_feat_grouped[featsets['all']])),
        #                       index=featsets['all'])
        _feat_cov['worm_gene'] = g
        cov_by_date.append(_feat_cov.to_frame().transpose())

    cov_by_date = pd.concat(cov_by_date)
    cov_by_date.set_index('worm_gene',inplace=True)
    
    #find mean and std for cov for N2 by resampling
    n_per_strains = meta.query('@genes in worm_gene').groupby('worm_gene').apply(len)
    
    #10 iterations of resampling
    control_cov = []
    for i in range(0,10):
        _meta = meta_df.query('@CONTROL_STRAIN in worm_gene').sample(int(n_per_strains.mean()))
        
        _feat_grouped = pd.concat([feat_nonan.loc[list(_meta.index),:],
                             _meta], axis=1).groupby('date_yyyymmdd').mean()
        
        _feat_cov = _feat_grouped[featsets['all']].apply(stats.variation).apply(np.abs)
        # pd.Series(np.abs(stats.variation(_feat_grouped[featsets['all']])),
        #                       index=featsets['all'])
        control_cov.append(_feat_cov.to_frame().transpose())
    
    control_cov = pd.concat(control_cov)
    
    fig, ax = plt.subplots(figsize=(10,8))
    # plt.errorbar(cov_by_date.index, cov_by_date.mean(axis=1), cov_by_date.std(axis=1)/(cov_by_date.shape[1]**-2))
    cov_by_date.mean(axis=1).plot(color=(0.2, 0.2, 0.2))
    ax.set_xticks(range(0, len(genes)))  
    ax.set_xticklabels(labels=genes,
                       rotation=90)
    ax.set_ylabel('Coefficient of Variation')
    
    plt.plot(np.zeros(len(genes))+control_cov.mean(axis=1).mean(),
             '--', color=strain_lut[CONTROL_STRAIN])
    plt.fill_between([x.get_text() for x in ax.axes.xaxis.get_ticklabels()],
                     np.zeros(len(genes))+control_cov.mean(axis=1).mean() - control_cov.mean(axis=1).std(),
                     np.zeros(len(genes))+control_cov.mean(axis=1).mean() + control_cov.mean(axis=1).std(),
                     color=strain_lut[CONTROL_STRAIN],
                     alpha=0.7)
    plt.tight_layout()
    plt.savefig(saveto / 'coefficient_of_variation.png', dpi=300)
    plt.savefig(saveto / 'coefficient_of_variation.svg', dpi=300)

    # %% PCA plots by strain - prestim, bluelight and poststim separately
    
    # do PCA on entire space and plot worms as they travel through
    # long form featmat
    long_featmat = []
    for stim,fset in featsets.items():
        if stim != 'all':
            _featmat = pd.DataFrame(data=feat.loc[:,fset].values,
                                    columns=['_'.join(s.split('_')[:-1])
                                               for s in fset],
                                    index=feat.index)
            _featmat['bluelight'] = stim
            _featmat = pd.concat([_featmat,
                                  meta.loc[:,'worm_gene']],
                                 axis=1)
            long_featmat.append(_featmat)
    long_featmat = pd.concat(long_featmat,
                             axis=0)
    long_featmat.reset_index(drop=True,
                             inplace=True)
    
    full_fset = list(set(long_featmat.columns) - set(['worm_gene', 'bluelight']))
    long_feat_nonan = impute_nan_inf(long_featmat[full_fset])

    long_meta = long_featmat[['worm_gene', 'bluelight']]
    long_featmatZ = pd.DataFrame(data=stats.zscore(long_feat_nonan[full_fset], axis=0),
                                 columns=full_fset,
                                 index=long_feat_nonan.index)
    
    assert long_featmatZ.isna().sum().sum() == 0
    #%% PCA
    pca = PCA()
    X2=pca.fit_transform(long_featmatZ.loc[:,full_fset])

    cumvar = np.cumsum(pca.explained_variance_ratio_)
    thresh = cumvar <= 0.95 #set 95% variance threshold
    cut_off = int(np.argwhere(thresh)[-1])

    #make a plot
    plt.figure()
    plt.plot(range(0, len(cumvar)), cumvar*100)
    plt.plot([cut_off,cut_off], [0, 100], 'k')
    plt.xlabel('Number of Principal Components')
    plt.ylabel('variance explained')
    plt.tight_layout()
    plt.savefig(saveto / 'long_df_variance_explained.png', dpi =150)
    
    #now put the 1:cut_off PCs into a dataframe
    PCname = ['PC_%d' %(p+1) for p in range(0,cut_off+1)]
    PC_df = pd.DataFrame(data=X2[:,:cut_off+1],
                         columns=PCname,
                         index=long_featmatZ.index)

    PC_plotting = pd.concat([PC_df,
                             long_meta[['worm_gene',
                                            'bluelight']]],
                             axis=1)
    
    # groupby worm gene to see the trajectory through PC space
    PC_plotting_grouped = PC_plotting.groupby(['worm_gene',
                                               'bluelight']).mean().reset_index()
    PC_plotting_grouped['stimuli_order'] = PC_plotting_grouped['bluelight'].map(STIMULI_ORDER)
    PC_plotting_grouped.sort_values(by=['worm_gene',
                                        'stimuli_order'],
                                    ascending=True,
                                    inplace=True)

    #%% nice figures
    plt.figure(figsize = [14,12])
    sns.lineplot(x='PC_1',
                    y='PC_2',
                    data=PC_plotting_grouped,
                    hue='worm_gene',
                    hue_order=strain_lut.keys(),
                    palette=strain_lut,
                    alpha=0.8,
                    legend=False,
                    sort=False)
    sns.scatterplot(x='PC_1',
                    y='PC_2',
                    data=PC_plotting_grouped,
                    hue='worm_gene',
                    style='bluelight',
                    style_order=STIMULI_ORDER.keys(),
                    hue_order=strain_lut.keys(),
                    palette=strain_lut,
                    linewidth=0,
                    s=70)  
    plt.xlim(-15,45)
    plt.ylim(-15,45)
    # plt.axis('equal')
    plt.legend(loc='right', bbox_to_anchor=(0.75, 0.25, 0.5, 0.5), fontsize='large')
    plt.xlabel('PC_1 ({}%)'.format(np.round(pca.explained_variance_ratio_[0]*100,2)))
    plt.ylabel('PC_2 ({}%)'.format(np.round(pca.explained_variance_ratio_[1]*100,2)))                                 
    plt.tight_layout()
    plt.savefig(saveto / 'PC1PC2_trajectory_space.png', dpi=400)
    # plt.savefig(saveto / 'PC1PC2_trajectory_space.svg', dpi=400)

    
    plt.figure(figsize = [14,12])
    sns.lineplot(x='PC_3',
                y='PC_4',
                data=PC_plotting_grouped,
                hue='worm_gene',
                hue_order=strain_lut.keys(),
                palette=strain_lut,
                alpha=0.8,
                legend=False,
                sort=False)
    sns.scatterplot(x='PC_3',
                    y='PC_4',
                    data=PC_plotting_grouped,
                    hue='worm_gene',
                    style='bluelight',
                    style_order=STIMULI_ORDER.keys(),
                    hue_order=strain_lut.keys(),
                    palette=strain_lut,
                    linewidth=0,
                    s=70)  
    plt.xlim(-10,20)
    plt.ylim(-20,10)
    # plt.axis('equal')
    plt.legend(loc='right', bbox_to_anchor=(0.75, 0.25, 0.5, 0.5), fontsize='large')
    plt.xlabel('PC_3 ({}%)'.format(np.round(pca.explained_variance_ratio_[2]*100,2)))
    plt.ylabel('PC_4 ({}%)'.format(np.round(pca.explained_variance_ratio_[3]*100,2)))                                 
    plt.tight_layout()
    plt.savefig(saveto / 'PC3PC4_trajectory_space.png', dpi=400)
    # plt.savefig(saveto / 'PC3PC4_trajectory_space.svg', dpi=400)
    
     #%% plot the example features
    # and make nice plots of the selected figures
    
    for k,v in EXAMPLES.items():
        examples_feat_df, examples_meta_df, dx, gene_list = select_strains(v,
                                                                            CONTROL_STRAIN,
                                                                            feat_df=feat,
                                                                            meta_df=meta)
    
        # filter features
        examples_feat_df, examples_meta_df, featsets = filter_features(examples_feat_df,
                                                     examples_meta_df)


        examples_strain_lut = make_colormaps(gene_list,
                                                featlist=featsets['all'],
                                                idx=dx,
                                                candidate_gene=v
                                                )
        examples_strain_lut = examples_strain_lut[0]
    
        feature_box_plots(k,
                          examples_feat_df,
                          examples_meta_df,
                          examples_strain_lut,
                          show_raw_data=False,
                          add_stats=True)
        plt.savefig(saveto / '{}_boxplot.png'.format(k),
                    dpi=200)
        # plt.close('all')
