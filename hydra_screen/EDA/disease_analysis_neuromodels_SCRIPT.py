#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:11:21 2020

@author: ibarlow
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from tierpsytools.read_data.hydra_metadata import read_hydra_metadata, align_bluelight_conditions
from sklearn.decomposition import PCA

# analysis type
subgroup = 'neuromodels'
k_feat = False
pairwise_stats = False

FEAT_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')
SAVETO = FEAT_FILE.parent.parent.parent / 'Figures' / subgroup
SAVETO.mkdir(exist_ok=True)
feat256_fname = Path('/Users/ibarlow/tierpsy-tools-python/tierpsytools/extras/feat_sets/tierpsy_256.csv')

# features set types
feat256 = []
with open(feat256_fname, 'r', encoding='utf-8-sig') as fid:
    for l in fid.readlines():
        feat256.append(l.rstrip().lstrip())
featsets = {'prestim': [],
            'bluelight': [],
            'poststim': []}
bluelight_feat256 = ['{}_{}'.format(f, b) for f in feat256 for b in featsets.keys()]
featset256 ={}
for b in featsets.keys():
    featset256[b] = ['{}_{}'.format(f,b) for f in feat256]
    
stimuli_order = {'prestim':1,
                 'bluelight':2,
                 'poststim':3}

# vars
CONTROL_STRAIN = 'N2'
DATES_TO_DROP = '20200626'
BAD_FEAT_THRESH = 3 # 3 standard deviations away from the mean
BAD_FEAT_FILTER = 0.1 # 10% threshold for removing bad features
BAD_WELL_FILTER = 0.3 # 30% threshold for bad well
# NAN_CUTOFF = 7000


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

    #%%
    # strain sets
    genes = [g for g in meta.worm_gene.unique() if g != CONTROL_STRAIN]

    #myosin genes
    myosin_genes = [s for s in genes if 'myo' in s]
    myosin_genes.extend([s for s in genes if 'unc-54' in s])
    myosin_genes.sort()
    # myosin_genes.append(CONTROL_STRAIN)
    myo_locs = meta.query('@myosin_genes in worm_gene').index

    if subgroup == 'neuromodels':
        # Neuro disease strains
        neuro_genes = list(set(genes) - set(myosin_genes))
        neuro_genes.sort()
    
        neuro_locs = list(meta.query('@neuro_genes in worm_gene').index)
        N2_locs = meta.query('@CONTROL_STRAIN in worm_gene').sample(300).index
    
        neuro_locs.extend(N2_locs)
        neuro_genes.append(CONTROL_STRAIN)
    
        cmap = list(np.flip(sns.color_palette('cubehelix', len(neuro_genes)-1)))
        cmap.append((0.6, 0.6, 0.6))
    
         #Only do analysis on the disease strains
        meta = meta.loc[neuro_locs,:]
        feat = feat.loc[neuro_locs,:]

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
    for f in featsets.keys():
        featsets[f] = [x for x in feat.columns if f in x]
        featsets[f] = list(set(featsets[f]) - set(pathcurvature_feats))
    
    
    #%%  fillnans and normalise - used later for clustering and PCA/LDA
    feat_nonan = feat.copy()
    feat_nonan.fillna(feat_nonan.mean(axis=0),
                      inplace=True)

    featZ = pd.DataFrame(data=stats.zscore(feat_nonan, axis=0),
                         columns=feat.columns,
                         index=feat_nonan.index)

    #%% sig features
    if k_feat:
        from tierpsytools.analysis.significant_features import k_significant_feat  
        (SAVETO / 'k_sig_feats').mkdir(exist_ok=True)
        neuro_feats = {}
        for fset in featsets.keys():
            neuro_feats[fset], scores, support = k_significant_feat(
                feat_nonan[featsets[fset]],
                meta.worm_gene,
                k=20,
                plot=False)
    
        for f in neuro_feats:
            (SAVETO / 'k_sig_feats' / f).mkdir(exist_ok=True)
            for s in neuro_feats[f]:
                plt.figure()
                sns.boxplot(x=feat[s],
                            y=meta['worm_gene'],
                            order=neuro_genes,
                            orient='h',
                            showfliers=False,
                            palette = cmap)
                plt.savefig(SAVETO / 'k_sig_feats' / f / '{}.png'.format(s),
                            dpi=200)
                plt.close('all')

#%% paired stats tests
    if pairwise_stats:
        
        from tierpsytools.analysis.paired_stats_tests import paired_stats_tests
    
        pVals, bhP_values, group_classes = paired_stats_tests(feat,
                                                               meta.loc[:, 'worm_gene'],
                                                               control_group='N2')
        bhP_values['worm_gene'] = group_classes
        bhP_values.set_index('worm_gene', inplace=True)
        
        if subgroup == 'neuromodels':
            n2_index=-1
            
        (SAVETO / 'pairwise_stats_test').mkdir(exist_ok=True)
        for c,g in enumerate(neuro_genes):
            if g == 'N2':
                continue
            else:
                genes_to_plot = [g, CONTROL_STRAIN]
    
                _gene = bhP_values.columns[bhP_values.loc[g,:].notna()]
                _sort = np.argsort(bhP_values.loc[g, _gene])
                _top = _gene[_sort]
    
                (SAVETO / 'pairwise_stats_test' / g).mkdir(exist_ok=True)
                (SAVETO / 'pairwise_stats_test' / g / 'top25').mkdir(exist_ok=True)
                for f in _top[:25]:
                    plt.figure(figsize=[10,5])
                    ax1 = plt.subplot(1,2, 1)
                    ax1 = sns.boxplot(y=meta.loc[:, 'worm_gene'],
                                x=feat.loc[:, f],
                                orient='h',
                                order=neuro_genes,
                                showfliers=False,
                                palette=cmap)
                    ax2 = plt.subplot(1,2,2)
                    ax2 = sns.boxplot(y=meta.loc[meta.query('@genes_to_plot in worm_gene').index,
                                                  'worm_gene'],
                                      x=feat.loc[meta.query('@genes_to_plot in worm_gene').index,
                                                  f],
                                      orient='h',
                                      showfliers=False,
                                      palette = [cmap[c],
                                                  cmap[n2_index]],
                                      order = genes_to_plot
                                      )
                    plt.savefig(SAVETO / 'pairwise_stats_test' / g / 'top25' / \
                                '{}.png'.format(f))
                    plt.close('all')
    
                #plot bluelight comparisons
                (SAVETO / 'pairwise_stats_test' / g / 'bluelight_comparisons').mkdir(exist_ok=True)
                _bluelight = ['_'.join(f.split('_')[:-1]) for f in _top[:40] if 'bluelight' in f]
                
                id_to_merge = ['worm_gene',
                                'well_name',
                                'imaging_plate_id',
                                'imaging_date_yyyymmdd']
                for f in _bluelight:
                    flist = [f+'_{}'.format(s) for s in featsets.keys()]
                    plot_df = pd.concat([feat.loc[:, flist], meta],
                                        axis=1)
                    to_plot = plot_df.melt(id_vars=id_to_merge,
                                            value_vars=flist,
                                            value_name=f)
                    to_plot['bluelight'] = to_plot.variable.apply(lambda x: x.split('_')[-1])
    
    
                    plt.figure()
                    sns.pointplot(data=to_plot.query('@genes_to_plot in worm_gene'),
                                x='bluelight',
                                y=f,
                                hue='worm_gene',
                                ci=95,
                                palette='cubehelix')
                    plt.savefig(SAVETO / 'pairwise_stats_test' / g / 'bluelight_comparisons' /\
                                '{}_bluelightcomparisons.png'.format(f))
                    plt.close('all')

    # %% hierachical clustering
    
    lut = dict(zip(neuro_genes, cmap))
    row_colors = meta['worm_gene'].map(lut)

    # average linkage clustering of 256 features
    for fset in featset256:
        g = sns.clustermap(feat_nonan[featset256[fset]],
                           z_score=1,
                           row_colors=row_colors,
                           vmin=-1.5,
                           vmax=1.5,
                           method='average',
                           xticklabels=False,
                           yticklabels=False)
        plt.title(fset)
        # plt.tight_layout()
        plt.savefig(SAVETO / '{}_256feats_neurogenes_clustermap.png'.format(fset),
                    dpi=400)

    # all 256 together
    feat256_allconditions = list(featset256.values())
    feat256_allconditions = [i for subfeat in feat256_allconditions for i in subfeat]
    g = sns.clustermap(feat_nonan[feat256_allconditions],
                        z_score=1,
                        row_colors=row_colors,
                        vmin=-1.5,
                        vmax=1,
                        method='average',
                        xticklabels=False,
                        yticklabels=False)
    plt.savefig(SAVETO / 'allconditions_256feats_neurogenes_clustermap.png',
                dpi=400)
    
    # colorbar to map colors to strains
    sns.set_style('dark')
    plt.figure(figsize=[10,4])
    ax = plt.imshow([cmap])
    ax.axes.set_xticks(range(0, len(cmap), 1))
    ax.axes.set_xticklabels(lut.keys(), rotation=90)
    plt.savefig(SAVETO / 'colormap.png')

    # %% hierachical clustering by avearging the strains for each plate
    
    groupby_vars = ['worm_gene',
                    # 'imaging_plate_id',
                    'imaging_date_yyyymmdd']
    
    feat_nonan_grouped = pd.concat([feat, meta[groupby_vars]],
                                   axis=1)
    feat_nonan_grouped = feat_nonan_grouped.groupby(groupby_vars)
    
    feat_av_grouped = feat_nonan_grouped.mean().reset_index().drop(columns=groupby_vars[1:])
    row_colors_agg = feat_av_grouped['worm_gene'].map(lut)
    
    # average linkage clustering of 256 features
    for fset in featset256:
        g = sns.clustermap(feat_av_grouped[featset256[fset]],
                           z_score=1,
                           row_colors=row_colors_agg,
                           vmin=-1.5,
                           vmax=1.5,
                           method='average',
                           xticklabels=False,
                           yticklabels=False)
        plt.title(fset)
        # plt.tight_layout()
        plt.savefig(SAVETO / '{}_averaged_by_day_256feats_neurogenes_clustermap.png'.format(fset),
                    dpi=400)

    # all 256 together
    feat256_allconditions = list(featset256.values())
    feat256_allconditions = [i for subfeat in feat256_allconditions for i in subfeat]
    g = sns.clustermap(feat_av_grouped[feat256_allconditions],
                        z_score=1,
                        row_colors=row_colors_agg,
                        vmin=-1.5,
                        vmax=1,
                        method='average',
                        xticklabels=False,
                        yticklabels=False)
    plt.savefig(SAVETO / 'allconditions_averaged_by_day_256feats_neurogenes_clustermap.png',
                dpi=400)
    
#%% PCA of all the neuro strains

    # 256 features across bluelight, prestim and poststim
    pca = PCA()
    X2 = pca.fit_transform(featZ.loc[:,feat256_allconditions])
    cumvar = np.cumsum(pca.explained_variance_ratio_)
    thresh = cumvar <= 0.95 #set 95% variance threshold
    cut_off = int(np.argwhere(thresh)[-1])

    #make a plot
    sns.set_style('whitegrid')
    plt.figure()
    plt.plot(range(0, len(cumvar)), cumvar*100)
    plt.plot([cut_off,cut_off], [0, 100], 'k')
    plt.xlabel('Number of Principal Components', fontsize=16)
    plt.ylabel('variance explained', fontsize=16)
    plt.savefig(SAVETO / 'feat256all_variance_explained.png', dpi =150)

    #now put the 1:cut_off PCs into a dataframe
    PCname = ['PC_%d' %(p+1) for p in range(0,cut_off+1)]
    PC_df = pd.DataFrame(data=X2[:,:cut_off+1],
                         columns=PCname,
                         index=featZ.index)
    PC_plotting = pd.concat([PC_df,
                             meta],
                             axis=1)
    plt.figure(figsize = [12,12])
    sns.scatterplot(x='PC_1',
                    y='PC_2',
                    data=PC_plotting,
                    hue='worm_gene',
                    hue_order=neuro_genes,
                    palette=cmap,
                    linewidth=0,
                    alpha=0.8)
    plt.xlim(-20,60)
    plt.ylim(-30,50)
    plt.tight_layout()
    plt.xlabel('PC 1: {}%'.format(np.round(pca.explained_variance_ratio_[0]*100,2)))
    plt.ylabel('PC 2: {}%'.format(np.round(pca.explained_variance_ratio_[1]*100,2)))
    plt.savefig(SAVETO / 'feat256_allconditions_pc1pc2.png', dpi=400)

    plt.figure(figsize = [12,12])
    sns.scatterplot(x='PC_3',
                    y='PC_4',
                    data=PC_plotting,
                    hue='worm_gene',
                    hue_order=neuro_genes,
                    palette=cmap,
                    linewidth=0,
                    alpha=0.8)
    plt.xlim(-40,40)
    plt.ylim(-40,40)
    plt.tight_layout()
    plt.xlabel('PC 3: {}%'.format(np.round(pca.explained_variance_ratio_[2]*100,2)))
    plt.ylabel('PC 4: {}%'.format(np.round(pca.explained_variance_ratio_[3]*100,2)))
    plt.savefig(SAVETO / 'feat256_allconditions_pc2pc3.png', dpi=400)

    #%% pca on feature sets separated by prestim poststim and bluelight

    # do PCA on entire space and plot worms as they travel through
    long_featmat = []
    for f in featset256:
        _featmat = pd.DataFrame(data=feat.loc[:,featset256[f]].values,
                                columns=['_'.join(s.split('_')[:-1])
                                           for s in featset256[f]],
                                index=feat.index)
        _featmat['bluelight'] = f
        _featmat = pd.concat([_featmat,
                              meta.loc[:,'worm_gene']],
                             axis=1)
        long_featmat.append(_featmat)

    long_featmat = pd.concat(long_featmat,
                             axis=0)
    long_featmat.reset_index(drop=True,
                             inplace=True)
    long_featmat.loc[:,feat256] = long_featmat.fillna(long_featmat[feat256].mean(axis=0))

    long_featmatZ = pd.DataFrame(stats.zscore(long_featmat[feat256]),
                                 columns=feat256)

    
    pca = PCA()

    X2=pca.fit_transform(long_featmatZ.loc[:,feat256])

    cumvar = np.cumsum(pca.explained_variance_ratio_)
    thresh = cumvar <= 0.95 #set 95% variance threshold
    cut_off = int(np.argwhere(thresh)[-1])

    #make a plot
    sns.set_style('whitegrid')
    plt.figure()
    plt.plot(range(0, len(cumvar)), cumvar*100)
    plt.plot([cut_off,cut_off], [0, 100], 'k')
    plt.xlabel('Number of Principal Components', fontsize=16)
    plt.ylabel('variance explained', fontsize=16)
    plt.savefig(SAVETO / 'long_df_variance_explained.png', dpi =150)

    #now put the 1:cut_off PCs into a dataframe
    PCname = ['PC_%d' %(p+1) for p in range(0,cut_off+1)]
    PC_df = pd.DataFrame(data=X2[:,:cut_off+1],
                         columns=PCname,
                         index=long_featmatZ.index)

    PC_plotting = pd.concat([PC_df,
                             long_featmat[['worm_gene',
                                            'bluelight']]],
                             axis=1)

    for cond in featset256.keys():
        plt.figure(figsize = [12,12])
        ax = sns.scatterplot(x='PC_1',
                        y='PC_2',
                        data=PC_plotting[PC_plotting['bluelight']==cond],
                        hue='worm_gene',
                        hue_order=neuro_genes,
                        palette=cmap,
                        linewidth=0,
                        alpha=0.8)
        plt.title(cond)
        ax.axes.set_xlim(-40,40)
        ax.axes.set_ylim(-40,40)
        # plt.savefig(SAVETO / 'PC1_PC3_{}_allcondstogether.png'.format(cond),
        #             dpi=400)

    # groupby worm gene to see the trajectory through PC space
    PC_plotting_grouped = PC_plotting.groupby(['worm_gene',
                                               'bluelight']).mean().reset_index()
    PC_plotting_grouped['stimuli_order'] = PC_plotting_grouped['bluelight'].map(stimuli_order)
    PC_plotting_grouped.sort_values(by=['worm_gene', 'stimuli_order'],
                                    ascending=True,
                                    inplace=True)

    # nice figures
    plt.figure(figsize = [12,12])
    sns.lineplot(x='PC_1',
                    y='PC_2',
                    data=PC_plotting_grouped,
                    hue='worm_gene',
                    hue_order=neuro_genes,
                    palette=cmap,
                    alpha=0.8,
                    legend=False,
                    sort=False)
    sns.scatterplot(x='PC_1',
                    y='PC_2',
                    data=PC_plotting_grouped,
                    hue='worm_gene',
                    style='bluelight',
                    style_order=featsets.keys(),
                    hue_order=neuro_genes,
                    palette=cmap,
                    linewidth=0,
                    s=70)  
    plt.xlim(-10,20)
    plt.ylim(-10,10) 
    plt.xlabel('PC_1 ({}%)'.format(np.round(pca.explained_variance_ratio_[0]*100,2)))
    plt.ylabel('PC_2 ({}%)'.format(np.round(pca.explained_variance_ratio_[1]*100,2)))                                 
    plt.tight_layout()
    plt.savefig(SAVETO / 'PC1PC2_256feats_trajectory_space.png', dpi=400)
    
    plt.figure(figsize = [12,12])
    sns.lineplot(x='PC_3',
                    y='PC_4',
                    data=PC_plotting_grouped,
                    hue='worm_gene',
                    hue_order=neuro_genes,
                    palette=cmap,
                    alpha=0.8,
                    legend=False,
                    sort=False)
    sns.scatterplot(x='PC_3',
                    y='PC_4',
                    data=PC_plotting_grouped,
                    hue='worm_gene',
                    style='bluelight',
                    style_order=featsets.keys(),
                    hue_order=neuro_genes,
                    palette=cmap,
                    linewidth=0,
                    linestyle='-',
                    s=70)  
    plt.xlim(-10,20)
    plt.ylim(-10,10) 
    plt.xlabel('PC_3 ({}%)'.format(np.round(pca.explained_variance_ratio_[2]*100,2)))
    plt.ylabel('PC_4 ({}%)'.format(np.round(pca.explained_variance_ratio_[3]*100,2)))                                 
    plt.tight_layout()
    plt.savefig(SAVETO / 'PC3PC4_256feats_trajectory_space.png', dpi=400)
    
    # %% Do LDA on longform space
    
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

    y_classes = ['{}_{}'.format(r.worm_gene, r.bluelight) for i,r in long_featmat.iterrows()]
    
    clf = LinearDiscriminantAnalysis()
    X3 = clf.fit_transform(long_featmatZ.loc[:,feat256],
                           long_featmat['worm_gene'])
    
    lda_plotting = pd.DataFrame(data=X3,
                                columns=['LD_{}'.format(l) for l in range(1,X3.shape[1]+1)])
    lda_plotting = pd.concat([lda_plotting,
                              long_featmat[['worm_gene',
                                            'bluelight']]],
                             axis=1)
    
    for cond in featset256.keys():
        # plt.figure(figsize = [12,12])
        # ax = sns.scatterplot(x='LD_1',
        #                 y='LD_2',
        #                 data=lda_plotting[lda_plotting['bluelight']==cond],
        #                 hue='worm_gene',
        #                 hue_order=neuro_genes,
        #                 palette=cmap,
        #                 linewidth=0,
        #                 alpha=0.8)
        # plt.title(cond)
        # ax.axes.set_xlim(-15,15)
        # ax.axes.set_ylim(-15,15)
        # plt.savefig(SAVETO / 'LD1_LD2_{}_feat256LDA.png'.format(cond))
        
        plt.figure(figsize = [12,12])
        ax = sns.scatterplot(x='LD_9',
                        y='LD_10',
                        data=lda_plotting[lda_plotting['bluelight']==cond],
                        hue='worm_gene',
                        hue_order=neuro_genes,
                        palette=cmap,
                        linewidth=0,
                        alpha=0.8)
        plt.title(cond)
        ax.axes.set_xlim(-15,15)
        ax.axes.set_ylim(-15,15)
        plt.savefig(SAVETO / 'LD3_LD4_{}_feat256LDA.png'.format(cond))

    lda_grouped = lda_plotting.groupby(['worm_gene', 'bluelight']).mean().reset_index()
    lda_grouped['stimuli_order'] = lda_grouped.bluelight.map(stimuli_order)
    lda_grouped = lda_grouped.sort_values(by=['worm_gene',
                                               'stimuli_order'],
                                          ascending=True)
    # nice figures
    plt.figure(figsize = [12,12])
    sns.lineplot(x='LD_1',
                    y='LD_2',
                    data=lda_grouped,
                    hue='worm_gene',
                    hue_order=neuro_genes,
                    palette=cmap,
                    alpha=0.8,
                    legend=False,sort=False)
    sns.scatterplot(x='LD_1',
                    y='LD_2',
                    data=lda_grouped,
                    hue='worm_gene',
                    style='bluelight',
                    style_order=featsets.keys(),
                    hue_order=neuro_genes,
                    palette=cmap,
                    linewidth=0,
                    s=70)  
    plt.xlim(-5,10)
    plt.ylim(-10,5)
    plt.tight_layout()
    plt.xlabel('LD_1 ({}%)'.format(np.round(clf.explained_variance_ratio_[0]*100,2)))
    plt.ylabel('LD_2 ({}%)'.format(np.round(clf.explained_variance_ratio_[1]*100,2)))
    plt.savefig(SAVETO / 'LD1_LD2_256feats_trajectory_space.png', dpi=400)
    
        # nice figures
    plt.figure(figsize = [12,12])
    sns.lineplot(x='LD_3',
                    y='LD_4',
                    data=lda_grouped,
                    hue='worm_gene',
                    hue_order=neuro_genes,
                    palette=cmap,
                    alpha=0.8,
                    legend=False,sort=False)
    sns.scatterplot(x='LD_3',
                    y='LD_4',
                    data=lda_grouped,
                    hue='worm_gene',
                    style='bluelight',
                    style_order=featsets.keys(),
                    hue_order=neuro_genes,
                    palette=cmap,
                    linewidth=0,
                    s=70)  
    plt.xlim(-5,10)
    plt.ylim(-10,5)
    plt.tight_layout()
    plt.xlabel('LD_3 ({}%)'.format(np.round(clf.explained_variance_ratio_[2]*100,2)))
    plt.ylabel('LD_4 ({}%)'.format(np.round(clf.explained_variance_ratio_[3]*100,2)))
    plt.savefig(SAVETO / 'LD3_LD4_256feats_trajectory_space.png', dpi=400)
