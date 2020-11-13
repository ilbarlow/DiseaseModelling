#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 16:59:05 2020

@author: ibarlow

Analyse the timeseries windows in the bluelight and prestim and poststim to 
evaluate:
    
    Time windows (seconds):
        0 - 50:60
        1 - 65:75 BLUELIGHT
        2 - 75:85
        3 - 150:160
        4 - 165:175 BLUELIGHT
        5 - 175:185 
        6 - 250:260
        7 - 265:275 BLUELIGHT
        8 - 275:285
    
    1. differences that can be seen in the summary statistics between the 
    strains that look like they 'freeze' after the bluelight
    
    2. see if differences between strains can be observed in only a 10sec 
    time window
    
    3. Comparison of each strain between timewindows - like a dose response
    
    

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

PROJECT_DIR = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/Results/window_summaries')

METADATA_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')
SAVETO = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/Figures/bluelight_windows') / subgroup
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
    
# vars
CONTROL_STRAIN = 'N2'
DATES_TO_DROP = '20200626'
BAD_FEAT_THRESH = 3 # 3 standard deviations away from the mean
BAD_FEAT_FILTER = 0.1 # 10% threshold for removing bad features
BAD_WELL_FILTER = 0.3 # 30% threshold for bad well
# NAN_CUTOFF = 7000

bluelight_window_dict = {0: 'prelight',
                         1: 'bluelight',
                         2: 'postlight',
                         3: 'prelight',
                         4: 'bluelight',
                         5: 'postlight',
                         6: 'prelight',
                         7: 'bluelight',
                         8: 'postlight'}

stimuli_order = {'prestim':1,
                 'bluelight':2,
                 'poststim':3}
#%%
def find_window(fname):
    import re
    window_regex = r"(?<=_window_)\d{0,9}"
    window = int(re.search(window_regex, str(fname))[0])
    return window

# %%
if __name__ =='__main__':
    
    meta = pd.read_csv(METADATA_FILE)
    
    summary_files = list(PROJECT_DIR.rglob('*_window_*'))   
    feat_files = [f for f in summary_files if 'features' in str(f)]
    feat_files.sort(key=find_window)
    fname_files = [f for f in summary_files if 'filenames' in str(f)]
    fname_files.sort(key=find_window)

    assert (find_window(f[0]) == find_window(f[1]) for f in list(zip(feat_files, fname_files)))

    feat_df = []
    meta_df = []
    for c,f in enumerate(list(zip(feat_files, fname_files))):
        _feat = pd.read_csv(f[0],
                            comment='#')
        _fname = pd.read_csv(f[1],
                             comment='#')
    
        _feat, _meta = read_hydra_metadata(_feat, _fname, meta)
    
        _feat, _meta = align_bluelight_conditions(_feat,
                                                  _meta,
                                                  how='inner')
        _meta['window'] = find_window(f[0])
        meta_df.append(_meta)
        feat_df.append(_feat)
    
    assert pd.concat(meta_df).shape[0] == pd.concat(feat_df).shape[0]

    meta = pd.concat(meta_df)
    meta.reset_index(drop=True, inplace=True)
    feat = pd.concat(feat_df)
    feat.reset_index(drop=True, inplace=True)
    
    #%% wells to check
    nan_worms = meta[meta.worm_gene.isna()][['featuresN_filename',
                                             'well_name',
                                             'imaging_plate_id',
                                             'instrument_name',
                                             'imaging_date_yyyymmdd']]
    # nan_worms.to_csv(SAVETO / 'nan_worms.csv',
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
        N2_locs = meta.groupby('window').apply(lambda x: x.query('@CONTROL_STRAIN in worm_gene').sample(300)).index
        N2_locs = [i[1] for i in N2_locs]
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

    # remove wells with too many nans
    feat = feat.loc[feat.isna().sum(axis=1)/feat.shape[0]<BAD_WELL_FILTER, :]
    meta = meta.loc[feat.index,:]

    # remove features with too many nans
    feat = feat.loc[:, feat.isna().sum(axis=0)/feat.shape[0]<BAD_FEAT_FILTER]

    # remove features with std=0
    feat = feat.loc[:, feat.std(axis=0) > 1e-6 ]

     # feature sets
     # abs features no longer in tierpsy
    pathcurvature_feats = [x for x in feat.columns if 'path_curvature' in x]
    for f in featsets.keys():
        featsets[f] = [x for x in feat.columns if f in x]
        featsets[f] = list(set(featsets[f]) - set(pathcurvature_feats))
        
    meta['window_light'] = meta.window.map(bluelight_window_dict)
        
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
    
        window_groups = meta.groupby('window').groups
        for w in window_groups: 
            neuro_feats = {}
            for fset in featsets.keys():
                neuro_feats[fset],scores,support =  k_significant_feat(
                    feat_nonan.loc[window_groups[w],
                                   featsets[fset]],
                    meta.loc[window_groups[w], 'worm_gene'],
                    k=10,
                    plot=False)
            
            # save some plots
            _window_type = str(meta.loc[window_groups[w], 'window_light'].unique()[0])
            (SAVETO / 'k_sig_feats' / 'window{}_{}'.format(w, _window_type)).mkdir(exist_ok=True)
            for f in neuro_feats:
                (SAVETO / 'k_sig_feats' / 'window{}_{}'.format(w, _window_type) / f).mkdir(exist_ok=True)
                for s in neuro_feats[f]:
                    plt.figure()
                    sns.boxplot(x=feat_nonan.loc[window_groups[w], s],
                                y=meta.loc[window_groups[w], 'worm_gene'],
                                order=neuro_genes,
                                orient='h',
                                showfliers=False,
                                palette = cmap)
                    plt.savefig(SAVETO / 'k_sig_feats' / 'window{}_{}'.format(w, _window_type) / f / '{}.png'.format(s),
                                dpi=200)
                    plt.close('all')

# %% make a long feature matrix
    long_featmat = []
    for f in featset256:
        _featmat = pd.DataFrame(data=feat.loc[:,featsets[f]].values,
                                columns=['_'.join(s.split('_')[:-1])
                                           for s in featsets[f]],
                                index=feat.index)
        _featmat['bluelight'] = f
        _featmat = pd.concat([_featmat,
                              meta.loc[:,['worm_gene', 'window']]],
                             axis=1)
        long_featmat.append(_featmat)

    long_featmat = pd.concat(long_featmat,
                             axis=0)
    long_featmat.reset_index(drop=True,
                             inplace=True)
    
    long_featmat['window_light'] = long_featmat.window.map(bluelight_window_dict)
    long_featlist = list(set(long_featmat.columns) - set(['worm_gene', 'window', 'bluelight', 'window_light']))
    long_featmat.loc[:,long_featlist] = long_featmat.fillna(long_featmat[long_featlist].mean(axis=0))

    # long_featmatZ = pd.DataFrame(stats.zscore(long_featmat[feat256]),
    #                              columns=feat256)

     # %% k significant features for each strain across each time window only on the 
    # bluelight videos
    import matplotlib.patches as patches
    #only look at the 'bluelight' videos for now
    (SAVETO / '{}_window_comparisons'.format('bluelight')).mkdir(exist_ok=True)
    bluelight_features = long_featmat.query('"bluelight" in bluelight')
    N2_locs = bluelight_features.query('@CONTROL_STRAIN in worm_gene').index
    n2_index=-1

    for c, g in enumerate(neuro_genes):
        
        gene_locs = bluelight_features.query('@g in worm_gene').index
        
        window_feats, scores, support = k_significant_feat(bluelight_features.loc[gene_locs, long_featlist],
                                                           bluelight_features.loc[gene_locs, 'window_light'].values,
                                                           k=20,
                                                           plot=False)

        (SAVETO / '{}_window_comparisons'.format('bluelight') / g).mkdir(exist_ok=True)
        for w in window_feats:
            plt.figure()
            ax = sns.pointplot(x='window',
                                y=w,
                                data=bluelight_features.loc[gene_locs.append(N2_locs),:],
                                order=bluelight_window_dict.keys(),
                                hue='worm_gene',
                                showfliers=False,
                                linestyle=['-'],
                                ci=95,
                                palette=[cmap[c],
                                         cmap[n2_index]])
            y_min = ax.axes.get_ylim()[0]
            y_max = ax.axes.get_ylim()[1]
            rects = (patches.Rectangle((0.5, y_min), 1, y_max, facecolor='b', alpha=0.3),
                     patches.Rectangle((3.5, y_min), 1, y_max, facecolor='b', alpha=0.3),
                     patches.Rectangle((6.5, y_min), 1, y_max, facecolor='b', alpha=0.3))

            [ax.add_patch(r) for r in rects]
            plt.savefig(SAVETO / '{}_window_comparisons'.format('bluelight') / g / '{}.png'.format(w))
            plt.close('all')
 
#%% plot motion mode fractions for each window - follow up from the Luigi's ts analysis
    motion_feats = [f for f in long_featlist if 'motion_mode' in f]
    
    (SAVETO / '{}_window_comparisons'.format('bluelight') / 'motion_modes').mkdir(exist_ok=True)

    for c, g in enumerate(neuro_genes):
        gene_locs = bluelight_features.query('@g in worm_gene').index
        (SAVETO / '{}_window_comparisons'.format('bluelight') / 'motion_modes' / g).mkdir(exist_ok=True)

        for m in motion_feats:
            plt.figure()
            ax = sns.pointplot(x='window',
                              y=m,
                              data=bluelight_features.loc[gene_locs.append(N2_locs)],
                              order=bluelight_window_dict.keys(),
                                hue='worm_gene',
                                showfliers=False,
                                linestyle=['-'],
                                ci=95,
                                palette=
                                [cmap[c],
                                         cmap[n2_index]])
            y_min = ax.axes.get_ylim()[0]
            y_max = ax.axes.get_ylim()[1]
            rects = (patches.Rectangle((0.5, y_min), 1, y_max, facecolor='b', alpha=0.3),
                     patches.Rectangle((3.5, y_min), 1, y_max, facecolor='b', alpha=0.3),
                     patches.Rectangle((6.5, y_min), 1, y_max, facecolor='b', alpha=0.3))

            [ax.add_patch(r) for r in rects]
            plt.savefig(SAVETO / '{}_window_comparisons'.format('bluelight') / 'motion_modes' / g / '{}.png'.format(m))
            plt.close('all')

#%% plot only w_paused features

    paused_feats = [f for f in long_featlist if 'w_paused' in f]

    (SAVETO / 'pairwise_stats_test').mkdir(exist_ok=True)
    (SAVETO / 'pairwise_stats_test' / 'paused_feats').mkdir(exist_ok=True)
    
    candidate_genes = ['unc-77',
                       'unc-80',
                       'nca-2']
    
    bluelight_fainters = long_featmat.query('"bluelight" in bluelight and @candidate_genes in worm_gene or @CONTROL_STRAIN in worm_gene')
    
    fainter_df = []
    for c, g in enumerate(candidate_genes):
        
        _test = [g, 'N2']
        gene_locs = bluelight_fainters.query('@_test in worm_gene').index
        
        for w in bluelight_window_dict:
            window_locs = bluelight_fainters.loc[gene_locs].query('window == @w').index
            window_feats, scores, support = k_significant_feat(bluelight_fainters.loc[window_locs, paused_feats],
                                                               bluelight_fainters.loc[window_locs, 'worm_gene'].values,
                                                               k=20,
                                                               plot=False)
            
            print(g, w, window_feats[:3])
            fainter_df.append(pd.DataFrame({'worm_gene':g,
                                    'window':w,
                                    'top20':list(window_feats),
                                   }))
            del window_feats, scores, support


    fainter_df = pd.concat(fainter_df)
    fainter_df['window_light'] = fainter_df.window.map(bluelight_window_dict)
    fainter_df.reset_index(drop=True,
                           inplace=True)
    
    # for c,g in candidate_genes:
    #     (SAVETO / 'pairwise_stats_test' / 'paused_feats' / g).mkdir(exist_ok=True)        
    for p in ['prelight', 'bluelight', 'postlight']:
        paused_feats = list(fainter_df.query('@p in window_light').drop_duplicates(subset='top20')['top20'])

        for w in paused_feats:
            plt.figure()
            ax = sns.pointplot(x='window',
                                y=w,
                                data=bluelight_fainters,
                                order=bluelight_window_dict.keys(),
                                hue='worm_gene',
                                showfliers=False,
                                linestyle=['-'],
                                ci=95,
                                alpha=0.8,
                                palette='pastel'
                                # palette=[cmap[n2_index]],
                                #         cmap[c]
                                          )
            y_min = ax.axes.get_ylim()[0]
            y_max = ax.axes.get_ylim()[1]
            rects = (patches.Rectangle((0.5, y_min), 1, y_max, facecolor='b', alpha=0.3),
                     patches.Rectangle((3.5, y_min), 1, y_max, facecolor='b', alpha=0.3),
                     patches.Rectangle((6.5, y_min), 1, y_max, facecolor='b', alpha=0.3))

            [ax.add_patch(r) for r in rects]
            plt.savefig(SAVETO / 'pairwise_stats_test' / 'paused_feats'/ '{}_{}.png'.format(p, w))
            plt.close('all')
 
    
    # if pairwise_stats:
    #     (SAVETO / 'pairwise_stats_test').mkdir(exist_ok=True)

    #     from tierpsytools.analysis.paired_stats_tests import paired_stats_tests
    #     for w in window_groups:
    #         _window_type = str(meta.loc[window_groups[w], 'window_light'].unique()[0])
    #         pVals, bhP_values, group_classes = paired_stats_tests(feat.loc[window_groups[w],:],
    #                                                                meta.loc[window_groups[w], 'worm_gene'],
    #                                                                control_group='N2')
    #         bhP_values['worm_gene'] = group_classes
    #         bhP_values.set_index('worm_gene', inplace=True)
            
    #         if subgroup == 'neuromodels':
    #             n2_index=-1
    #         (SAVETO / 'pairwise_stats_test' / 'window{}_{}'.format(w, _window_type)).mkdir(exist_ok=True)
    #         for c,g in enumerate(neuro_genes):
    #             if g == 'N2':
    #                 continue
    #             else:
    #                 genes_to_plot = [g, CONTROL_STRAIN]
        
    #                 _gene = bhP_values.columns[bhP_values.loc[g,:].notna()]
    #                 _sort = np.argsort(bhP_values.loc[g, _gene])
    #                 _top = _gene[_sort]
        
    #                 (SAVETO / 'pairwise_stats_test' / 'window{}_{}'.format(w, _window_type) / g).mkdir(exist_ok=True)
    #                 (SAVETO / 'pairwise_stats_test' / 'window{}_{}'.format(w, _window_type)/ g / 'top25').mkdir(exist_ok=True)
    #                 for f in _top[:25]:
    #                     plt.figure(figsize=[10,5])
    #                     ax1 = plt.subplot(1,2, 1)
    #                     ax1 = sns.boxplot(y=meta.loc[:, 'worm_gene'],
    #                                 x=feat.loc[:, f],
    #                                 orient='h',
    #                                 order=neuro_genes,
    #                                 showfliers=False,
    #                                 palette=cmap)
    #                     ax2 = plt.subplot(1,2,2)
    #                     ax2 = sns.boxplot(y=meta.loc[meta.query('@genes_to_plot in worm_gene').index,
    #                                                   'worm_gene'],
    #                                       x=feat.loc[meta.query('@genes_to_plot in worm_gene').index,
    #                                                   f],
    #                                       orient='h',
    #                                       showfliers=False,
    #                                       palette = [cmap[c],
    #                                                   cmap[n2_index]],
    #                                       order = genes_to_plot
    #                                       )
    #                     plt.savefig(SAVETO / 'pairwise_stats_test'/ 'window{}_{}'.format(w, _window_type) / g / 'top25' / \
    #                                 '{}.png'.format(f))
    #                     plt.close('all')
        
    #                 #plot bluelight comparisons
    #                 (SAVETO / 'pairwise_stats_test' / 'window{}_{}'.format(w, _window_type) / g / 'bluelight_comparisons').mkdir(exist_ok=True)
    #                 _bluelight = ['_'.join(f.split('_')[:-1]) for f in _top[:40] if 'bluelight' in f]
                    
    #                 id_to_merge = ['worm_gene',
    #                                 'well_name',
    #                                 'imaging_plate_id',
    #                                 'imaging_date_yyyymmdd']
    #                 for f in _bluelight:
    #                     flist = [f+'_{}'.format(s) for s in featsets.keys()]
    #                     plot_df = pd.concat([feat.loc[:, flist], meta],
    #                                         axis=1)
    #                     to_plot = plot_df.melt(id_vars=id_to_merge,
    #                                             value_vars=flist,
    #                                             value_name=f)
    #                     to_plot['bluelight'] = to_plot.variable.apply(lambda x: x.split('_')[-1])
        
        
    #                     plt.figure()
    #                     sns.pointplot(data=to_plot.query('@genes_to_plot in worm_gene'),
    #                                 x='bluelight',
    #                                 y=f,
    #                                 hue='worm_gene',
    #                                 ci=95,
    #                                 palette='cubehelix')
    #                     plt.savefig(SAVETO / 'pairwise_stats_test' / 'window{}_{}'.format(w, _window_type) / g / 'bluelight_comparisons' /\
    #                                 '{}_bluelightcomparisons.png'.format(f))
    #                     plt.close('all')



