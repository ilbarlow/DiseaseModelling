#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 09:15:35 2020

@author: ibarlow

Script that loops through all the strains and does first step of 
k sig features, bluelight features and ts plots 
"""

import pandas as pd
# import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from tierpsytools.analysis.significant_features import k_significant_feat   
import time   

sys.path.insert(0, '/Users/ibarlow/Documents/GitHub/pythonScripts/DiseaseModelling/hydra_screen/phenotype_summary')
from helper import (read_disease_data,
                    select_strains,
                    filter_features,
                    make_colormaps,
                    find_window,
                    strain_gene_dict,
                    BLUELIGHT_WINDOW_DICT,
                    DATES_TO_DROP)
from plotting_helper import  (plot_colormaps,
                              make_clustermaps,
                              clustered_barcodes,
                              feature_box_plots,
                              window_errorbar_plots)
from ts_helper import (align_bluelight_meta,
                       load_bluelight_timeseries_from_results,
                       make_feats_abs,
                       plot_strains_ts,
                       get_motion_modes,
                        # get_frac_motion_modes_with_ci,
                       plot_frac)

ANALYSIS_TYPE = ['all_stim', 'bluelight', 'timeseries']
exploratory = True
is_reload_timeseries_from_results = True
is_recalculate_frac_motion_modes = True

FEAT_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')

RAW_DATA_DIR = Path('/Volumes/Ashur Pro2/DiseaseScreen')
WINDOW_FILES = RAW_DATA_DIR / 'Results' /'window_summaries'

CONTROL_STRAIN = 'N2'
# CANDIDATE_GENE = 'dys-1'
# DISEASE_GENES = 
strains_done = [] #['mpz-1', 'snf-11', 'C43B7.2', 'unc-43'] add in myo and unc-54 here (.isin() method)

feat256_fname = Path('/Users/ibarlow/tierpsy-tools-python/tierpsytools/extras/feat_sets/tierpsy_256.csv')

selected_feats = ['curvature_std_head_norm_abs_50th_prestim',
                'curvature_mean_head_norm_abs_50th_bluelight',
                'minor_axis_50th_poststim']

#%%
if __name__ == '__main__':
    
    feat, meta = read_disease_data(FEAT_FILE,
                                   FNAME_FILE,
                                   METADATA_FILE,
                                   export_nan_worms=False)
    
    window_files = list(WINDOW_FILES.rglob('*_window_*'))   
    window_feat_files = [f for f in window_files if 'features' in str(f)]
    window_feat_files.sort(key=find_window)
    window_fname_files = [f for f in window_files if 'filenames' in str(f)]
    window_fname_files.sort(key=find_window)
    
    assert (find_window(f[0]) == find_window(f[1]) for f in list(zip(window_feat_files, window_fname_files)))
    
    genes = [g for g in meta.worm_gene.unique() if g != CONTROL_STRAIN]
    genes = [g for g in genes if 'myo' not in g and 'unc-54' not in g]
    
    genes = list(set(genes) - set(strains_done))
    
    for g in genes:
        
        if 'all_stim' in ANALYSIS_TYPE:
        
            CANDIDATE_GENE = g # all caps usually global so should set to lower
            print ('all stim plots for {}'.format(CANDIDATE_GENE))

            SAVETO = FEAT_FILE.parent.parent.parent / 'Figures' / 'paper_figures' / CANDIDATE_GENE
            SAVETO.mkdir(exist_ok=True)
        
            feat_df, meta_df, idx, gene_list = select_strains(CANDIDATE_GENE,
                                                              CONTROL_STRAIN,
                                                              feat_df=feat,
                                                              meta_df=meta)
            
            #TODO look at Eleni's preprocessing
            feat_df, meta_df, featsets = filter_features(feat_df,
                                                         meta_df)
            
            strain_lut, stim_lut, feat_lut = make_colormaps(gene_list,
                                                            idx,
                                                            CANDIDATE_GENE,
                                                            CONTROL_STRAIN,
                                                            featlist=featsets['all'])
            
            # colorbars to map colors to strains
            plot_colormaps(strain_lut, stim_lut, SAVETO)
            
            #%%  fillnans and normalise
            feat_nonan = feat_df.copy()
            feat_nonan.fillna(feat_nonan.mean(axis=0),
                              inplace=True)
        
            featZ = pd.DataFrame(data=stats.zscore(feat_nonan, axis=0),
                                 columns=featsets['all'],
                                 index=feat_nonan.index)
            
            assert featZ.isna().sum().sum() == 0
            #%% make a nice clustermap / heatmap
            (SAVETO / 'clustermaps').mkdir(exist_ok=True)
            
            #TODO should I use one of sklearns algorithms instead?
            clustered_features = make_clustermaps(featZ,
                                                  meta_df,
                                                  featsets,
                                                  strain_lut,
                                                  feat_lut,
                                                  saveto=SAVETO / 'clustermaps')
            plt.close('all')
                    
            #%% k significant features for each prestim, bluelight and poststim
            (SAVETO / 'ksig_feats').mkdir(exist_ok=True)
            sns.set_style('white')
            label_format = '{:.4f}'
            kfeats = {}
            for stim, fset in featsets.items():
                kfeats[stim], scores, support = k_significant_feat(
                    feat_nonan[fset],
                    meta_df.worm_gene,
                    k=100,
                    plot=False)
                
                for i in range(0,5):
                    fig, ax = plt.subplots(4, 5, sharex=True, figsize = [20,20])
                    for c, axis in enumerate(ax.flatten()):
                        counter=(20*i)-(20-c)
                        sns.boxplot(x=meta_df['worm_gene'],
                                    y=feat_df[kfeats[stim][counter]],
                                    palette=strain_lut.values(),
                                    ax=axis)
                        axis.set_ylabel(fontsize=8,
                                        ylabel=kfeats[stim][counter])
                        axis.set_yticklabels(labels=[label_format.format(x) for x in axis.get_yticks()],
                                             fontsize=6)
                        axis.set_xlabel('')
                    plt.tight_layout()
                    fig.fontsize=11
                    plt.savefig(SAVETO / 'ksig_feats' / '{}_{}_ksig_feats.png'.format(i*20, stim),
                                dpi=400)
                    plt.close('all')
            
            #%% pairwise statistics to find features that are different from N2
            if not exploratory:
                (SAVETO / 'heatmaps').mkdir(exist_ok=True)
                from tierpsytools.analysis.paired_stats_tests import paired_stats_tests
                
                pVals, bhP_values, group_classes = paired_stats_tests(feat_df,
                                                                       meta_df.loc[:, 'worm_gene'],
                                                                       control_group='N2')
                bhP_values['worm_gene'] = group_classes
                bhP_values.rename(mapper={0:'p<0.05'},
                                  inplace=True)
                # bhP_values.set_index('worm_gene', inplace=True)
                
                #%%
                # nice figure with single barcode for each strains, and asterisks of signficicantly different features
                
                clustered_barcodes(clustered_features, selected_feats,
                                   featZ,
                                   meta_df,
                                   bhP_values,
                                   SAVETO / 'heatmaps')
            
                # and make nice plots of the selected figures
                for f in selected_feats:
                    feature_box_plots(f, feat, meta, bhP_values, strain_lut)
                    plt.savefig(SAVETO / '{}_boxplot.png'.format(f), dpi=200)

        # %% bluelight
        if 'bluelight' in ANALYSIS_TYPE:
            print ('all window_plots for {}'.format(CANDIDATE_GENE))
            
            feat_windows = []
            meta_windows = []
            for c,f in enumerate(list(zip(window_feat_files, window_fname_files))):
                _feat, _meta = read_disease_data(f[0],
                                                 f[1],
                                                 METADATA_FILE,
                                                 drop_nans=True)
                _meta['window'] = find_window(f[0])
                # feat_df, meta_df = drop_nan_worms(feat_df, meta_df, saveto=None)
                _feat, _meta, idx, gene_list = select_strains(CANDIDATE_GENE,
                                                              CONTROL_STRAIN,
                                                              _meta,
                                                              _feat)
                meta_windows.append(_meta)
                feat_windows.append(_feat)
             
            meta_windows = pd.concat(meta_windows)
            meta_windows.reset_index(drop=True,
                                     inplace=True)
            feat_windows = pd.concat(feat_windows)
            feat_windows.reset_index(drop=True,
                                     inplace=True)
           
            #only need the bluelight features
            bluelight_feats = [f for f in feat_windows.columns if 'bluelight' in f]
            feat_windows = feat_windows[bluelight_feats]
    
            feat_windows, meta_windows, featsets = filter_features(feat_windows,
                                                                   meta_windows)
            
            bluelight_feats = [f for f in feat_windows.columns if 'bluelight' in f]
            strain_lut, stim_lut, feat_lut = make_colormaps(gene_list,
                                                            idx,
                                                            CANDIDATE_GENE,
                                                            CONTROL_STRAIN,
                                                            featlist=bluelight_feats)
            
            #%% fillnans and normalize
            feat_nonan = feat_windows.copy()
            feat_nonan.fillna(feat_nonan.mean(axis=0),
                              inplace=True)
        
            featZ = pd.DataFrame(data=stats.zscore(feat_nonan[bluelight_feats], axis=0),
                                 columns=bluelight_feats,
                                 index=feat_nonan.index)
            
            assert featZ.isna().sum().sum() == 0
            
            #%%      
            # find the k sig feats that differentiate between prelight and bluelight
            # for N2 vs cat-2
            
            (SAVETO / 'windows_features').mkdir(exist_ok=True)
            meta_windows['light'] = meta_windows['window'].map(BLUELIGHT_WINDOW_DICT)
            y_classes = ['{}, {}'.format(r.worm_gene, r.light) for i,r in meta_windows.iterrows()]
            
            kfeats, scores, support = k_significant_feat(
                    feat_nonan,
                    y_classes,
                    k=100,
                    plot=False,
                    score_func='f_classif')
            
            for f in kfeats[:50]:
                window_errorbar_plots(f, 
                                      feat_windows,
                                      meta_windows,
                                      strain_lut,
                                      window_order=BLUELIGHT_WINDOW_DICT.keys())   
                plt.savefig(SAVETO / 'windows_features' / '{}'.format(f))
                plt.close('all')
                
        #%% make nice ts plots
        
        if 'timeseries' in ANALYSIS_TYPE:
            print ('timeseries plots for {}'.format(CANDIDATE_GENE))

            timeseries_fname = RAW_DATA_DIR/ 'Results' / '{}_timeseries.hdf5'.format(CANDIDATE_GENE)
            
            meta_ts = pd.read_csv(METADATA_FILE, index_col=None)  
            assert meta_ts.worm_strain.unique().shape[0] == meta_ts.worm_gene.unique().shape[0]
            
            meta_ts.loc[:,'date_yyyymmdd'] = meta_ts['date_yyyymmdd'].apply(lambda x: str(int(x)))     
            #drop nan wells
            meta_ts.dropna(axis=0,
                        subset=['worm_gene'],
                        inplace=True)   
            # remove data from dates to exclude
            good_date = meta_ts.query('@DATES_TO_DROP not in imaging_date_yyyymmdd').index
            # bad wells
            good_wells_from_gui = meta_ts.query('is_bad_well == False').index
            meta_ts = meta_ts.loc[good_wells_from_gui & good_date,:]
            
            # only select strains of interest
            meta_ts, idx, gene_list = select_strains(CANDIDATE_GENE,
                                                  CONTROL_STRAIN,
                                                  meta_df=meta_ts)        
            strain_lut, stim_lut = make_colormaps(gene_list,
                                                    idx,
                                                    CANDIDATE_GENE,
                                                    CONTROL_STRAIN,
                                                    featlist=[])
            #strain to gene dictionary
            strain_dict = strain_gene_dict(meta_ts)
            gene_dict = {v:k for k,v in strain_dict.items()}
    
            meta_ts = align_bluelight_meta(meta_ts)
            
            if is_reload_timeseries_from_results:
                # this uses tierpytools under the hood
                timeseries_df, hires_df  = load_bluelight_timeseries_from_results(
                                    meta_ts,
                                    RAW_DATA_DIR / 'Results')
                                    # save to disk
                try:
                    timeseries_df.to_hdf(timeseries_fname, 'timeseries_df', format='table')
                    hires_df.to_hdf(timeseries_fname, 'hires_df', format='table')
                except Exception:
                    print ('error creating {} HDF5 file'.format(CANDIDATE_GENE))
            else:  # from disk, then add columns
                # dataframe from the saved file
                timeseries_df = pd.read_hdf(timeseries_fname, 'timeseries_df')
                hires_df = pd.read_hdf(timeseries_fname, 'hires_df')
                
            #%% add in information about the replicates
            date_to_repl = pd.DataFrame({'date_yyyymmdd': [
                                                          '20200730',
                                                         '20200801',
                                                         '20200806',
                                                         '20200808',
                                                         '20200811',
                                                         '20200813',
                                                        '20200902',
                                                        '20200903',
                                                        '20200904'
                                                           ],
                                     'replicate': [1, 1, 2, 2, 3, 3, 1, 2, 3]})
            timeseries_df = pd.merge(timeseries_df, date_to_repl,
                                     how='left',
                                     on='date_yyyymmdd')
            
            timeseries_df['worm_strain'] = timeseries_df['worm_gene'].map(gene_dict)
            hires_df['worm_strain'] = hires_df['worm_gene'].map(gene_dict)
    
            #make d/v signed features absolute as in hydra d/v is not assigned 
            timeseries_df = make_feats_abs(timeseries_df)
                
       
            # %%% hand-picked features from the downsampled dataframe
            
            plt.close('all')
            (SAVETO / 'ts_plots').mkdir(exist_ok=True)
            feats_toplot = ['speed',
                            'abs_speed',
                            'angular_velocity',
                            'abs_angular_velocity',
                            'relative_to_body_speed_midbody',
                            'abs_relative_to_body_speed_midbody',
                            'abs_relative_to_neck_angular_velocity_head_tip',
                            'speed_tail_base',
                            'length',
                            'major_axis',
                            'd_speed',
                            'head_tail_distance',
                            'abs_angular_velocity_neck',
                            'abs_angular_velocity_head_base',
                            'abs_angular_velocity_hips',
                            'abs_angular_velocity_tail_base',
                            'abs_angular_velocity_midbody',
                            'abs_angular_velocity_head_tip',
                            'abs_angular_velocity_tail_tip']
            
            plot_strains_ts(timeseries_df, 
                            strain_lut,
                            CONTROL_STRAIN,
                            feats_toplot,
                            SAVETO / 'ts_plots')
            plt.close('all')
            #%% motion modes
            # get motion_mode stats
            tic = time.time()    
            if is_recalculate_frac_motion_modes:
                motion_modes, frac_motion_modes_with_ci = get_motion_modes(hires_df,
                                                                            saveto=timeseries_fname
                                                                           )
            else:
                frac_motion_modes_with_ci = pd.read_hdf(timeseries_fname,
                                                       'frac_motion_mode_with_ci')
          
            frac_motion_modes_with_ci['worm_strain'] = frac_motion_modes_with_ci['worm_gene'].map(gene_dict)
    
            fps = 25
            frac_motion_modes_with_ci = frac_motion_modes_with_ci.reset_index()
            frac_motion_modes_with_ci['time_s'] = (frac_motion_modes_with_ci['timestamp']
                                                  / fps)
            print('Time elapsed: {}s'.format(time.time()-tic))
            
            # %% make some nicer plotes with confindence intervals using luigi's bootstrapping
            for ii, (strain, df_g) in enumerate(frac_motion_modes_with_ci.groupby('worm_gene')):
                plot_frac(df_g, strain, strain_lut)
                plt.savefig(SAVETO / '{}_motion_modes.png'.format(strain), dpi=200)
    
