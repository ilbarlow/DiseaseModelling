#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 16:27:02 2021

@author: ibarlow

script to make figures associated with mental retardation disorders


"""
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from itertools import chain
# from tierpsytools.read_data.hydra_metadata import align_bluelight_conditions

sys.path.insert(0, '/Users/ibarlow/Documents/GitHub/pythonScripts/DiseaseModelling/hydra_screen/phenotype_summary')
from helper import (read_disease_data,
                    select_strains,
                    filter_features,
                    make_colormaps,
                    find_window,
                    # long_featmap,
                    BLUELIGHT_WINDOW_DICT,
                    DATES_TO_DROP,
                    STIMULI_ORDER)
from plotting_helper import  (plot_colormap,
                              plot_cmap_text,
                              feature_box_plots,
                              window_errorbar_plots,
                              CUSTOM_STYLE)

from ts_helper import (MODECOLNAMES,
                       plot_frac_by_mode)

ROOT_DIR = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen')
FEAT_FILE = list(ROOT_DIR.rglob('*filtered/features_summary_tierpsy_plate_20200930_125752.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = list(ROOT_DIR.rglob('*filtered/filenames_summary_tierpsy_plate_20200930_125752.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = list(ROOT_DIR.rglob('*wells_annotated_metadata.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')

RAW_DATA_DIR = Path('/Volumes/Ashur Pro2/DiseaseScreen')
WINDOW_FILES = RAW_DATA_DIR / 'Results' / 'window_summaries_newfilters'

LEARNING_FEATURES = ['curvature_head_norm_abs_50th_prestim',
                    'relative_to_neck_radial_velocity_head_tip_10th_prestim',
                    'motion_mode_backward_frequency_prestim',
                    'motion_mode_backward_duration_50th_prestim',
                    'speed_midbody_norm_50th_prestim',
                    'width_midbody_norm_50th_prestim',
                    ]

LEARNING_BLUELIGHT = ['motion_mode_backward_frequency_bluelight',
                      'motion_mode_backward_duration_50th_bluelight',
                     'angular_velocity_neck_abs_50th_bluelight'
                     ]

strain_list = ['glr-1',
               'glr-4',
               'add-1']

CONTROL_STRAIN = 'N2'

#%%

#%%
if __name__ == '__main__':

    #set style for all figures
    plt.style.use(CUSTOM_STYLE)
    sns.set_style('ticks')

    saveto = (ROOT_DIR / 'Figures' / 'learning')
    saveto.mkdir(exist_ok=True) 
    feat, meta = read_disease_data(FEAT_FILE,
                                   FNAME_FILE,
                                   METADATA_FILE,
                                   export_nan_worms=False)
    
    meta.worm_gene.replace({'C43B7.2':'figo-1'},
                           inplace=True)

    window_files = list(WINDOW_FILES.rglob('*_window_*'))
    window_feat_files = [f for f in window_files if 'features' in str(f)]
    window_feat_files.sort(key=find_window)
    window_fname_files = [f for f in window_files if 'filenames' in str(f)]
    window_fname_files.sort(key=find_window)

    assert (find_window(f[0]) == find_window(f[1]) for f in list(zip(window_feat_files, window_fname_files)))

    genes = [g for g in meta.worm_gene.unique() if g != CONTROL_STRAIN]
    # genes = [g.replace('C43B7.2', 'figo-1') for g in genes]
    genes = [g for g in genes if 'myo' not in g and 'unc-54' not in g]

    #%%
    feat_df, meta_df, idx, gene_list = select_strains(strain_list,
                                                    CONTROL_STRAIN,
                                                    feat_df=feat,
                                                    meta_df=meta)
    # filter features
    feat_df, meta_df, featsets = filter_features(feat_df,
                                                 meta_df)

    # strain_numbers.append(meta_df.groupby('worm_strain')['file_id_prestim'].describe()['count'])


    strain_lut, stim_lut, feat_lut = make_colormaps(gene_list,
                                                    featlist=featsets['all'],
                                                    idx=idx,
                                                    candidate_gene=strain_list
                                                    )

    # colorbars to map colors to strains
    plot_colormap(strain_lut)
    plt.savefig(saveto / 'strain_cmap.png')
    plot_cmap_text(strain_lut)
    plt.savefig(saveto / 'strain_cmap_text.png')

    plot_colormap(stim_lut, orientation='horizontal')
    plt.savefig(saveto / 'stim_cmap.png')
    plot_cmap_text(stim_lut)
    plt.savefig(saveto / 'stim_cmap_text.png')

    plt.close('all')
    
    #%% plot the features
    
    # and make nice plots of the selected figures
    for f in  LEARNING_FEATURES:
        feature_box_plots(f,
                          feat_df,
                          meta_df,
                          strain_lut,
                          show_raw_data=False,
                          add_stats=True)
        plt.savefig(saveto / '{}_boxplot.png'.format(f),
                    dpi=200)
        plt.close('all')
        
        #%% plot for the windows
    feat_windows = []
    meta_windows = []
    for c,f in enumerate(list(zip(window_feat_files, window_fname_files))):
        _feat, _meta = read_disease_data(f[0],
                                         f[1],
                                         METADATA_FILE,
                                         drop_nans=True)
        _meta['window'] = find_window(f[0])
        
        meta_windows.append(_meta)
        feat_windows.append(_feat)

    meta_windows = pd.concat(meta_windows)
    meta_windows.reset_index(drop=True,
                             inplace=True)
    meta_windows.worm_gene.replace({'C43B7.2':'figo-1'},
                                   inplace=True)
    
    feat_windows = pd.concat(feat_windows)
    feat_windows.reset_index(drop=True,
                             inplace=True)
    
    feat_windows_df, meta_windows_df, idx, gene_list = select_strains(strain_list,
                                                  CONTROL_STRAIN,
                                                  meta_windows,
                                                  feat_windows)

    # #only need the bluelight features
    bluelight_feats = [f for f in feat_windows_df.columns if 'bluelight' in f]
    feat_windows_df = feat_windows_df.loc[:,bluelight_feats]

    feat_windows_df, meta_windows_df, featsets = filter_features(feat_windows_df,
                                                           meta_windows_df)
    
    bluelight_feats = list(feat_windows_df.columns)


    (saveto / 'windows_features').mkdir(exist_ok=True)
    meta_windows_df['light'] = [x[1] for x in meta_windows_df['window'].map(BLUELIGHT_WINDOW_DICT)]
    meta_windows_df['window_sec'] = [x[0] for x in meta_windows_df['window'].map(BLUELIGHT_WINDOW_DICT)]
    meta_windows_df['stim_number'] = [x[2] for x in meta_windows_df['window'].map(BLUELIGHT_WINDOW_DICT)]

    stim_groups = meta_windows_df.groupby('stim_number').groups
    
    #%%
    for f in LEARNING_BLUELIGHT:
        (saveto / 'windows_features' / f).mkdir(exist_ok=True)
        window_errorbar_plots(f,
                              feat_windows_df,
                              meta_windows_df,
                              strain_lut,
                              plot_legend=True)
        plt.savefig(saveto / 'windows_features' / f / 'allwindows_{}'.format(f), dpi=200)
        plt.close('all')

        for stim,locs in stim_groups.items():
            window_errorbar_plots(f,
                                  feat_windows_df.loc[locs],
                                  meta_windows_df.loc[locs],
                                  strain_lut,
                                  plot_legend=True)
            plt.savefig(saveto / 'windows_features' / f / 'window{}_{}'.format(stim,f),
                        dpi=200)
            plt.close('all')
    

