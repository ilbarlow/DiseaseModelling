#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 14:17:23 2021

@author: ibarlow

look for shrinker defects

"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
from scipy import signal
# from tierpsytools.preprocessing.preprocess_features import impute_nan_inf

sys.path.insert(0, '/Users/ibarlow/Documents/GitHub/pythonScripts/DiseaseModelling/hydra_screen/phenotype_summary')
from helper import (
                    select_strains,
                    make_colormaps,
                    strain_gene_dict,
                    BLUELIGHT_WINDOW_DICT,
                    DATES_TO_DROP,
                    STIMULI_ORDER)

from plotting_helper import CUSTOM_STYLE
from luigi_helper import plot_stimuli

ROOT_DIR = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen')
FEAT_FILE = list(ROOT_DIR.rglob('*filtered/features_summary_tierpsy_plate_20200930_125752.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = list(ROOT_DIR.rglob('*filtered/filenames_summary_tierpsy_plate_20200930_125752.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = list(ROOT_DIR.rglob('*wells_annotated_metadata.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')

RAW_DATA_DIR = Path('/Volumes/Ashur Pro2/DiseaseScreen')

saveto = ROOT_DIR / 'bluelight_ts_figures'
saveto.mkdir(exist_ok=True)

#vars and filters
CONTROL_STRAIN = 'N2'
n_subsample = 10 #10 wells for each strain
seed=20210205
feature = 'length' #'d_relative_to_neck_radial_velocity_head_tip'#'length' 
frameRate = 25
minTrajDuration = 100 # mininum trajectory in seconds
smoothWindow = 1 # smooth feature over the window, in seconds

groupby_vars = ['imaging_plate_id', 'date_yyyymmdd', 'well_name']

#%%
if __name__ == "__main__":

    plt.style.use(CUSTOM_STYLE)
    sns.set_style('ticks')

    meta_ts = pd.read_csv(METADATA_FILE,
                          index_col=None)
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

    meta_ts.replace({'C43B7.2':'figo-1'},
                    inplace=True)

    genes = [g for g in meta_ts.worm_gene.unique() if g != CONTROL_STRAIN]
    genes = [g for g in genes if 'myo' not in g and 'unc-54' not in g]

    #%%

    for g in genes:
        print ('making ts stimuli plots for {}'.format(g))
        candidate_gene = g
        timeseries_fname = RAW_DATA_DIR / 'Results' / '{}_timeseries.hdf5'.format(candidate_gene)
        (saveto / candidate_gene).mkdir(exist_ok=True)

        # only select strains of interest
        meta, idx, gene_list = select_strains(candidate_gene,
                                              CONTROL_STRAIN,
                                              meta_df=meta_ts)
        strain_lut, stim_lut = make_colormaps(gene_list,
                                              [],
                                                idx,
                                                candidate_gene)
        
        if g == 'figo-1':
            strain_lut['C43B7.2'] = strain_lut['figo-1']

        timeseries_df = pd.read_hdf(timeseries_fname, 'timeseries_df')
        # hires_df = pd.read_hdf(timeseries_fname, 'hires_df')
        
        # test_feats = [f for f in timeseries_df.columns if 'radial_velocity' in f]

        # t = ['d_relative_to_neck_radial_velocity_head_tip',
        #      'relative_to_neck_radial_velocity_head_tip']
             
        for k in strain_lut.keys():

            _ts_df = timeseries_df.query('@k in worm_gene')
            
            if not _ts_df.empty:
                _ts_df = _ts_df.groupby(groupby_vars).filter(lambda x: x.shape[0]>0)
    
                random_groups = _ts_df.groupby(groupby_vars).apply(len).sample(5,
                                                                               random_state = 20210205).index
    
                for well in random_groups:
    
                    _widx = _ts_df.groupby(groupby_vars).get_group(well).groupby('worm_index').filter(lambda x: x.shape[0] > minTrajDuration)
    
                    _worm_grouped = _widx.groupby('worm_index')
    
                    for worm in _worm_grouped.groups.keys():
    
                        # for t in test_feats:
                        y = _worm_grouped.get_group(worm)[feature].rolling(smoothWindow).median()
    
                        x = _worm_grouped.get_group(worm)['time_binned_s']
    
                        fig, ax = plt.subplots(figsize=[8,5])
    
                        plt.plot(x,y,
                                 color=strain_lut[k])
                        plot_stimuli(ax, units='s')
                        # plt.ylabel(t)
    
                        ax.set_xticks(list(range(0,int(round(x.max()))+1,30)))
                        plt.ylabel('length (microns)')
                        plt.xlabel('time (s)')
                        plt.tight_layout()
                        plt.savefig(saveto / candidate_gene / '{}_{}_{}_bluelight_{}_ts.png'.format(k, '_'.join(well), worm, feature))
                        plt.close('all')
