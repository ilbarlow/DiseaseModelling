#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 10:21:21 2021

@author: ibarlow

Have a look at the trajectories to see if can detect defecation cycle of 
N2s and mutants

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


ROOT_DIR = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen')
FEAT_FILE = list(ROOT_DIR.rglob('*filtered/features_summary_tierpsy_plate_20200930_125752.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/features_summary_tierpsy_plate_20200930_125752.csv')
FNAME_FILE = list(ROOT_DIR.rglob('*filtered/filenames_summary_tierpsy_plate_20200930_125752.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_results_files/filtered/filenames_summary_tierpsy_plate_20200930_125752.csv')
METADATA_FILE = list(ROOT_DIR.rglob('*wells_annotated_metadata.csv'))[0] #Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/wells_annotated_metadata.csv')

RAW_DATA_DIR = Path('/Volumes/Ashur Pro2/DiseaseScreen')

saveto = ROOT_DIR / 'ts_figures'
saveto.mkdir(exist_ok=True)

#vars and filters
CONTROL_STRAIN = 'N2'
n_subsample = 10 #10 wells for each strain
seed=20210205
feature = 'length'
frameRate = 25
minTrajDuration = 100 # mininum trajectory in seconds
smoothWindow = 1 # smooth feature over the window, in seconds 

strains_done = ['snf-11',
                'unc-49',
                'unc-43',
                'kcc-2',
                'snn-1',
                'tub-1',
                'glc-2',
                'avr-14',
                'unc-25',
                'unc-80',
                'nca-2',
                'bbs-1',
                'dys-1',
                'add-1',
                'cat-2']

#%%

if __name__ == "__main__":
    #set style for all figures
    plt.style.use(CUSTOM_STYLE)
    sns.set_style('ticks')
    
    meta = pd.read_csv(METADATA_FILE,
                       index_col=None)
    
    meta.dropna(axis=0,
                    subset=['worm_gene'],
                    inplace=True)
    meta.worm_gene.replace({'C43B7.2':'figo-1'}, inplace=True)
    
    # remove data from dates to exclude
    good_date = meta.query('@DATES_TO_DROP not in imaging_date_yyyymmdd').index
    # bad wells
    good_wells_from_gui = meta.query('is_bad_well == False').index
    meta = meta.loc[good_wells_from_gui & good_date,:]
    prestim_ix = [i for i,r in meta.iterrows() if 'prestim' in r.imgstore_name]
    meta = meta.loc[prestim_ix,:]
        
    genes = [g for g in meta.worm_gene.unique() if g != CONTROL_STRAIN]
    genes = [g.replace('C43B7.2', 'figo-1') for g in genes]
    genes = [g for g in genes if 'myo' not in g and 'unc-54' not in g]
    
    genes = list(set(genes) - set(strains_done))
    # strain_numbers = []
    
    #%%
    for count,g in enumerate(genes):
        print('Analysing {} {}/{}'.format(g, count+1, len(genes)))
        candidate_gene = g
        
        (saveto / candidate_gene).mkdir(exist_ok=True)
        
        meta_df, idx, gene_list = select_strains(candidate_gene,
                                                        CONTROL_STRAIN,
                                                        meta_df=meta)
        
        strain_lut, stim_lut = make_colormaps(gene_list,
                                            featlist=[],
                                            idx=idx,
                                            candidate_gene=candidate_gene,
                                            )
        
        subsample = meta_df.groupby('worm_gene').sample(n_subsample,
                                                        random_state=seed)
        
        
        for k,v in strain_lut.items():
            _meta = subsample.query('@k in `worm_gene`')
            
            for i,r in _meta.iterrows():
                _imgstore = RAW_DATA_DIR / 'Results' / r['imgstore_name'] / 'metadata_featuresN.hdf5'
                
                with pd.HDFStore(_imgstore) as fid:
                    ts = fid['/timeseries_data'].query("`well_name` == @r.well_name")
                    
                _worm_idx = ts.groupby('worm_index'
                                       ).filter(lambda x: x.shape[0] >= frameRate*minTrajDuration)
                
                _worm_grouped = _worm_idx.groupby('worm_index')
                
                for widx in _worm_grouped.groups.keys():
                    
                    y = _worm_grouped.get_group(widx)['length'].rolling(smoothWindow*frameRate).median()
                    
                    x = _worm_grouped.get_group(widx)['timestamp']/frameRate
                    x = x-x.min()
                    
                    #add in the peaks
                    peaks, _ = signal.find_peaks(-y.values, distance = 40*25)
                    
                    fig, ax = plt.subplots(figsize=[7,5])
                    plt.plot(x, y, color=v)
                    plt.plot(x.values[peaks], y.values[peaks], 'xk')
                    ax.set_xticks(list(range(0,int(round(x.max()))+1,30)))
                    plt.ylabel('length (microns)')
                    plt.xlabel('time (s)')
                    plt.tight_layout()
                    plt.savefig(saveto / candidate_gene / '{}_{}_length_ts.png'.format(k, widx))
                    plt.close('all')
                    