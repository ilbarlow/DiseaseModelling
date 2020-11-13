#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 08:53:36 2020

@author: ibarlow

Script to do initial first glance of the Disesase screen data and output the 
results files for the candidate mutants that may be interesting for Madeleine

"""

import pandas as pd
from pathlib import Path


FEAT_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/features_summary_tierpsy_plate_20200912_224547.csv')
FNAME_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/filenames_summary_tierpsy_plate_20200912_224547.csv')

METADATA_FILE = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/metadata.csv')


SAVETO = FEAT_FILE.parent.parent / 'Figures'
SAVETO.mkdir(exist_ok=True)

# %%
if __name__ == '__main__':

    feat = pd.read_csv(FEAT_FILE,
                       comment='#')
    fname = pd.read_csv(FNAME_FILE,
                        comment='#')

    meta = pd.read_csv(METADATA_FILE, index_col=None)
    # meta = meta[meta.worm_strain.notna()]
    
    assert meta.worm_strain.unique().shape[0] == meta.worm_gene.unique().shape[0]
    
    meta = meta.loc[meta['imaging_date_yyyymmdd']!=float(DATES_TO_DROP), :]
    
    #%% find files with gpb-2, unc-25, bbs-1, bbs-2, N2
    selected_genes = ['gpb-2', 'unc-25', 'bbs-1', 'bbs-2', 'N2']
    
    select_rows = meta[meta.worm_gene.isin(selected_genes)] 
    select_rows = select_rows[select_rows.imgstore_name.notna()]
   
    #prestim  only
    select_rows = select_rows[select_rows.imgstore_name.str.contains('prestim')]
    
    #select only the ones 
    select_rows_grouped = select_rows.groupby('imgstore_name')
    
    to_export = select_rows_grouped.apply(lambda x: x.drop_duplicates(
                subset='worm_strain'))

    N2_to_keep = to_export[to_export.worm_strain != 'N2']['imaging_date_yyyymmdd'].unique()

    to_export = to_export[to_export.imaging_date_yyyymmdd.isin(N2_to_keep)]
    to_export.reset_index(drop=True,
                          inplace=True)
    
    to_export.to_csv(SAVETO / 'files_for_madeleine.csv',
                     index=False)
    
    to_export[to_export.worm_gene == 'gpb-2'].to_csv(SAVETO.parent / 'gpb2_files.csv',
                                                     index=False)
    