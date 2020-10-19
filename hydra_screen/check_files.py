#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:16:26 2020

@author: ibarlow

script to find which files are missing from my hard drive

"""

from pathlib import Path
import pandas as pd
from tierpsytools.hydra.hydra_filenames_helper import find_imgstore_videos, raw_to_masked, raw_to_featuresN
import datetime

HD_DIR = Path('/Volumes/AshurPro2/DiseaseScreen')
BEHAV_DIR = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen')

#%%
def check_masked_results(df):
    df['masked_vid'] = df.full_path.apply(raw_to_masked)
    df['missing_masked'] = df.masked_vid.apply(lambda x: not x.exists())
    df['results'] = df.full_path.apply(raw_to_featuresN)
    df['missing_results'] = df.results.apply(lambda x: not x.exists())
    return df

#%%
# df = find_imgstore_videos(BEHAV_DIR)
# behav_genom_df = check_masked_results(df)
# behav_genom_df.to_csv(BEHAV_DIR / 'AuxiliaryFiles' / \
#           '{}_checked_files.csv'.format(
#               datetime.datetime.today().strftime('%Y%m%d'))
#           )
# missing_files = behav_genom_df[behav_genom_df.all(axis=1)]
# missing_files.to_csv(BEHAV_DIR / 'AuxiliaryFiles' /\
#                      '{}_missing_files.csv'.format(
#                          datetime.datetime.today().strftime('%Y%m%d')))

#%% load behavgenom results
behavgenom_df = pd.read_csv('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/AuxiliaryFiles/20200914_checked_files.csv',
                            index_col=False)
df = behavgenom_df.drop(columns = ['Unnamed: 0',
                                   'masked_vid',
                                   'missing_masked',
                                   'results',
                                   'missing_results'])

#%% and do checks against external hard drive
SAVE_HD = HD_DIR / 'AuxiliaryFiles'
SAVE_HD.mkdir(exist_ok=True)

hd_df = df.copy()
hd_df['full_path'] = hd_df.full_path.apply(lambda x: str(x).replace(str(BEHAV_DIR),
                                                                    str(HD_DIR)))
hd_df = check_masked_results(hd_df)
hd_df.to_csv(SAVE_HD / '{}_checked_files.csv'.format(
                 datetime.datetime.today().strftime('%Y%m%d')),
             index=False
             )

missing_files_hd = hd_df[hd_df.all(axis=1)]
missing_files_hd.to_csv(HD_DIR / 'AuxiliaryFiles' /\
                     '{}_missing_files.csv'.format(
                         datetime.datetime.today().strftime('%Y%m%d')),
                     index=False
                     )
