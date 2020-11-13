#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 15:13:58 2020

@author: ibarlow

Script to add the wells annotations to the metadata.

Wells were annotated using Luigi's WellAnnotator software and results
are saved as *wells_annotations.hdf5 with two tables:
    /filenames_df
    /wells_annotations_df

annotations are coded so that 1 = good well, >1 are bad wells for range of
reasons (precipitation, wet), and 0 = not annotated

This script will:
    - extract the annotations from the hdf5 file
    - propagate 'prestim' annotations to bluelight and poststim
    - add the wells_annotations onto the metadata.csv and then
add another column is_bad as a boolean for filtering during analysis

"""


import pandas as pd
from pathlib import Path
from tierpsytools.hydra.match_bluelight_videos import match_bluelight_videos_in_folder


PROJECT_DIR = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen')
BEHAVGENOM_DIR = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen')

search_string = '*wells_annotations.hdf5'
table_names = ['/filenames_df',
                '/wells_annotations_df']
merge_on = ['imgstore_name',
            'well_name']

#find all masked video files:
matched_rawvids = match_bluelight_videos_in_folder(BEHAVGENOM_DIR / 'RawVideos')

bluelight_names = ['imgstore_prestim',
                   'imgstore_bluelight',
                   'imgstore_poststim']

#%%
annotations_files = list(PROJECT_DIR.rglob(search_string))
annotations = []
for f in annotations_files:
    with pd.HDFStore(f) as fid:
        _fnames = fid[table_names[0]]
        _markings = fid[table_names[1]]

    # find matching bluelight and poststim files
    # _fnames = expand_filename(_fnames)
    annotations.append(_fnames.merge(_markings[['file_id',
                                                'well_label',
                                                'well_name']],
                                    on='file_id',
                                    right_index=True,
                                    validate='one_to_many')
                                      )
annotations = pd.concat(annotations)
annotations.reset_index(drop=True,
                        inplace=True)
annotations['imgstore'] = annotations.filename.apply(lambda x: x.split('/')[0])

matched_vids = matched_rawvids.merge(annotations,
                                 how='outer',
                                 left_on=['imgstore_prestim'],
                                 right_on=['imgstore'],
                                 validate='one_to_many')
matched_vids.drop(columns='imgstore',
                  inplace=True)
matched_long = matched_vids.melt(id_vars = ['well_name', 'well_label'],
                                value_vars = bluelight_names,
                                value_name = 'imgstore')
matched_long.drop(columns = ['variable'],
                  inplace=True)
matched_long['is_bad_well'] = matched_long['well_label'] != 1

#%%

metadata_fname = list(PROJECT_DIR.rglob('metadata.csv'))[0]
metadata = pd.read_csv(metadata_fname)
if metadata['imgstore_name'].isna().sum() > 0:
    print('Nan values in imgstore names')
    metadata = metadata[metadata['imgstore_name'].notna()]
metadata.loc[:,'imgstore'] = metadata['imgstore_name'].apply(
    lambda x: x.split('/')[1])

update_metadata = metadata.merge(matched_long,
                                 on=['imgstore',
                                     'well_name'],
                                 how='outer')
update_metadata.drop(columns='imgstore',
                    inplace=True)
# drop unannotated wells
update_metadata = update_metadata[update_metadata.well_label.notna()
                                  ].reset_index(drop=True)

update_metadata.to_csv(PROJECT_DIR / 'AuxiliaryFiles' / 'wells_annotated_metadata.csv',
                       index=False)
