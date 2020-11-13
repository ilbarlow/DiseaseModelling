#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 16:12:35 2020

@author: ibarlow

Script to add the missing filenames to the Annotations.hdf5 so that all
the well annotations can be completed

use tierpsy_dev environment

"""

import pandas as pd
from pathlib import Path
import h5py

PROJECT_DIR = Path('/Volumes/Ashur Pro2/DiseaseScreen')
N_PRESTIM_VIDEOS = 60

annotations_files = list(PROJECT_DIR.rglob('*annotations.hdf5'))

for c,f in enumerate(annotations_files):
    with pd.HDFStore(f) as fid:    
        annotated_files = fid['/filenames_df']
    
    if annotated_files.shape[0] == N_PRESTIM_VIDEOS:
        continue
        
    prestim_files = list(
        Path(str(f.parent).replace('AuxiliaryFiles',
                                   'MaskedVideos')).rglob('metadata.hdf5')
        )
    prestim_files = [Path(f) for f in prestim_files if 'prestim' in str(f)]
    prestim_files = ['{}/{}'.format(f.parent.name, f.name)
                     for f in prestim_files] 
        
    missing_files = list(set(prestim_files) - set(annotated_files['filename']))
    print(missing_files)
    
    last_fileid = annotated_files['file_id'].max()
    
    # add to annotation files
    for m in missing_files:
        last_fileid += 1
        _to_add = pd.Series({'file_id':last_fileid,
                             'filename':m})
        
        with pd.HDFStore(f,'r+') as fid:
            update = fid['/filenames_df'].copy()
            update = update.append(_to_add,
                                   ignore_index=True)
            update.to_hdf(
                fid,
                key='/filenames_df',
                index=False,
                mode='r+')              
            # print(update)
            wdir = fid.get_storer('/filenames_df').attrs.working_dir             
        
        with h5py.File(f, 'r+') as fid:
            fid["/filenames_df"].attrs["working_dir"] = str(wdir)
            
            
#%% 
# fix broken annotations files with file_id 59 at head with wrong well names
bad_file = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/AuxiliaryFiles/20200808/DiseaseScreen_20200904_164411_wells_annotations.hdf5')

with pd.HDFStore(bad_file,'r+') as fid:
    update = fid['/wells_annotations_df'].copy()
    update = update[update.file_id!=59].reset_index(drop=True)
    update.to_hdf(
        fid,
        key='/wells_annotations_df',
        index=False,
        mode='r+')   


bad_file2 = Path('/Users/ibarlow/Desktop/backup_update_annotations/AuxiliaryFiles/20200808/DiseaseScreen_20200904_164411_wells_annotations.hdf5')           


    # print(update)
#     wdir = fid.get_storer('/filenames_df').attrs.working_dir             

# with h5py.File(f, 'r+') as fid:
#     fid["/filenames_df"].attrs["working_dir"] = str(wdir)


# test_file = '/Volumes/Ashur Pro2/DiseaseScreen/AuxiliaryFiles/20200813/DiseaseScreen_20200908_134329_wells_annotations.hdf5'
# fid1 = pd.HDFStore(test_file, 'r+')
# well_annotations = fid1['/wells_annotations_df']



# # test_file2 = '/Users/ibarlow/Desktop/backup_update_annotations/AuxiliaryFiles/20200813/DiseaseScreen_20200908_134329_wells_annotations.hdf5'
# # fid2 = pd.HDFStore(test_file2, 'r+')
# #     # test_files = Path('/Users/ibarlow/Desktop/DiseaseScreen_20200904_164411_wells_annotations.hdf5')

    

    