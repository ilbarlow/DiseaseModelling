#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:16:26 2020

@author: ibarlow

script to find which files are missing from my hard drive

"""

from pathlib import Path
import pandas as pd
from tierpsytools.hydra.hydra_filenames_helper import find_imgstore_videos
from tierpsytools.hydra.masked_videos_checks import *

N_VIDEOS = 180
id_string = 'disease_models*'

HD_DIR = Path('/Volumes/Ashur Pro2/DiseaseScreen/MaskedVideos')

BEHAV_DIR = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/MaskedVideos')

date_dirs = list(HD_DIR.glob('2020*'))

stim_types = ['bluelight', 'prestim', 'poststim']

for d in date_dirs:
    for s in stim_types:
        prestim = get_prestim_videos(d, bluelight=s)
        check_all_videos(prestim, d, bluelight=s)
        
        prestim.to_csv(d.parent.parent / 'AuxiliaryFiles' / d.stem / \
                       '{}_{}_masked_videos_to_check.csv'.format(d.stem, s),
                       index=False)


# missing_files = []
# for d in date_dirs:
#     if len(list(d.glob('disease_models*'))) != N_VIDEOS:
#         missing_files.append(d)

# files_to_find = [Path(str(d).replace('Ashur Pro2', 'behavgenom$/Ida/Data/Hydra'))
#                  for d in missing_files]

# missing = []
# for f in missing_files:
#     behav_dir = Path(str(f).replace('Ashur Pro2', 'behavgenom$/Ida/Data/Hydra'))
#                             )
#     behav_files = set(Path(str(f).replace('Ashur Pro2',
#                                            'behavgenom$/Ida/Data/Hydra')
#                             ).glob(id_string))
    
#     hd_files = set(f.glob(id_string))
#     hd_files = [str(h) for h in hd_files]
    
#     missing.append(set(behav_files).symmetric_difference(hd_files))

    
    