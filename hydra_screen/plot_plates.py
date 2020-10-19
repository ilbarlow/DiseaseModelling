#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 10:30:13 2020

@author: ibarlow


Script for plotting the plate trajectories for all the DiseaseScreen experiments

Uses Saul's function

Run in tierpsytools environment
"""

from pathlib import Path
import sys
import matplotlib.pyplot as plt

sys.path.insert(0,'/Users/ibarlow/tierpsy-tracker')#'/tierpsy') #/analysis/split_fov

# import tierpsy
from tierpsy.analysis.split_fov.helper import CAM2CH_df, serial2channel, parse_camera_serial
from tierpsy.analysis.split_fov.FOVMultiWellsSplitter import FOVMultiWellsSplitter

from tierpsytools.plot.plot_plate_trajectories_with_raw_video_background import plot_plate_trajectories

PROJECT_DIR = Path('/Volumes/AshurPro2/DiseaseScreen')

saveDir = Path('/Users/ibarlow/Desktop/test_plates')
featfiles = list(PROJECT_DIR.rglob('*featuresN.hdf5'))

#%%
filestem_list = []
featurefile_list = []
for c,f in enumerate(featfiles):
    # obtain file stem
    filestem = f.parent.parent / f.parent.stem

    # only record featuresN filepaths with a different file stem as we only
    # need 1 camera's video per plate to find the others
    if filestem not in filestem_list:
        filestem_list.append(filestem)
        featurefile_list.append(f)
#%%
#loop through to make the plate images
bad_files = []
for c, f in enumerate(featurefile_list):
    saveto = saveDir / f.parent.parent.stem
    saveto.mkdir(exist_ok=True)

    saveto = saveto / str(f).split('_')[3]
    saveto.mkdir(exist_ok=True)

    try:
        if Path(str(saveto / f.parent.stem) + '.png').exists():
            print (c)
            # continue
        else:
            print('makiing {}'.format(f))
            plot_plate_trajectories(f, saveto)
            plt.close('all')

    except Exception as error:
        print(error)
        bad_files.append((f, error))

with open(saveDir / 'bad_files.csv', 'w+') as fid:
    for b in bad_files:
        fid.write('{}, {}'.format(str(b[0]), b[1]) + '\n')
