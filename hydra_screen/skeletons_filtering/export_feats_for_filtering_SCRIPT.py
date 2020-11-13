#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 11:57:18 2020

@author: ibarlow

Script for subsampling videos to get parameters for filtering the
disease models data

"""
import pandas as pd
from pathlib import Path
from random import sample

def _traj_length(grouped_series):
    traj_len =grouped_series['timestamp'].nunique()
    return traj_len

def _distance_traveled(grouped_series):
    dist = []
    for body_part in ['body', 'tail', 'midbody', 'head']:
        coords = ['_'.join([x, body_part])  for x in ['coord_x', 'coord_y']]
        dist.append(grouped_series[coords].apply(
            lambda x: x.diff().pow(2).sum(1).pow(0.5).sum())) # applies per group
    dist = pd.concat(dist, axis=1)
    return dist

def _mean_timeseries_values(grouped_series, names):
    mean_values = grouped_series[names].mean()
    return mean_values

#%% Input
root_dir = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/Results')

timeseries_names = ['width_midbody', 'length']
n_videos = 100

SAVETO = Path('/Users/ibarlow/OneDrive - Imperial College London/'+\
              'Documents/behavgenom_copy/DiseaseScreen/samples_feats_for_filtering.csv')


#%% Sample videos
files = [file for file in root_dir.rglob('*featuresN.hdf5')]
files = sample(files,
               n_videos)

#%% Get values
traj_lengths = []
dist_traveled = []
mean_values = []
for i,fname in enumerate(files):
    print('Reading file {}/{}...'.format(i+1, len(files)))
    print(fname)
    with pd.HDFStore(fname, 'r') as f:
        series = f['timeseries_data']
        coord_cols = [col for col in series.columns if 'coord' in col]
        id_cols = ['well_name',
                   'worm_index',
                   'timestamp']
        series = series[coord_cols + timeseries_names + id_cols]
        grouped_series = series.groupby(by='worm_index')
        traj_lengths.append(_traj_length(grouped_series))
        dist_traveled.append(_distance_traveled(grouped_series))
        mean_values.append(
            _mean_timeseries_values(grouped_series, timeseries_names))

fname_ids = pd.DataFrame({'file_id':list(range(len(files))),
                          'filename':files})
traj_lengths = pd.concat(
    [pd.DataFrame(
        {'file_id':[i]*x.shape[0],
         'worm_index':x.index,
         'traj_length':x.values}
        )
        for i,x in enumerate(traj_lengths)])

dist_traveled = pd.concat(
    [pd.DataFrame(
        {'file_id':[i]*x.shape[0],
         'worm_index':x.index,
         'max_dist_traveled':x.max(axis=1).values}
        )
        for i,x in enumerate(dist_traveled)])

mean_values =  pd.concat(
    [pd.DataFrame(
        {'file_id':[i]*x.shape[0],
         'worm_index':x.index,
         'width_midbody': x['width_midbody'],
         'length': x['length']}
        )
        for i,x in enumerate(mean_values)]).reset_index(drop=True)

all_data = pd.merge(fname_ids, traj_lengths, on='file_id', how='outer')
all_data = pd.merge(all_data, dist_traveled, on=['file_id','worm_index'], how='outer')
all_data = pd.merge(all_data, mean_values, on=['file_id','worm_index'], how='outer')

all_data.to_csv(SAVETO, index=False)
