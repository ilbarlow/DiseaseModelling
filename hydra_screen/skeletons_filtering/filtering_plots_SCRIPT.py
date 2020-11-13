#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:00:36 2020

@author: ibarlow

Plot the distributions of the samples Disease models data to decide if
filtering needs to be applied for extracting out the feature summaries


"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import statsmodels.api as sm
from pathlib import Path
import pandas as pd
from random import sample
import tables

ROOT_DIR = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen')
SAVETO = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/Figures')

all_data = pd.read_csv('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/samples_feats_for_filtering.csv')
files = list(all_data.filename.unique())

filtering_variables = {'traj_length': np.log10,
                       'max_dist_traveled': np.log10,
                       'width_midbody': np.log10,
                       'length': np.log10
                       }

filtering_variable_df = pd.DataFrame(data=filtering_variables.values(),
                                     index=filtering_variables.keys(),
                                     columns=['log_function'])

#%% functions

def summary_plots(x_array,
                  y_array,
                  log_base=np.log10,
                  variable_name=None,
                  distribution_type=stats.distributions.norm):
    """
    Function for making some histograms and cdfs from the sampled videos

    Parameters
    ----------
    x_array : TYPE, optional
        DESCRIPTION..
    y_array : TYPE, optional
        DESCRIPTION..
        
    log_base: TYPE
        DESCRIPTION. default in log10
    distribution_type : TYPE, optional
        DESCRIPTION. The default is stats.distributions.expon.
    variable_name : TYPE, optional
        DESCRIPTION. The default is 'traj_length'.

    Returns
    -------
    None.

    """

    x = np.sort(x_array) # assert to check the data; better to use argsort of x and y to make robust against if x and y are linked
    y = np.sort(y_array)
    if not (x == 0).any():
        logx = log_base(x)
    else:
        print('{} has zero values'.format(variable_name))
        xl = x.copy()
        xl[xl == 0] = xl[np.nonzero(x)][0]-10**-4 # subtract small value from all the data?
        logx = log_base(xl)
    
    log_type = str(log_base).split(' ')[-1]
    
    fig = plt.figure(figsize = [10,10])   
    # ECDF
    ax1 = fig.add_subplot(2,2,1) # set logscale on xaxes
    ax1 = plt.scatter(x,
                      y, 
                      s=0.5,
                      alpha=0.7,
                      c='k')
    ax1.axes.yaxis.set_label_text('ECDF')
    ax1.axes.xaxis.set_label_text(variable_name)
    
    # logged ECDF
    ax2 = fig.add_subplot(2,2,2)
    ax2 = plt.scatter(logx,
                      y,
                      s=0.5,
                      alpha=0.7,
                      c='k')
    ax2.axes.yaxis.set_label_text('ECDF')    
    ax2.axes.xaxis.set_label_text('{}_{}'.format(
        log_type,
        variable_name
        ))
    
    # histogram
    ax3 = fig.add_subplot(2,2,3)
    ax3 = sns.distplot(logx,
                       norm_hist=True,
                       color='k')
    ax3.axes.xaxis.set_label_text('{}_{}'.format(
        log_type,
        variable_name
        ))
    ax3.axes.yaxis.set_label_text('probability density')
    
    # qq plot
    ax4 = fig.add_subplot(2,2,4)
    sm.qqplot(logx,
              distribution_type,
              line='q',
              fit=True,
              ax=ax4,
              markersize=0.5,
              color='k')
    ax4.axes.xaxis.set_label_text('Theoretical quantiles ({})'.format(
        str(distribution_type).split('.')[-1]
        ))
    ax4.axes.yaxis.set_label_text('Sample quantiles ({}_{})'.format(
        log_type,
        variable_name
        ))

    return

# %% skeleton plot functions
def find_longest_traj(grouped_dataframe,
                      raw_dataframe,
                      nlimit=20,
                      mergeon=('filename','traj_length')):
    max_traj = pd.DataFrame(grouped_dataframe['traj_length'].nlargest(nlimit))
    max_traj.reset_index(inplace=True)
    max_worm_ids = pd.merge(max_traj,
                            raw_dataframe,
                            on=mergeon)[['filename',
                                         'worm_index',
                                         'traj_length']]
    return max_worm_ids


def plot_skeletons(skelInds, skeletons_array, color, axes=None):
    """
    
    Plot sekected skeletons from loaded data
    
    Parameters
    ----------
    skelInds : TYPE
        DESCRIPTION.
    skeletons_array : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if axes:
        axes.plot(skeletons_array[skelInds['skeleton_id'].values,:,0],
                 skeletons_array[skelInds['skeleton_id'].values,:,1],
                 alpha=0.3,
                 c=color)
    else:
        plt.plot(skeletons_array[skelInds['skeleton_id'].values,:,0],
                 skeletons_array[skelInds['skeleton_id'].values,:,1],
                 alpha=0.3,
                 c=color)
    return


def select_skeletons(wormIndex,
                     fname,
                     skeleton_df,
                     skeletons_array,
                     color,
                     subsample_rate = 10,
                     plot=True,
                     axes=None):
    """
    
    Select the skeletons with the corresponding worm id
    
    Parameters
    ----------
    wormIndex : TYPE
        DESCRIPTION.
    skeleton_df : TYPE
        DESCRIPTION.
    subsample_rate : TYPE, optional
        DESCRIPTION. The default is 10.

    Returns
    -------
    None.

    """
    skelInds = skeleton_df[
                    (skeleton_df['filename'] == fname) &
                    (skeleton_df['worm_index_joined'] == wormIndex)
                    ]
    skelInds = skelInds[skelInds['skeleton_id'] != -1]
    skelInds = skelInds[skelInds['was_skeletonized'] == 1]
    skelInds = skelInds.sort_values(by='frame_number')
    if skelInds.frame_number.diff().sum() > skelInds.shape[0]:
        skelInds_subsampled = np.split(skelInds,np.where(skelInds.frame_number.diff()>25)[0])
        for s in skelInds_subsampled:
            s = s[::subsample_rate]
            if plot:
                plot_skeletons(s, skeletons_array, color, axes=axes)
            else:
                return skelInds
    else:
        skelInds = skelInds[::subsample_rate]
        if plot:
            plot_skeletons(skelInds, skeletons_array, color, axes=axes)
        else:
            return skelInds
        
    return skelInds

#%%
raw_save_to = SAVETO / 'unfiltered'
raw_save_to.mkdir(exist_ok=True)
    
for i,r in filtering_variable_df.iterrows():
    x = np.sort(all_data[i].dropna())
    y = np.arange(1, x.size+1) / x.size #linearly spaced numbers between 0 and 1
    summary_plots(x,
                  y,
                  variable_name=i,
                  log_base=r.log_function)
    plt.savefig(raw_save_to / 'unfiltered_distplots_{}.png'.format(i))

#%%
# look at the plots to decide on filtering thresholds
threshold_dict = {'traj_length': [0, np.inf],
                  'width_midbody': [np.round(np.power(10, 1),2),
                                      np.round(np.power(10, 3),2)],
                  'max_dist_traveled': [np.round(np.power(10, 0),2),
                                        np.inf],
                  'length': [0,
                             np.round(np.power(10, 3.5), 2)]
                  }

for k,v in threshold_dict.items():
    filtering_variable_df.loc[k, 'min_thresh'] = v[0]
    filtering_variable_df.loc[k, 'max_thresh'] = v[1]
    filtering_variable_df.loc[k,
                              'min_max_ix'
                              ] = [set(all_data[(all_data[k] > threshold_dict[k][0]) &
             (all_data[k] < threshold_dict[k][1])].index.values)]

width_length_ix= set.intersection(filtering_variable_df.loc['width_midbody',
                                                        'min_max_ix'],
                                  filtering_variable_df.loc['length',
                                                        'min_max_ix'],
                                  filtering_variable_df.loc['traj_length',
                                                             'min_max_ix']
                                  )
print(len(width_length_ix)/all_data.shape[0])

width_length_data = all_data.loc[width_length_ix, :]

# make summary plots of all the variables with these filters applied
width_length_saveto = SAVETO / 'width_length_filters'
width_length_saveto.mkdir(exist_ok=True)
for i,r in filtering_variable_df.iterrows():
    x = np.sort(width_length_data[i].dropna())
    y = np.arange(1, x.size+1) / x.size #linearly spaced numbers between 0 and 1
    summary_plots(x,
                  y,
                  variable_name=i,
                  log_base=r.log_function)
    plt.savefig(width_length_saveto / 'filtered_distplots_{}.png'.format(i))



#find the longest trajectories
all_grouped = all_data.groupby('filename')
max_worm_ids_all = find_longest_traj(all_grouped,
                                     all_data,
                                     )

#get the worm indices of the timeseries that have been kept
width_length_data_grouped = width_length_data.groupby('filename')
max_worm_ids_width_length = find_longest_traj(width_length_data_grouped,
                                              width_length_data)

#make a look up table for all the kept worm ids
lut = pd.DataFrame((width_length_data_grouped['worm_index'].unique()))
lut.reset_index(inplace=True)

#%% import the skeletons
skel_unfiltered=[]
skeletons = {}
skel_cols= ['frame_number',
            'worm_index_joined',
            'skeleton_id',
            'coord_x',
            'coord_y',
            'was_skeletonized']
for i, f in enumerate(files[:10]):
    skelfile = f.replace('featuresN', 'skeletons')
    with tables.File(f, 'r') as fid:
        skeletons[skelfile] = fid.get_node('/coordinates/skeletons')[:]
    with pd.HDFStore(f, 'r') as fid:
        skel_unfiltered.append(fid['trajectories_data'][skel_cols])
    skel_unfiltered[i]['filename'] = f

skel_unfiltered = pd.concat(skel_unfiltered)    
skel_unfiltered.reset_index(drop=True, inplace=True)

#%% plots
font_params = {'fontsize':8}
for f in files[1:11]:
    skel_fname = f.replace('featuresN', 'skeletons')
    
    # unfiltered - sample 10% of the skeletons
    wormInds = list(
        skel_unfiltered[
            skel_unfiltered['filename'] == f
            ]['worm_index_joined'].unique()
                    )
    # wormInds = sample(wormInds, 
    #                   int(0.10 * len(wormInds)))
    
    #filtered - sample 5% of the skeletons
    filtered_wormInds = list(
        lut[lut['filename'] == f]['worm_index'].values[0]
                            )
    # filtered_wormInds = sample(filtered_wormInds,
    #                            int(0.1 * len(filtered_wormInds))
    #                            )   
    
    fig, ax = plt.subplots(2, 2, sharey=True, sharex=True, figsize = [20,20])

    seq_cmap = plt.get_cmap('plasma',
                            len(wormInds))
    for c, i in enumerate(wormInds):
        skelInds = select_skeletons(i,
                                    f,
                                    skel_unfiltered,
                                    skeletons[skel_fname],
                                    seq_cmap(c),
                                    subsample_rate = 10,
                                    plot=False)
 
        plot_skeletons(skelInds,
                        skeletons[skel_fname],
                        seq_cmap(c),
                        ax[0,0])
    
    ax[0,0].set_title('10% unfiltered worm indices; {} worm ids'.format(len(wormInds)),
              fontdict=font_params)
    ax[0,0].set_ylim([0, 35000])
    ax[0,0].set_xlim([0,35000])

    longest = max_worm_ids_all[max_worm_ids_all.filename == f].sort_values(by='traj_length',
                                                                    ascending=False)#['worm_index']).values
    skel_lengths = []
    for c, l in longest.iterrows():
        print (l)
        skelInds = select_skeletons(l.worm_index,
                                    f,
                                     skel_unfiltered,
                                     skeletons[skel_fname],
                                     color='r',
                                     subsample_rate=10,
                                     axes=ax[1,0])
        if skelInds.shape[0] > 0:
            skel_length = (l.traj_length, skelInds.shape[0])
            break

    ax[1,0].set_xlabel('unfiltered longest skel: {} frames'.format(skel_length),
              fontdict=font_params)  
    del skel_length


# filtered
    seq_cmap = plt.get_cmap('plasma',
                            len(filtered_wormInds))

    
    for c,i in enumerate(filtered_wormInds):
        skelInds = select_skeletons(i,
                                    f,
                         skel_unfiltered,
                         skeletons[skel_fname],
                         seq_cmap(c),
                         subsample_rate = 10,
                         plot=False)
        plot_skeletons(skelInds,
                       skeletons[skel_fname],
                       seq_cmap(c),
                       ax[0,1])   

    ax[0,1].set_title('10% filtered worm indices; {} worm ids'.format(len(filtered_wormInds)),
              fontdict=font_params)
    
    longest = max_worm_ids_width_length[max_worm_ids_width_length.filename == f].sort_values(by='traj_length',
                                                                    ascending=False)#['worm_index']).values
    skel_lengths = []
    for c, l in longest.iterrows():
        skelInds = select_skeletons(l.worm_index,
                                    f,
                                     skel_unfiltered,
                                     skeletons[skel_fname],
                                     color='r',
                                     subsample_rate=10,
                                     axes=ax[1,1])
        if skelInds.shape[0] > 0:
            print(l)
            skel_length = (l.traj_length, skelInds.shape[0])
            break
    # plt.ylim([600, 1600])
    # plt.axis('equal')
    ax[1,1].set_xlabel('filtered longest skel: {} frames'.format(skel_length),
              fontdict=font_params)
    del skel_length
    
    # fig.tight_layout()
    
    plt.savefig(SAVETO / 'skeletons_coords_{}.png'.format(f.split('/')[-2].replace('.','_')),
                dpi=400)
    
    plt.close('all')





