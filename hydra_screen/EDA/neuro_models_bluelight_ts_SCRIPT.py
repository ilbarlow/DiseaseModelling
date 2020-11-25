#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 17:11:28 2020

@author: ibarlow

Get timeseries data from the disease screen for all the neurodisease mutants only
- most of this script is copied from Luigi's original bluelight analysis

Useful mostly for eyeballing the phenotypes


"""


import time
import numpy as np
import pandas as pd
import seaborn as sns

from tqdm import tqdm
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from luigi_helper import (
    # find_motion_changes,
    # read_metadata,
    plot_stimuli,
    load_bluelight_timeseries_from_results,
    # just_load_one_timeseries,
    count_motion_modes,
    get_frac_motion_modes,
    HIRES_COLS
    )

HD = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen')
METADATA_FNAME = HD / 'AuxiliaryFiles/wells_annotated_metadata.csv'

is_reload_timeseries_from_results = True
is_recalculate_frac_motion_mode = True
timeseries_fname = HD / 'timeseries.hdf5'


SAVETO = HD / 'Figures/bluelight_timeseries'
SAVETO.mkdir(exist_ok=True)

ROOT_DIR = Path('/Volumes/Ashur Pro2/DiseaseScreen')

#analysis variables
CONTROL_STRAIN = 'N2'
DATES_TO_DROP = ['20200626',
                '20200730',
                '20200801',
                '20200806',
                '20200808',
                '20200811',
                '20200813',
                # '20200902',
                # '20200903',
                # '20200904',
                ]

#%%
def align_bluelight_meta(metadata_df):
    # reshape the dataframe slightly so we only have one line per (plate,well)
    pre_df = metadata_df[metadata_df['imgstore_name'].str.contains('prestim')]
    blue_df = metadata_df[metadata_df['imgstore_name'].str.contains('bluelight')]
    post_df = metadata_df[metadata_df['imgstore_name'].str.contains('poststim')]

    cols_to_join = list(set(pre_df.columns) - set(['imgstore_name']))

    md_df = pd.merge(pre_df, blue_df,
                     how='outer',
                     on=cols_to_join,
                     suffixes=('_pre', '_blue'))
    md_df = pd.merge(md_df, post_df,
                     how='outer',
                     on=cols_to_join,
                     suffixes=('', '_post'))
    md_df.rename(columns={'imgstore_name': 'imgstore_name_post'},
                 inplace=True)
    
    return md_df


# %%uses Luigi's own bootstrapping with ci functions
def my_sum_bootstrap(data, alpha=95):
    if isinstance(data, pd.core.series.Series):
        data_clean = data.values
        data_clean = data_clean[~np.isnan(data_clean)]
    else:
        data_clean = data[~np.isnan(data)]
    stafunction = np.sum
    n_samples = 1000
    funevals = np.ones(n_samples) * np.nan
    maxint = len(data_clean)
    sampling = np.random.randint(0, maxint, size=(n_samples, maxint))
    for evalc, sampling_ind in enumerate(sampling):
        funevals[evalc] = stafunction(data_clean[sampling_ind])
    pctiles = np.percentile(funevals, (50 - alpha/2, 50 + alpha/2))
    return tuple(pctiles)

def get_frac_motion_modes_with_ci(df, is_for_seaborn=False):
    """get_frac_motion_modes_with_ci
    divide number of worms in a motion mode by
    the total number of worms in that frame.
    Does *not* require that data have already been consolidated
    by a .groupby('timestamp')
    """
    # sum all n_worms across wells
    total_n_worms_in_frame = df.groupby(
        ['worm_gene', 'timestamp'], observed=True)['n_worms'].transform(sum)
    tmp = df.drop(columns='n_worms')
    # transform n_worms in frac_worms
    tmp = tmp.divide(total_n_worms_in_frame, axis=0)
    tmp.rename(lambda x: x.replace('n_worms_', 'frac_worms_'),
               axis='columns',
               inplace=True)
    if is_for_seaborn:
        # will use seaborn to get confidence intervals (estimator=sum),
        # don't need to calculate them here
        return tmp

    # implied else
    # now use bootstrap to get CIs
    out = tmp.groupby(['worm_gene', 'timestamp'],
                      observed=True).agg([np.sum, my_sum_bootstrap])

    def col_renamer(cname):
        cname_out = cname[0]+'_ci' if ('ci' in cname[1] or
                                       'bootstrap' in cname[1]) else cname[0]
        return cname_out
    out.columns = [col_renamer(c) for c in out.columns]

    return out

#%%
if __name__ == '__main__':
    
    meta = pd.read_csv(METADATA_FNAME, index_col=None)  
    # meta = meta[meta.worm_strain.notna()]
    
    assert meta.worm_strain.unique().shape[0] == meta.worm_gene.unique().shape[0]
    #check files in summary vs files in metadata
    
    meta.loc[:,'date_yyyymmdd'] = meta['date_yyyymmdd'].apply(lambda x: str(int(x)))
    
    #drop nan wells
    meta.dropna(axis=0,
                subset=['worm_gene'],
                inplace=True)

    # remove data from dates to exclude
    good_date = meta.query('@DATES_TO_DROP not in imaging_date_yyyymmdd').index

    # bad wells
    good_wells_from_gui = meta.query('is_bad_well == False').index

    meta = meta.loc[good_wells_from_gui & good_date,:]
    # strain sets
    genes = [g for g in meta.worm_gene.unique()]
    
    # Neuro disease strains
    # neuro_genes = list(set(genes) - set(myosin_genes))
    neuro_genes = [g for g in genes if g != CONTROL_STRAIN]
    neuro_genes.sort()
    
    neuro_locs = list(meta.query('@neuro_genes in worm_gene').index)
    N2_locs = meta.query('@CONTROL_STRAIN in worm_gene').sample(300).index
    
    neuro_locs.extend(N2_locs)
    neuro_genes.append(CONTROL_STRAIN)
    
    cmap = list(np.flip(sns.color_palette('cubehelix', len(neuro_genes)-1)))
    cmap.append((0.6, 0.6, 0.6))
    
    neuro_dict =  {r.worm_strain : r.worm_gene for
                   i,r in meta[['worm_strain',
                                'worm_gene']].drop_duplicates().iterrows()
                   }
    neuro_dict_inv = {v:k for k,v in neuro_dict.items()}
    
     #%% Only do analysis on the disease strains    
    meta = meta.loc[neuro_locs,:]
    
    meta = align_bluelight_meta(meta)
    
    neuro_genes = list(set(meta.worm_gene.unique()).intersection(neuro_genes))
    # %% read timeseries from results or disk

    if is_reload_timeseries_from_results:
        # this uses tierpytools under the hood
        timeseries_df, hires_df = load_bluelight_timeseries_from_results(
            meta,
            ROOT_DIR / 'Results')
        # save to disk
        timeseries_df.to_hdf(timeseries_fname, 'timeseries_df', format='table')
        hires_df.to_hdf(timeseries_fname, 'hires_df', format='table')
    else:  # from disk, then add columns
        # dataframe from the saved file
        timeseries_df = pd.read_hdf(timeseries_fname, 'timeseries_df')
        hires_df = pd.read_hdf(timeseries_fname, 'hires_df')
      
    #%% add in information about the replicates
    date_to_repl = pd.DataFrame({'date_yyyymmdd': [
                                                   #  '20200730',
                                                   # '20200801',
                                                   # '20200806',
                                                   # '20200808',
                                                   # '20200811',
                                                   # '20200813',
                                                   '20200902',
                                                   '20200903',
                                                   '20200904'
                                                   ],
                             'replicate': [1,2,3]}) #[1, 1, 2, 2, 3, 3]
    timeseries_df = pd.merge(timeseries_df, date_to_repl,
                             how='left',
                             on='date_yyyymmdd')
    hires_df = pd.merge(hires_df, date_to_repl,
                        how='left',
                        on='date_yyyymmdd')

    #map worm gene on worm strain
    # timeseries_df['worm_gene'] = timeseries_df['worm_strain'].map(neuro_dict)
    # hires_df['worm_gene'] = hires_df['worm_strain'].map(neuro_dict)
    
    timeseries_df['worm_strain'] = timeseries_df['worm_gene'].map(neuro_dict_inv)
    hires_df['worm_strain'] = hires_df['worm_gene'].map(neuro_dict_inv)

    #%% make d/v signed features absolute as in hydra d/v is not assigned 
    
    # other feats to abs:
    # a few features only make sense if we know ventral/dorsal
    feats_to_abs = ['speed',
                    'relative_to_body_speed_midbody',
                    'd_speed',
                    'relative_to_neck_angular_velocity_head_tip']
    feats_to_abs.extend([f for f in timeseries_df.columns
                         if f.startswith(('path_curvature',
                                          'angular_velocity'))])
    for feat in feats_to_abs:
        timeseries_df['abs_' + feat] = timeseries_df[feat].abs()
        
   
    # %%% hand-picked features from the downsampled dataframe
    
    plt.close('all')
    
    feats_toplot = ['speed',
                    'abs_speed',
                    'angular_velocity',
                    'abs_angular_velocity',
                    'relative_to_body_speed_midbody',
                    'abs_relative_to_body_speed_midbody',
                    'abs_relative_to_neck_angular_velocity_head_tip',
                    'speed_tail_base',
                    'length',
                    'major_axis',
                    'd_speed',
                    'head_tail_distance',
                    'abs_angular_velocity_neck',
                    'abs_angular_velocity_head_base',
                    'abs_angular_velocity_hips',
                    'abs_angular_velocity_tail_base',
                    'abs_angular_velocity_midbody',
                    'abs_angular_velocity_head_tip',
                    'abs_angular_velocity_tail_tip']
    from pandas.api.types import CategoricalDtype

    for ng in neuro_genes:
        if ng != CONTROL_STRAIN:
            to_plot = [CONTROL_STRAIN, ng]
            _plot_df = timeseries_df.query('motion_mode != 0 and @to_plot in worm_gene')
            _plot_df['worm_gene'] = _plot_df['worm_gene'].astype(
                CategoricalDtype(categories=to_plot, ordered=False)
                )
            with PdfPages(SAVETO / '{}_downsampled_feats.pdf'.format(ng), keep_empty=False) as pdf:
                for feat in tqdm(feats_toplot):
                    fig, ax = plt.subplots()
                    sns.lineplot(x='time_binned_s',
                                 y=feat,
                                 hue='motion_mode',
                                 style='worm_gene',
                                 size_order=to_plot,
                                 data=_plot_df,
                                 estimator=np.mean,
                                 ci='sd',
                                 legend='full')
                    plot_stimuli(ax=ax, units='s')
                    pdf.savefig(fig)
                    plt.close(fig)

# high resolution plots
    for ng in neuro_genes:
        if ng != CONTROL_STRAIN:
            to_plot = [CONTROL_STRAIN, ng]
            _plot_df = timeseries_df.query('motion_mode != 0 and @to_plot in worm_gene')
            _plot_df['worm_gene'] = _plot_df['worm_gene'].astype(
                CategoricalDtype(categories=to_plot, ordered=False)
                )                    
            
        with PdfPages(SAVETO / '{}_hires_feats_noerrs.pdf'.format(ng), keep_empty=False) as pdf:
            for col in tqdm(HIRES_COLS):
                if hires_df[col].dtype == 'float32':
                    fig, ax = plt.subplots()
                    sns.lineplot(x='timestamp', y=col,
                                 hue='motion_mode',
                                 style='worm_gene',
                                 data=_plot_df,
                                 estimator='median', ci=None,
                                 legend='full', ax=ax)
                    plot_stimuli(ax=ax, units='frames')
                    pdf.savefig(fig)
                    plt.close(fig)
                    
    # v noisy hard to see anything
    
    # %%% fraction motion modes

    # get motion_mode stats
    motion_mode_by_well = count_motion_modes(hires_df)
    
    # plots with no error bars:
    
    # aggregate data from all different wells, but keep strains separate
    motion_mode = motion_mode_by_well.groupby(['worm_gene', 'timestamp'],
                                              observed=True).sum()
    # compute fraction of worms in each motion mode (works with 'worm_strain' too)
    frac_motion_mode = get_frac_motion_modes(motion_mode)

    # %% motion modes no error bars
    
    plt.close('all')
    for strain_name, frmotmode_strain in frac_motion_mode.groupby(
                'worm_gene'):
    
        with PdfPages(SAVETO/ '{}_motion_mode_frac_no_errs.pdf'.format(strain_name),
                      keep_empty=False) as pdf:
        

            # look at area between curves here
            frmotmode_strain[['frac_worms_fw',
                              'frac_worms_st',
                              'frac_worms_bw',
                              'frac_worms_nan']].cumsum(axis=1).droplevel(
                                  'worm_gene').plot()
            plt.gca().set_ylim([0, 1])
            plt.gca().set_ylabel('cumulative fraction')
            plt.gca().set_title(strain_name)
            plot_stimuli(units='frames', ax=plt.gca())
            pdf.savefig()
            plt.close()
    
            # just timeseries here instead
            frmotmode_strain.droplevel(
                'worm_gene').plot(
                    y=['frac_worms_fw', 'frac_worms_st', 'frac_worms_bw'])
            plt.gca().set_ylim([0, 1])
            plt.gca().set_ylabel('fraction')
            plot_stimuli(units='frames', ax=plt.gca())
            plt.gca().set_title(strain_name)
            pdf.savefig()
            plt.close()
            
    # %% make some nicer plotes with confindence intervals using luigi's bootstrapping
    tic = time.time()

    if is_recalculate_frac_motion_mode:
        frac_motion_mode_with_ci = get_frac_motion_modes_with_ci(
            motion_mode_by_well)
        for col in ['frac_worms_bw_ci', 'frac_worms_st_ci',
                    'frac_worms_fw_ci', 'frac_worms_nan_ci']:
            
            frac_motion_mode_with_ci[col+'_lower'] = \
                frac_motion_mode_with_ci[col].apply(lambda x: x[0])
            
            frac_motion_mode_with_ci[col+'_upper'] = \
                frac_motion_mode_with_ci[col].apply(lambda x: x[1])
            frac_motion_mode_with_ci.drop(columns=col, inplace=True)
    
        frac_motion_mode_with_ci.to_hdf(timeseries_fname,
                                        'frac_motion_mode_with_ci',
                                        format='table')
    else:
        frac_motion_mode_with_ci = pd.read_hdf(timeseries_fname,
                                               'frac_motion_mode_with_ci')
    
    fps = 25
    frac_motion_mode_with_ci = frac_motion_mode_with_ci.reset_index()
    frac_motion_mode_with_ci['time_s'] = (frac_motion_mode_with_ci['timestamp']
                                          / fps)
    print('Time elapsed: {}s'.format(time.time()-tic))

    # add the strain to motion mode df
    
    # TODO check worm strain
    frac_motion_mode_with_ci['worm_strain'] = frac_motion_mode_with_ci['worm_gene'].map(neuro_dict_inv)
    
    # tic = time.time()
    
    # foo = get_frac_motion_modes_with_ci(motion_mode_by_well, is_for_seaborn=True)
    
    # fig, ax = plt.subplots()
    # sns.lineplot(x='timestamp', y='frac_worms_fw', style='worm_gene',
    #              estimator='sum', ci=95,
    #              data=foo.reset_index(), ax=ax)
    # # frac_motion_mode.plot(y='frac_worms_fw', linestyle='--', ax=ax)  # check
    # sns.lineplot(x='timestamp', y='frac_worms_bw', style='worm_strain',
    #              estimator='sum', ci=95,
    #              data=foo.reset_index(), ax=ax)
    # # frac_motion_mode.plot(y='frac_worms_bw', linestyle='--', ax=ax)  # check
    # plot_stimuli(ax=ax, units='frames')
    # fig.savefig(figures_dir / 'frac_modes_sns.pdf')
    
    # print('Time elapsed: {}s'.format(time.time()-tic))    
#%%
    def plot_frac(df, modecolnames, ax=None, **kwargs): #style_dict={'N2': '--','CB4856': '-'}
        """plot_frac
        plots modecolname of df with shaded errorbar
        example:
            plot_frac(frac_motion_mode_with_ci, 'frac_worms_fw', ax=ax)
        """
        if ax is None:
            ax = plt.gca()
        coldict = {'frac_worms_fw': 'tab:green',
                   'frac_worms_st': 'tab:purple',
                   'frac_worms_bw': 'tab:orange',
                   'frac_worms_nan': 'tab:gray'}
        namesdict = {'frac_worms_fw': 'forwards',
                     'frac_worms_st': 'stationary',
                     'frac_worms_bw': 'backwards',
                     'frac_worms_nan': 'not defined'}
        # styledict = {'N2': '--',
        #              'CB4856': '-'}
        if len(ax) != 1:
            assert(len(ax) == df['worm_gene'].nunique())
    
        for ii, (strain, df_g) in enumerate(df.groupby('worm_gene')):
            # df_g = df_g.droplevel('worm_strain')
            if len(ax) != 1:
                this_ax = ax[ii]
            else:
                this_ax = ax
            for col in modecolnames:
                df_g.plot(x='time_s',
                          y=col, ax=this_ax,
                          color=coldict[col],
                          label='_nolegend_',
                          **kwargs)
                          # linestyle=styledict[strain],
                          # label=strain+' '+col.split('_')[-1],
    
                lower = df_g[col+'_ci_lower']
                upper = df_g[col+'_ci_upper']
                this_ax.fill_between(x=df_g['time_s'],
                                     y1=lower.values,
                                     y2=upper.values,
                                     alpha=0.3,
                                     facecolor=coldict[col])
            # plt.plot([-1, -1], [-1, -1], color='black',
            #          linestyle=styledict[strain],
            #          label=strain)
            this_ax.set_ylabel('fraction of worms')
            this_ax.set_xlabel('time, (s)')
            this_ax.set_title(strain)
            this_ax.set_ylim((0, 1))
    
        # plt.legend(frameon=False, loc='upper left')
            this_ax.get_legend().remove()
        for i, col in enumerate(modecolnames):
            xm, xM = this_ax.get_xlim()
            x = xm + 0.99 * (xM - xm)
            y = 0.95 - i * 0.05
            this_ax.text(x, y, namesdict[col], color=coldict[col],
                         fontweight='heavy',
                         horizontalalignment='right')
        return

# %%
    for ng in neuro_genes:
        if ng != CONTROL_STRAIN:
            to_plot = [CONTROL_STRAIN, ng]
            _plot_df = frac_motion_mode_with_ci.query('@to_plot in worm_gene')
            _plot_df['worm_gene'] = _plot_df['worm_gene'].astype(
                CategoricalDtype(categories=to_plot, ordered=False)
                )
            fig, axs = plt.subplots(1, 2, figsize=(12.8, 4.8), sharey=True)
            axs = [axs[-1], axs[0]]
            plot_frac(_plot_df.reset_index(),
                      ['frac_worms_fw', 'frac_worms_st', 'frac_worms_bw'],
                      ax=axs)
        for ax in axs:
            plot_stimuli(ax=ax, units='s')
        plt.tight_layout()
        fig.subplots_adjust()
        fig.savefig(SAVETO / '{}_frac_worms_motion_ci.pdf'.format(ng))
        plt.close('all')
 
