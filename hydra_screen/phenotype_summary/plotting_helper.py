#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 09:44:49 2020

@author: ibarlow
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import numpy as np
import sys
sys.path.insert(0, '/Users/ibarlow/Documents/GitHub/pythonScripts/DiseaseModelling/hydra_screen/phenotype_summary')

from helper import STIMULI_ORDER, BLUELIGHT_WINDOW_DICT

def plot_colormaps(strain_lut, stim_lut, saveto):
    """
    

    Parameters
    ----------
    strain_lut : TYPE
        DESCRIPTION.
    stim_lut : TYPE
        DESCRIPTION.
    saveto : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    sns.set_style('dark')
    fig, ax = plt.subplots(1,2, figsize = [10,2])
    # ax=ax.flatten()
    for c, (axis, lut) in enumerate([(ax[0], strain_lut), (ax[1], stim_lut)]):
        axis.imshow([[v for v in lut.values()]])
        axis.axes.set_xticks(range(0, len(lut), 1))
        axis.axes.set_xticklabels(lut.keys(), rotation=90, fontsize=12)
        axis.axes.set_yticklabels([])
    fig.tight_layout()
    if saveto==None:
        return
    else:
        fig.savefig(saveto / 'colormaps.png')
        
    return

def make_clustermaps(featZ, meta, featsets, strain_lut, feat_lut, saveto, group_vars=['worm_gene','imaging_date_yyyymmdd']):
    """
    Parameters
    ----------
    featZ : TYPE
        DESCRIPTION.
    meta : TYPE
        DESCRIPTION.
    featsets : TYPE
        DESCRIPTION.
    strain_lut : TYPE
        DESCRIPTION.
    feat_lut : TYPE
        DESCRIPTION.
    saveto : TYPE
        DESCRIPTION.
    group_vars : TYPE, optional
        DESCRIPTION. The default is ['worm_gene','imaging_date_yyyymmdd'].

    Returns
    -------
    clustered_features : TYPE
        DESCRIPTION.

    """
    featZ_grouped = pd.concat([featZ,
                               meta],
                              axis=1
                              ).groupby(group_vars).mean()
    featZ_grouped.reset_index(inplace=True)
    
    row_colors = featZ_grouped['worm_gene'].map(strain_lut)
    
    # make clustermaps
    clustered_features = {}
    
    for stim, fset in featsets.items():
        col_colors = featZ_grouped[fset].columns.map(feat_lut)
        cg = sns.clustermap(featZ_grouped[fset],
                        row_colors=row_colors,
                        col_colors=col_colors,
                        vmin=-2,
                        vmax=2
                        )
        cg.ax_heatmap.axes.set_xticklabels([])
        cg.ax_heatmap.axes.set_yticklabels([])
        if saveto!=None:          
            cg.savefig(Path(saveto) / '{}_clustermap.png'.format(stim))
        
        clustered_features[stim] = np.array(featsets[stim])[cg.dendrogram_col.reordered_ind]
        plt.close('all')
        
    return clustered_features

def make_barcode(heatmap_df, selected_feats, cm=['inferno', 'inferno', 'Greys', 'Pastel1'], vmin_max = [(-2,2), (-2,2), (0,20), (1,3)]):
    """

    Parameters
    ----------
    heatmap_df : TYPE
        DESCRIPTION.
    selected_feats : TYPE
        DESCRIPTION.
    cm : TYPE, optional
        DESCRIPTION. The default is ['inferno', 'inferno', 'gray', 'Pastel1'].
    vmin_max : TYPE, optional
        DESCRIPTION. The default is [(-2,2), (-2,2), (-20, 0), (1,3)].

    Returns
    -------
    f : TYPE
        DESCRIPTION.

    """
    from matplotlib.gridspec import GridSpec
    font_settings = {'fontsize':14}
    
    f = plt.figure(figsize= (20,5))
    gs = GridSpec(4, 1,
                  wspace=0,
                  hspace=0,
                  height_ratios=[3,3,1,1])
    cbar_ax = f.add_axes([.91, .3, .03, .4])

    for n, ((ix,r), c, v) in enumerate(zip(heatmap_df.iterrows(), cm, vmin_max)):  
        axis = f.add_subplot(gs[n])
        sns.heatmap(r.to_frame().transpose().astype(float),
                    yticklabels=[ix],
                    xticklabels=[],
                    ax=axis,
                    cmap=c,
                    cbar=n==0, #only plots colorbar for first plot
                    cbar_ax=None if n else cbar_ax,
                    vmin=v[0],
                    vmax=v[1])
        axis.set_yticklabels(labels=[ix],
                             rotation=0,
                             fontsize=18)
        
        if n>2:
            c = sns.color_palette('Pastel1',3)
            sns.heatmap(r.to_frame().transpose(),
                    yticklabels=[ix],
                    xticklabels=[],
                    ax=axis,
                    cmap=c,
                    cbar=n==0, 
                    cbar_ax=None if n else cbar_ax,
                    vmin=v[0],
                    vmax=v[1])
            axis.set_yticklabels(labels=[ix],
                                 rotation=0,
                                 fontsize=18)
        cbar_ax.set_yticklabels(labels = cbar_ax.get_yticklabels(), fontdict=font_settings)
        f.tight_layout(rect=[0, 0, .9, 1])

    for sf in selected_feats:
        try:
            axis.text(heatmap_df.columns.get_loc(sf), 1, '*', fontdict=font_settings)
        except KeyError:
            print('{} not in featureset'.format(sf))
    return f

def make_heatmap_df(fset, featZ, meta, p_vals):
    """

    Parameters
    ----------
    fset : TYPE
        DESCRIPTION.
    featZ : TYPE
        DESCRIPTION.
    meta : TYPE
        DESCRIPTION.
    p_vals : TYPE
        DESCRIPTION.

    Returns
    -------
    heatmap_df : TYPE
        DESCRIPTION.

    """
    heatmap_df = [pd.concat([featZ,
                        meta],
                       axis=1
                       ).groupby('worm_gene').mean()[fset]]
                                                         
    heatmap_df.append(-np.log10(p_vals[fset])) 
    
    _stim = pd.DataFrame(data=[i.split('_')[-1] for i in fset],
                         columns=['stim_type'])
    _stim['stim_type'] = _stim['stim_type'].map(STIMULI_ORDER)
    _stim = _stim.transpose()
    _stim.rename(columns={c:v for c,v in enumerate(fset)}, inplace=True)
    heatmap_df.append(_stim)
    
    heatmap_df = pd.concat(heatmap_df)
    return heatmap_df

def clustered_barcodes(clustered_feats_dict, selected_feats, featZ, meta, p_vals, saveto):
    """
    
    Parameters
    ----------
    clustered_feats_dict : TYPE
        DESCRIPTION.
    selected_feats : TYPE
        DESCRIPTION.
    featZ : TYPE
        DESCRIPTION.
    meta : TYPE
        DESCRIPTION.
    p_vals : TYPE
        DESCRIPTION.
    saveto : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    for stim, fset in clustered_feats_dict.items():
        
        heatmap_df = make_heatmap_df(fset, featZ, meta, p_vals)
        
        f = make_barcode(heatmap_df, selected_feats)
        f.savefig(saveto / '{}_heatmap.png'.format(stim))
    return

def feature_box_plots(feature, feat_df, meta_df, bhP_values_df, strain_lut, add_stats=True):
    """

    Parameters
    ----------
    feature : TYPE
        DESCRIPTION.
    feat_df : TYPE
        DESCRIPTION.
    meta_df : TYPE
        DESCRIPTION.
    bhP_values_df : TYPE
        DESCRIPTION.
    strain_lut : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    from statannot import add_stat_annotation
    label_format = '{0:.4g}'  
    
    plt.figure(figsize=(5,10))
    ax = sns.boxplot(y=feature,
                x='worm_gene',
                data=pd.concat([feat_df, meta_df],
                               axis=1),
                order=strain_lut.keys(),
                palette=strain_lut.values())
    ax.set_ylabel(fontsize=20, ylabel=feature)
    ax.set_yticklabels(labels=[label_format.format(x) for x in ax.get_yticks()], fontsize=14) #labels = ax.get_yticks(), 
    ax.set_xlabel('')
    ax.set_xticklabels(labels = '')
    if add_stats:
        add_stat_annotation(ax,
                            data=pd.concat([feat_df,
                                            meta_df], axis=1),
                            x='worm_gene',
                            y=feature,
                            order=strain_lut.keys(),
                            box_pairs=[strain_lut.keys()],
                            perform_stat_test=False,
                            pvalues=[bhP_values_df[feature].values[0]],
                            test=None,
                            text_format='star',
                            loc='inside',
                            verbose=2,
                            fontsize=18)     
    plt.tight_layout()
    
    return

def window_errorbar_plots(feature, feat, meta, cmap_lut, window_order, plot_legend=False):
    import matplotlib.patches as patches
    from textwrap import wrap
    label_format = '{0:.4g}' 
    plt.figure()
    ax = sns.pointplot(x='window',
                        y=feature,
                        data=pd.concat([feat, meta], axis=1),
                        order=window_order,
                        hue='worm_gene',
                        showfliers=False,
                        linestyle=['-'],
                        ci=95,
                        alpha=0.8,
                        palette=cmap_lut)
    if not plot_legend:
        ax.legend().remove()
    ax.set_ylabel(fontsize=14,
                  ylabel='\n'.join(wrap(feature, 30)))
    ax.set_yticklabels(labels=[label_format.format(x) for x in ax.get_yticks()],
                       fontsize=14)
    ax.set_xlabel(fontsize=14,
                  xlabel='window')
    ax.set_xticklabels(labels=[BLUELIGHT_WINDOW_DICT[x] for x in ax.get_xticks()],
                       fontsize=12,
                       rotation=45)
    y_min = ax.axes.get_ylim()[0]
    y_max = ax.axes.get_ylim()[1]
    if y_min<y_max:
        rects = (patches.Rectangle((0.5, y_min), 1, abs(y_min-y_max), facecolor='b', alpha=0.3),
                 patches.Rectangle((3.5, y_min), 1, abs(y_min-y_max), facecolor='b', alpha=0.3),
                 patches.Rectangle((6.5, y_min), 1, abs(y_min-y_max), facecolor='b', alpha=0.3))

    [ax.add_patch(r) for r in rects]

    plt.tight_layout()
    return