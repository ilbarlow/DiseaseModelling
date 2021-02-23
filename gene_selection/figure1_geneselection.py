#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 09:11:51 2021

@author: ibarlow

Script for making venn diagrams for disease model paper
"""

import pandas as pd
import numpy as np
from pathlib import Path
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
import sys

sys.path.insert(0, '/Users/ibarlow/Documents/GitHub/pythonScripts/DiseaseModelling/hydra_screen/phenotype_summary')

from plotting_helper import CUSTOM_STYLE
plt.style.use(CUSTOM_STYLE)

ORTHOLIST_FNAME = Path('/Users/ibarlow/Documents/MATLAB/WormDiseaseModelling/ortholist_master.txt')
WORMINE_FNAME = Path('/Users/ibarlow/Documents/MATLAB/WormDiseaseModelling/wormbase_simplemine_results.txt')

WORM_KEYWORDS = ['neural', 'neuron', 'muscle', 'muscular']

DISEASE_KEYWORDS = ['Epilepsy', 'Epileptic', 'Autism', 'Rett',
    'mental', 'McArdle', 'Cerebral palsy', 'Schizophrenia',
    'Major depressive disorder', 'personality', 'Bipolar', 'ADHD',
    'Attention deficit-hyperactivity disorder', 'seizure', 'parkinson',
   'dopamin', 'serotonin', '5-HT',
    'GPCR', 'antidepressant', 'antipsychotic']

SAVETO = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/behavgenom_copy/DiseaseScreen/summary_figures')

#%%
if __name__ == '__main__':
    ortho_df = pd.read_table(ORTHOLIST_FNAME)
    wormine_df = pd.read_table(WORMINE_FNAME)
    
    combi_df = pd.merge(wormine_df,
                        ortho_df,
                        left_on='Your Input',
                        right_on='WormBase ID',
                        how='outer')
    # plot number of programs in agreemen 
    combi_df.groupby('No. of Programs')['OMIM Phenotypes'].count().plot.bar()
    plt.ylabel('Count')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(SAVETO / 'OMIM_number_of_programs.png', dpi=200)
    
    #only use genes with agreement >=2
    true_ortho = combi_df.query('`No. of Programs` >=2')
    
    #only take genes with OMIM phenotype 
    phenotypeFlag = combi_df['OMIM Phenotypes'].notna()
    true_omim = combi_df[phenotypeFlag]

    # look in the wormbase desciptive test for worm keywords
    worm_kwords = []
    for k in WORM_KEYWORDS:
        worm_kwords.append(combi_df[combi_df["Description Text"].str.contains(k,
                                                                              na=False)])

    worm_kwords = pd.concat(worm_kwords).drop_duplicates()
  
    # make a venn diagram of overlap of omim and ortholist selection criteria  
    plt.figure()
    venn2(subsets = (true_ortho['WormBase Gene ID'].unique().shape[0],
                     true_omim['WormBase Gene ID'].unique().shape[0],
                     len(set(true_ortho['WormBase Gene ID'].unique()).intersection(true_omim['WormBase Gene ID'].unique()))),
                     set_labels = ('Ortholist,\n number of programs >= 2',
                                'OMIM diseases'),
                     set_colors = ('m', 'g'),
                     alpha=0.3)
    plt.title('Number of worm genes')
    plt.tight_layout()
    plt.savefig(SAVETO / 'Ortholist_OMIM_venn.png', dpi=200)
    
    #make figure of overlap with all three selection criteria
    plt.figure()
    venn3(subsets = (true_ortho['WormBase Gene ID'].unique().shape[0],
                     true_omim['WormBase Gene ID'].unique().shape[0],
                      len(set(true_ortho['WormBase Gene ID'].unique()).intersection(true_omim['WormBase Gene ID'].unique())),
                      worm_kwords['WormBase Gene ID'].unique().shape[0],
                      len(set(true_ortho['WormBase Gene ID'].unique()).intersection(worm_kwords['WormBase Gene ID'].unique())),
                      len(set(true_omim['WormBase Gene ID'].unique()).intersection(worm_kwords['WormBase Gene ID'].unique())),
                       len(set(true_ortho['WormBase Gene ID'].unique()).intersection(true_omim['WormBase Gene ID'].unique()).intersection(worm_kwords['WormBase Gene ID'].unique()))
                      ),
                  set_labels = ('Ortholist,\n number of programs >= 2',
                                'OMIM phenotype',
                                'Wormbase neuro/musc gene'),
                  alpha = 0.3,
                  set_colors=('m', 'g', 'b'))
    plt.tight_layout()
    plt.savefig(SAVETO / 'ortholist_OMIM_wormbase_venn.png', dpi=200)
    plt.savefig(SAVETO / 'ortholist_OMIM_wormbase_venn.svg', dpi=200)
    
    #get indices of intersection of all three
    selected_wbids = set(true_ortho['WormBase Gene ID'].unique()).intersection(true_omim['WormBase Gene ID'].unique()).intersection(worm_kwords['WormBase Gene ID'].unique())
    
    first_filter = combi_df.query('@selected_wbids in `WormBase Gene ID`')
    
    print('Number of unique worm genes: {}'.format(first_filter['WormBase Gene ID'].unique().shape[0]))
    print('Number of unique human genes: {}'.format(first_filter['HGNC Symbol'].unique().shape[0]))
    
    #%%
    omim_kwords = []
    for dk in DISEASE_KEYWORDS:
        omim_kwords.append(true_omim[true_omim['OMIM Phenotypes'].str.contains(dk,
                                                                     na=False)])
    omim_kwords = pd.concat(omim_kwords).drop_duplicates()
    # ortho_omim = set(true_ortho.index) - set(true_omim.index)
