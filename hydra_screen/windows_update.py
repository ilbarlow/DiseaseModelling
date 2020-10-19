#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:40:02 2020

@author: ibarlow

prepare windows data for analysis
"""

import pandas as pd
from pathlib import Path
import re

PROJECT_DIR = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen/Results')

windows_files = list(PROJECT_DIR.rglob('*window_*'))

fname_files = [f for f in windows_files if 'filename' in (str(f))]
feat_files = [f for f in windows_files if 'features' in str(f)]

#%%
def find_window(fname):
    window_regex = r"(?<=_window_)\d{0,9}"
    window = int(re.search(window_regex, str(fname))[0])
    return window
   
#%%
if __name__ == "__main__":
    # make tuples of the paired files
    feat_start = [f for f in feat_files if len(f.parts)>9]
    feat_start.sort(key=find_window)
    feat_end = [f for f in feat_files if len(f.parts)==9]
    feat_end.sort(key=find_window)
    feat_tuples = list(zip(feat_start, feat_end))
    
    fname_start = [f for f in fname_files if len(f.parts)>9]
    fname_start.sort(key=find_window)
    fname_end = [f for f in fname_files if len(f.parts)==9]
    fname_end.sort(key=find_window)
    fname_tuples = list(zip(fname_start, fname_end))

  #%%  
    for f in feat_tuples:
        assert find_window(f[0]) == find_window(f[1])
        
        with open(f[0],'a') as fd_start:
            df = pd.read_csv(f[1], comment='#',
                             header = None)   
            df.to_csv(fd_start, header=False, index=False)
            
    for f in fname_tuples:
        assert find_window(f[0]) == find_window(f[1])
        
        with open(f[0],'a') as fd_start:
            df = pd.read_csv(f[1],
                             comment='#',
                             header = None)            
            df.to_csv(fd_start,
                      header=False,
                      index=False)