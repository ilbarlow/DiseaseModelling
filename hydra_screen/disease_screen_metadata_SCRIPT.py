#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 09:33:02 2020

@author: ibarlow

Script for generating metadata for disease models

"""

import pandas as pd
from pathlib import Path
import re

from tierpsytools.hydra.compile_metadata import populate_96WPs,\
    merge_robot_metadata, merge_robot_wormsorter, get_day_metadata,\
    concatenate_days_metadata, day_metadata_check, number_wells_per_plate

date_regex = r"\d{8}"

PROJECT_DIRECTORY = Path('/Volumes/behavgenom$/Ida/Data/Hydra/DiseaseScreen')

#%%
if __name__ == '__main__':

    day_root_dirs = [d for d in (PROJECT_DIRECTORY /
                                 'AuxiliaryFiles').glob("*")
                     if d.is_dir() and re.search(date_regex, str(d))
                     is not None]

    print('Calculating metadata for {} days of experiments'.format(
            len(day_root_dirs)))

    for count, day in enumerate(day_root_dirs):
        exp_date = re.findall(date_regex, str(day))[0]
        manualmetadata_file = list(day.rglob('*_manual_metadata.csv'))[0]
        assert (exp_date in str(manualmetadata_file))
        wormsorter_file = list(day.rglob('*_wormsorter.csv'))[0]
        assert (exp_date in str(wormsorter_file))

        print('Collating manual metadata files: {}'.format(wormsorter_file))

        plate_metadata = populate_96WPs(wormsorter_file,
                                        del_if_exists=True,
                                        saveto='default')

#%%
        metadata_file = day / '{}_day_metadata.csv'.format(exp_date)

        print('Generating day metadata: {}'.format(
                metadata_file))

        day_metadata = get_day_metadata(plate_metadata,
                                        manualmetadata_file,
                                        saveto=metadata_file,
                                        del_if_exists=True)

   
        files_to_check = day_metadata_check(day_metadata, day, plate_size=48)
        number_wells_per_plate(day_metadata, day)

# %%
    # combine all the metadata files
    concatenate_days_metadata(PROJECT_DIRECTORY / 'AuxiliaryFiles',
                              list_days=None,
                              saveto=None)


