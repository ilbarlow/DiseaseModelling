#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:29:12 2021

@author: ibarlow

Script for making disease model metadata for the data that Tom collected

"""

import pandas as pd
from pathlib import Path
import re

from tierpsytools.hydra.compile_metadata import populate_96WPs,\
    get_day_metadata, concatenate_days_metadata, day_metadata_check, \
        number_wells_per_plate

date_regex = r"\d{8}"

ROOT_DIR = Path('/Volumes/behavgenom$/Tom/Data/Hydra/DiseaseModel/RawData')
AUX_DIR = ROOT_DIR / 'AuxiliaryFiles'

DRUG_PLATES = Path('/Volumes/behavgenom$/Tom/Data/Hydra/DiseaseModel/RawData/AuxiliaryFiles/manual_metadata/source_plates.xlsx')

# %% 
if __name__=='__main__':
    day_root_dirs = [d for d in AUX_DIR.glob("*")
                     if d.is_dir() and re.search(date_regex, str(d))
                     is not None]
    
    print('Calculating metadata for {} days of experiments'.format(
            len(day_root_dirs)))
    
    drug_plates = pd.read_excel(DRUG_PLATES)
    
    for count, day in enumerate(day_root_dirs):
        exp_date = re.findall(date_regex, str(day))[0]
        manualmetadata_file = list(day.rglob('*_manual_metadata.csv'))[0]
        assert (exp_date in str(manualmetadata_file))
        wormsorter_file = list(day.rglob('*_wormsorter.csv'))[0]
        assert (exp_date in str(wormsorter_file))
        # bad_wells_file = list(day.rglob('*_robot_bad_imaging_wells.csv'))[0]
        # assert (exp_date in str(bad_wells_file))
        
        print('Collating wormsorter files: {}'.format(wormsorter_file))
        
        plate_metadata = populate_96WPs(wormsorter_file,
                                        del_if_exists=True,
                                        saveto='default')
        # print(day)
        # print(number_wells_per_plate(plate_metadata, day))
        
        # bad_wells_df = convert_bad_wells_lut(bad_wells_file)
        
        # plate_metadata = pd.merge(plate_metadata,
        #                           bad_wells_df,
        #                           on=['imaging_plate_id', 'well_name'],
        #                           how='outer')
        # plate_metadata['is_bad_well'].fillna(False,
        #                                      inplace=True)
        metadata_file = day / '{}_day_metadata.csv'.format(exp_date)
    
        print('Generating day metadata: {}'.format(
                metadata_file))
        
        try:
            day_metadata = get_day_metadata(plate_metadata,
                                            manualmetadata_file,
                                            saveto=metadata_file,
                                            del_if_exists=True,
                                            include_imgstore_name=True)
        except ValueError:
            print('imgstore error')
            day_metadata = get_day_metadata(plate_metadata,
                                            manualmetadata_file,
                                            saveto=metadata_file,
                                            del_if_exists=True,
                                            include_imgstore_name=False)
        
        day_metadata['source_plate_id'] = day_metadata['imaging_plate_id'
                                                       ].apply(lambda x: '_'.join(x.split('_')[:-1]) if len(x.split('_'))>3 else x)
        day_metadata = pd.merge(day_metadata,
                               drug_plates,
                               on=['source_plate_id',
                                        'well_name'],
                               suffixes=('_day', '_robot'),
                               how='outer')
        day_metadata.drop(
            day_metadata[day_metadata['imaging_plate_id'].isna()].index,
                        inplace=True)
           
        files_to_check = day_metadata_check(day_metadata, day, plate_size=48)
        number_wells_per_plate(day_metadata, day)
        
        day_metadata.to_csv(metadata_file, index=False)
#%%        
    import datetime
    # combine all the metadata files
    concat_meta = concatenate_days_metadata(AUX_DIR,
                                            list_days=None,
                                            saveto=None)
    
    
    concat_meta_grouped = concat_meta.groupby('worm_gene')

    strains = pd.DataFrame(concat_meta_grouped.apply(lambda x: x.drop_duplicates(subset='worm_strain')))
    strains.reset_index(drop=True,
                        inplace=True)
    
    strains = strains[['worm_gene',
                        'worm_strain',
                        'worm_code',
                        'date_yyyymmdd']]
 
    strains.to_csv(AUX_DIR/ \
                   '{}_strain_name_errors.csv'.format(
                       datetime.datetime.today().strftime('%Y%m%d')
                       ),
                   index=False)
        