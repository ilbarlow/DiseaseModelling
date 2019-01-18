#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 14:06:19 2019

@author: ibarlow
"""

""" Create features files to annotate in Ape"""

import os
import sys
import pandas as pd


def FeaturesGen(ChopChopresults, outputDir, sgRNA_type):
    """ function to aggregate the sgRNA search results
    Input:
        ChopChopresult - results folder generated from ChopChop
        outputDir - directory into which the features.txt file should be saved
        sgRNA_type - 'General'/'GG'/'GA' to determine colouring of features
            General: pink + green
            GA: cyan + cornflower blue
            GG: yellow + plum
    Output:
        features.txt file containin the tab delimited features for importing into Ape"""
    
    #make output Directory if it does not already exist
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)
        
    #list the directory contents    
    for i,j,k in os.walk(ChopChopresults): #use walk to go through and find all directories
        
        if j == []: #no subdirectories
            saveDF = pd.DataFrame() #initiate dataframe
            for target in k: #loop through to find the sgRNA sequences
                if target.endswith('.offtargets'):
                    with open(os.path.join(i,target), 'r+') as f:
                        guide = f.readlines()
                    #add them to a dataframe
                    temp = pd.Series()
                    temp['guideNo'] = target.split('.')[0] + sgRNA_type
                    temp['guideSeq'] = guide.pop(0).rstrip()
                    
                    saveDF = saveDF.append(temp.to_frame().transpose())
            saveDF['type'] = 'sgRNA'
            
            if sgRNA_type == 'General' or sgRNA_type == None:
                saveDF['fwd'] = 'pink'
                saveDF['rev'] = 'green'
            elif sgRNA_type == 'GG':
                saveDF['fwd'] = 'yellow'
                saveDF['rev'] = 'plum'
            elif sgRNA_type == 'GA':
                saveDF['fwd'] = 'cyan'
                saveDF['rev'] = 'cornflower blue'
                
            
            #save to txt file with tab delimiter
            saveDF.to_csv(os.path.join(outputDir, os.path.basename(i) + '_features.txt'),\
                          index = False, header = False, sep = '\t')
        
            del saveDF

if __name__  == '__main__':
    ChopChopresults = sys.argv[1]
    outputDir = sys.argv[2]
    sgRNA_type = sys.argv[3]
    
    FeaturesGen(ChopChopresults, outputDir, sgRNA_type)
        
        