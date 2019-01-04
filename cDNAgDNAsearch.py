#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 11:16:16 2019

@author: ibarlow
"""

""" Go through the chopchop outputs and map them on to the cDNA and gDNA
sequences by adding a feature file"""

import os
from collections import defaultdict
import re

def find_sgRNAs(inputcDNA, inputgDNA, ChopChopresults, outputFeatures):
    """Input
    inputcDNA - cDNA.ape sequence
    inputgDNA - gDNA.fa sequence
    ChopChopResults - output folder from running Chopchop
    
    Output
    outputFeatures - .txt file of features to add on to the DNA sequences
    
    """


    DNA_COMPLEMENTS = defaultdict(lambda: 'x', str.maketrans('AaTtGgCc', 'TtAaCcGg'))
    
    with open (inputcDNA, 'r+') as fopen:
        cDNA = fopen.read()
        
    with open(inputgDNA, 'r+') as fopen:
        gDNA = fopen.readlines()
        #remove the first line
        gDNA.pop(0)
        #make into one long string
        gDNA2 = gDNA.pop(0).rstrip()
        for i in gDNA:
            gDNA2 += i.rstrip()
            
    chopchopDir = os.listdir(ChopChopresults)
    guideSites_cDNA = []
    guideSites_gDNA = []
    for i in chopchopDir:
        if i.endswith('.offtargets'):
            with open(os.path.join(ChopChopresults,i), 'r+') as f:
                guide = f.readlines()
            guide1 = guide.pop(0).rstrip()
            
            if re.search(guide1, cDNA.upper()):    
                guideSites_cDNA.append(re.search(guide1, cDNA.upper()).span())
                guideSites_cDNA.append('cyan')
            else:
                cDNA_rev= cDNA.translate(DNA_complements)
                try:
                    guideSites_cDNA.append(re.search(guide1, cDNA_rev).span())
                    guideSite_cDNA.append('green')
                
            
            
                
            
