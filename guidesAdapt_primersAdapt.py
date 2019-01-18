""" create function to add on 3' and 5' adapters for ordering guideRNA oligos
and to add adapters for the MiSeq primers"""

import sys
import pandas as pd

def guidesAdapt_primersAdapt(input_sgRNAs, outputfile):
    """ adapts the guides according to their 5' end:
    T7 polymerase if GG and
    SP6 polymerase if GA

    Input
    input_sgRNAs - .csv with columns 'target_site', 'fwd_primer', 'rev_primer',
    'size', 'Tm', 'Type', 'sgRNA_no'
    
    Output
    outputfile - .csv containing sgRNAs and MiSeq primers"""
    
    inputDF = pd.read_csv(input_sgRNAs)

    SP6promoter = 'ATTTAGGTGACACTATA'
    T7promoter = 'TAATACGACTCACTATA'
    
    crRNA = 'GTTTTAGAGCTAGAAATAGCAAG'
    
    MiSeqFwd = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
    MiSeqRev = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
   
    oligos = []
    for i,r in inputDF.iterrows():
        #first check the guide is the correct length 
        if len(r.target_site.upper()) == 23:
            print (r.gene + ' ' + r.target_site.upper() + ' guide correct length')
        else:
            print(r.gene + ' guide incorrect length!')
            continue
        
        if r.target_site.upper()[0:2] == 'GA':
            sgRNA = SP6promoter + r.target_site.upper()[:-3] + crRNA
            promoter = 'SP6'
        elif r.target_site.upper()[0:2] == 'GG':
            sgRNA = T7promoter + r.target_site.upper()[:-3] + crRNA
            promoter = 'T7'
        else:
            sgRNA = SP6promoter + r.target_site.upper()[:-3] + crRNA
            promoter = 'SP6'

        fwd_primer = MiSeqFwd + r.fwd_primer.upper()
        rev_primer = MiSeqRev + r.rev_primer.upper()
       
        oligos.append((r.gene.replace('-', '') + '_' + promoter, sgRNA, fwd_primer, rev_primer))
   
    #make into a dataframe to mave a two-column csv
    names = ['_guide', '_MiSeqFwd', '_MiSeqRev']
    oligo_order = pd.DataFrame()
    for g in range(1,4):
        temp = pd.DataFrame()
        temp['sequence'] = [i[g] for i in oligos]
        temp ['oligo_name'] = [i[0]+names[g-1] for i in oligos]
        oligo_order = oligo_order.append(temp)
        del temp
        
    oligo_order= oligo_order.reset_index(drop=True)
    
    #save to csv
    oligo_order.to_csv(outputfile, index=False)
    
if __name__ == '__main__':
    input_sgRNAs = sys.argv[1]
    outputfile = sys.argv[2]
    
    guidesAdapt_primersAdapt(input_sgRNAs, outputfile)