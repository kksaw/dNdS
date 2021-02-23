# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 17:08:20 2018

@author: khaik
"""

'''
GUI
- refseq
- seqList
- show alignment (default True) - show nMat, pMat (check)

- dirn (default v)
- windList (default [3,5,10])
- plot ()
- 


'''

import numpy as np
import pickle
from itertools import islice
from Bio import SeqIO, Entrez
import re
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

from Classes import Step1, Step2, Step3

pathRoot = 'C:\Users\khaik\Documents\camb\dNdS\icd'
refFile = 'icd_RefPA01'

acsFile = 'ID_pseudomonas4360.seq'
geneList = ['isocitrate dehydrogenase', 'icd']
gene1 = (1150,1350)
gene2 = (2150,2350)
nGene = 2
exclOrgn = 'aeruginosa'
    
def setup():
    with open(refFile,'rb') as f:
        refseq = pickle.load(f)
    

def main(exclude=False):
    step1 = Step1(acsFile, source=1, method=1, gene=geneList, 
                   size=(0,None), geneRange=gene1, geneRange2=gene2, nGene=nGene)
    step1.main()

    filename = pathRoot + r'\geneList'
    with open(filename,'wb') as f:
        pickle.dump(step1.gene_list,f)
    print('step1',end=' ')
    
    if exclude:
        sequence_list = []
        index_list = []
        for i in step1.gene_list:
            if exclOrgn not in step1.gene_list[i]['organism']:
                sequence_list.append(step1.gene_list[i]['sequence'])
                index_list.append(i)
    else:
        sequence_list = []
        for line in step1.gene_list:
            sequence_list.append(line['sequence'])
    
    
    step2 = Step2(refseq, sequence_list)
    step2.main()
    
    filename = pathRoot + r'\alignedNT'
    with open(filename,'wb') as f:
        pickle.dump(step2.alignedNT,f)
    print('step2',end=' ')
    
    sequence = []
    for line in step2.alignedNT:
        sequence.append(line['sequence'][0])
     
    step3 = Step3(refseq, sequence, dirn='v',
                 window=3)
    step3.main()
    
    filename = pathRoot + r'\values'
    with open(filename,'wb') as f:
        pickle.dump(step3.data,f)
    
    dnds = step3.data['dnds'][0][:-1].tolist()
    filename = pathRoot + r'\pbsdnds'
    with open(filename,'wb') as f:
        pickle.dump(dnds,f,protocol=2)

if __name__ == "__main__":
    main()