# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 10:08:35 2018

@author: khaik
"""

import numpy as np
import pickle
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

class UniOfCam:
    def __init__(self, refStr, obsList, 
                 gapOpen=-10, gapExtend=-0.5, reload=True):
        self.rnseq = refStr
        self.onseq = obsList

        self.gapOpen = gapOpen
        self.gapExtend = gapExtend
        self.reload = reload
        
        self.alignedNT = []

    def setup(self):        
        if self.reload:
            with open('dictionaries/pdict','rb') as f:
                self.pdict = pickle.load(f)
        
        self.rpseq = self.translate(self.rnseq)
        
    def okarin(self,i,oldseq):    
        pseq = self.translate(oldseq)
        if '*' in pseq[:-1]:
            return
        
        opInsertion = []
        score = []
        
        alignments = pairwise2.align.globalds(self.rpseq[:-1],pseq[:-1],blosum62,self.gapOpen,self.gapExtend)
        newseq_all = []
        
        for a in alignments:
            newseq=''
            k=0
            for j in range(len(a[0])):
                if (a[1][j]!='-') and (a[0][j]!='-'):
                    newseq += oldseq[k*3:k*3+3]
                    k+=1
                elif (a[1][j]=='-') and (a[0][j]!='-'):
                    newseq += '---'
                elif (a[0][j]=='-') and (a[1][j]!='-'):
                    opInsertion.append(k*3)
                    k+=1
            newseq += oldseq[-3:]
            newseq_all.append(newseq)
            score.append(a[2])
            
        self.alignedNT.append({'index':i,'sequence':newseq_all,
                               'insertion':opInsertion,'score':score})
    
    def matchAA(self,codon):
        aa = None
        for key,value in self.pdict.items():
            if codon in value:
                aa = key
        return aa

    def translate(self,seq):
        pseq = ''
        nAA = int(np.floor(len(seq)/3))
        for i in range(nAA):
            pseq += self.matchAA(seq[i*3:i*3+3])
        return pseq
    
    def format_alignment(self,alignment,view=False,size=75):
        lineThings = ''
        for i in range(len(alignment[0])):
            if alignment[0][i]==alignment[1][i]:
                lineThings+='|'
            else:
                lineThings+='-'
        FA = (alignment[0],lineThings,alignment[1])
        
        if view:
            for element in FA:
                print(element[:size])
        
        return FA
    
    def checkSeq(self, seq):
        if len(seq)%3!=0:
            return False
        for nt in seq:
            if nt not in ['T','C','A','G']:
                return False
        return True
        
        
    def main(self):
        self.setup()
        
        rlen = len(self.rnseq)
        for i in range(len(self.onseq)):
            seq = self.onseq[i]

            if not self.checkSeq(seq):
                continue
            
            if (len(seq)==rlen):
                self.alignedNT.append({'index':i,'sequence':[self.onseq[i]],
                                       'insertion':[],'score':[]})        
            else:
                self.okarin(i,seq)
            
            if i%100==0:
                print(i,end=' ')
