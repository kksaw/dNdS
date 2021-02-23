# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 17:08:20 2018

@author: khaik
"""

import numpy as np
import pickle
from itertools import islice
from Bio import SeqIO, Entrez
import re
import os

class getSeq:
    #output = [orgn, ID, seq, len]
    def __init__(self, filename, source=1, method=0, gene=None, size=None,
                 exclude=None, excludeFlag=False, excludeOrgn=None,
                 geneRange=(0,0), geneRange2=(0,0),
                 roughGene=None, nGene=1, keys=None):
        
        self.filename = filename
        self.source = source
        self.method = method
        
        self.geneName = []
        self.size = size
        self.exclude = exclude
        self.excludeFlag = excludeFlag
        self.excludeOrgn = []
        self.geneRange = geneRange
        self.geneRange2 = geneRange2
        self.roughGene = roughGene
        
        self.nGene = nGene
        self.keys = keys
        
        if self.method==0:
            for g in gene:
                self.geneName.append('/gene="'+g+'"')
        
        elif self.method==1:
            for g in gene:
                self.geneName.append(g)
        
        elif self.method==2:
            for g in gene:
                self.geneName.append(g)
            self.roughGene = []
            for g in roughGene:
                self.roughGene.append(g)
        
        if excludeOrgn:
            for orgn in excludeOrgn:
                self.excludeOrgn.append(orgn)
        
        if self.keys==None:
            self.keys = [str(i) for i in range(nGene)]
            
        self.gene_list = {}
        for i in range(nGene):
            self.gene_list.update({self.keys[i]:[]})
            
        #placeholders
        self.currID = ''
        self.currOrgn = ''
        self.currGene = ''
        self.currRegLine = ''
        self.currNameLine = ''
        
    def setup(self):
        Entrez.email= 'kks36@cam.ac.uk'

    def init_data_output(self):
        import time
        startTime = time.localtime()
        filename = 'genelist' + '-' + geneName + '-' + time.strftime('%y-%m-%d-%H-%M', startTime)
        
    def readNCBI0(self):
        with open(self.filename+'.txt','r') as f:
            data = f.readlines()

        for line in data:
            if '[' in line and ']' in line:
                orgn = line[line.find('[')+1:line.find(']')]
            elif 'Annotation: ' in line:
                geneID = line[line.find('NC'):line.find('(')-1]
                m = re.search('\d',line).start()
                geneStart = int(line[m:line.find('..')]) 
                m = re.search('\d',line[::-1]).start()
                geneEnd = int(line[line.find('..')+2:-m])
                dirn = 2 if 'complement' in line else 1
                
                with Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', 
                                   id=geneID, seq_start=geneStart, seq_stop=geneEnd, 
                                   strand=dirn) as handle:
                    new_seq = SeqIO.read(handle,'fasta').seq._data
                    len_seq = len(new_seq)
                
                self.gene_list.append([orgn,geneID,new_seq,len_seq])
            else:
                continue
    
    def readNCBI1(self):
        with open(self.filename,'r') as f:
            IDs = list(islice(f,self.size[0],self.size[1]))
        for line in IDs:
            self.currID = line[:-1]
            
            if self.excludeFlag:
                if self.currID.startswith(any(x for x in self.exclude)):
                    continue
            
            handle = Entrez.efetch(db='nucleotide', rettype='gbwithparts', retmode='text',
                                   id=self.currID)
            if self.method==0:
                self.byGene(handle)
            elif self.method==1:
                self.byProduct_short(handle)
            elif self.method==2:
                self.byProduct_long(handle)
                
            print(self.currID, end=' ')
    
    def byGene(self,handle):
        prevLine = ''
        for currLine in handle.readlines():
            if 'ORGANISM  ' in currLine:
                self.currOrgn = currLine[currLine.find('ORGANISM')+10:currLine.find('\n')]
                continue
        
            elif '     gene            ' in prevLine:
                self.currRegLine = prevLine
                self.currNameLine = currLine
                self.fetch('/gene="')
                pass
            
            prevLine = currLine

    def byProduct_short(self,handle):
        for currLine in handle.readlines():
            if 'ORGANISM  ' in currLine:
                self.currOrgn = currLine[currLine.find('ORGANISM')+10:currLine.find('\n')]
            
                if (any(x in self.currOrgn for x in self.excludeOrgn)):
                    return
                
                continue
        
            elif '     gene            ' in currLine:
                self.currRegLine = currLine
                continue
        
            elif ('/product="' in currLine):
                self.currNameLine = currLine
                self.fetch('/product="')
                

    def byProduct_long(self,handle):
        currLine0=''
        for currLine in handle.readlines():
            if 'ORGANISM  ' in currLine:
                self.currOrgn = currLine[currLine.find('ORGANISM')+10:currLine.find('\n')]
                continue
        
            elif '     gene            ' in currLine:
                self.currRegLine = currLine
                continue
        
            elif ('/product="' in currLine) and (any(x in currLine.lower() for x in self.roughGene)):
                if currLine.endswith('"\n'):
                    self.currNameLine = currLine
                    currLine0 = ''
                else:
                    currLine0 = currLine
                    continue
            
            elif currLine0:
                self.currNameLine = currLine0[:-1] + currLine[20:]
                currLine0 = ''
                    
            if self.currNameLine:
                self.fetch('/product="')
                self.currNameLine = ''        
        
    def fetch(self,label):                        
        if not self.currRegLine:
            return
        
        if (any(x in self.currNameLine.lower() for x in self.geneName)):
            if ('<' in self.currRegLine) or ('>' in self.currRegLine):
                return
            else:
                self.currGene = self.currNameLine[self.currNameLine.find(label)+len(label):self.currNameLine.find('"\n')]
        else:
            return
        
        m = re.search('\d',self.currRegLine).start()
        geneStart = int(self.currRegLine[m:self.currRegLine.find('..')])
        m = re.search('\d',self.currRegLine[::-1]).start()                
        geneEnd = int(self.currRegLine[self.currRegLine.find('..')+2:-m])
        geneLen = geneEnd-geneStart
        
        if self.geneRange[0]<geneLen<self.geneRange[1]:
            my_key=self.keys[0]
        elif self.geneRange2[0]<geneLen<self.geneRange2[1]:
            my_key=self.keys[1]
        else:
            return
        
        dirn = 2 if 'complement' in self.currRegLine else 1
        
        with Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', 
               id=self.currID, seq_start=geneStart, seq_stop=geneEnd, 
               strand=dirn) as handle:
            new_seq = SeqIO.read(handle,'fasta').seq._data
            len_seq = len(new_seq)
        
        data = {'organism':self.currOrgn,
                'accession':self.currID,
                'sequence':new_seq,
                'length':len_seq,
                'gene':self.currGene}
        
        self.gene_list[my_key].append(data)
        print('*',end='')


    def main(self):
        self.setup()
        if self.source==0:
            self.readNCBI0()
        elif self.source==1:
            self.readNCBI1()


        
#def main(ref, acsList='ID_pseudomonas4360.seq', source=1, geneList=['isocitrate dehydrogenase']):
#    
#    step1 = getSeq('ID_pseudomonas4360.seq', source=1, method=1, gene=['isocitrate dehydrogenase'], 
#                   size=(0,None), geneRange=(1150,1350), geneRange2=(2150,2350), excludeOrgn=['aeruginosa'])
#    step1.main()
#    gene_list = step1.gene_list
#    sequence_list = []
#    for line in gene_list:
#        sequence_list.append(line['sequence'])
#    print('step1',end=' ')
#    
#    ID_step1 = getSeq('ID_pseudomonas4360.seq', source=1, method=1, gene=['isocitrate dehydrogenase'], 
#               size=(42,None), geneRange=(1150,1350), geneRange2=(2150,2350),nGene=2)
#    sequence_list = []
#    index_list = []
#    for i in ID_step1.gene_list:
#        if 'aeruginosa' not in line['organism']:
#            sequence_list.append(line['sequence'])
#            index_list.append(i)
#   step1 = getSeq('icl_pseudomonas_xAE4212.seq',method=1,gene=['isocitrate lyase'],size=(0,None),geneRange=(1280,1400),geneRange2=(1550,1650),nGene=2)