# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 17:00:23 2018

@author: khaik
"""
'''get NdSd for all codons - list of dictionaries (4096 values)'''

import numpy as np
import pickle
       
with open('dictionaries/cdict','rb') as f:
    cdict = pickle.load(f)
clist = list(cdict.keys())
    
with open('dictionaries/pdict','rb') as f:
    pdict = pickle.load(f)

data = [None]*64

def matchAA(codon):
    for key,value in pdict.items():
        if codon in value:
            return key
        
def compareCodon(c1,c2):
    n=0
    s=0
    if matchAA(c1)==matchAA(c2):
        s+=1
    else:
        n+=1
    return n,s
    
def compareCodons(d,c1,c2,c3,c4):
    if d==2:
        n1,s1 = compareCodon(c1,c3)
        n2,s2 = compareCodon(c3,c2)
        n = np.sum([n1,n2])
        s = np.sum([s1,s2])
    elif d==3:
        n1,s1 = compareCodon(c1,c3)
        n2,s2 = compareCodon(c3,c4)
        n3,s3 = compareCodon(c4,c2)
        n = np.sum([n1,n2,n3])
        s = np.sum([s1,s2,s3])
    return n,s

#brute force one codon at a time
def yikye(c1,c2):
    nd=0
    sd=0
    d=0
    dLoc=[]
    for k in range(3):
        if c1[k]!=c2[k]:
            dLoc.append(k)
            d+=1
        else:
            continue
        
    if d==1:
        nd,sd = compareCodon(c1,c2)

    elif d==2:  #2 permutations
        nd_list = [None]*2
        sd_list = [None]*2
        
        if dLoc[0]==0 and dLoc[1]==1:
            c3a=c2[0]+c1[1:]
            c3b=c1[0]+c2[1]+c1[2]
        elif dLoc[0]==1 and dLoc[1]==2:
            c3a=c1[0]+c2[1]+c1[2]
            c3b=c1[:2]+c2[2]
        elif dLoc[0]==0 and dLoc[1]==2:
            c3a=c1[:2]+c2[2]
            c3b=c2[0]+c1[1:]
        
        nd_list[0],sd_list[0] = compareCodons(d,c1,c2,c3a,None)
        nd_list[1],sd_list[1] = compareCodons(d,c1,c2,c3b,None)
        
        nd = np.mean(nd_list)
        sd = np.mean(sd_list)

    elif d==3:  #6 permutations
        nd_list = [None]*6
        sd_list = [None]*6
        
        c3a=c1[:2]+c2[2]
        c3b=c1[0]+c2[1]+c1[2]
        c3c=c2[0]+c1[1:]
        
        c4a=c1[0]+c2[1:]
        c4b=c2[0]+c1[1]+c2[2]
        c4c=c2[:2]+c1[2]
        
        nd_list[0],sd_list[0] = compareCodons(d,c1,c2,c3a,c4a)
        nd_list[1],sd_list[1] = compareCodons(d,c1,c2,c3a,c4b)
        nd_list[2],sd_list[2] = compareCodons(d,c1,c2,c3b,c4a)
        nd_list[3],sd_list[3] = compareCodons(d,c1,c2,c3b,c4c)
        nd_list[4],sd_list[4] = compareCodons(d,c1,c2,c3c,c4b)
        nd_list[5],sd_list[5] = compareCodons(d,c1,c2,c3c,c4c)
        
        nd += np.mean(nd_list)
        sd += np.mean(sd_list)
    
    else:
        pass
    
    return dLoc,nd,sd

clist0 = clist[:]
clist0.extend(('TT-','TC-','TA-','TG-',
               'CT-','CC-','CA-','CG-',
               'AT-','AC-','AA-','AG-',
               'GT-','GC-','GA-','GG-',
                      
               'T-T','T-C','T-A','T-G',
               'C-T','C-C','C-A','C-G',
               'A-T','A-C','A-A','A-G',
               'G-T','G-C','G-A','G-G',
                      
               '-TT','-CT','-AT','-GT',
               '-TC','-CC','-AC','-GC',
               '-TA','-CA','-AA','-GA',
               '-TG','-CG','-AG','-GG',
                      
               'T--','C--','A--','G--',
               '-T-','-C-','-A-','-G-',
               '--T','--C','--A','--G','---'))

#data - dloc,Nd,Sd
for i in range(64):
    data[i] = {'codon':clist[i]}
    for j in range(64): 
        data[i].update({clist[j]:yikye(clist[i],clist[j])})
        data[i].update({'---':([0,1,2],0,0)})

import pickle
with open('NdSd_all','wb') as f:
    pickle.dump(data,f)
        