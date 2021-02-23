# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 12:03:48 2018

@author: khaik
"""
import numpy as np
import pickle

class dNdS:
    def __init__(self, rseq, oseq, dirn='v',
                 window=3, reload=True):
        self.dirn = dirn
        self.rseq = rseq
        self.oseq = oseq
        self.oSize = len(oseq)
        self.window = window

        self.reload = reload

        nseq = len(self.rseq)
        self.nAA = int(np.floor(nseq/3))

        if self.dirn=='h':
            self.n = 0
            self.s = 0
            self.nd = 0
            self.sd = 0
            self.dn = 0
            self.ds = 0
            self.dnds = 0
            self.mpos = []

        if self.dirn=='v':
            self.seqSize = len(oseq)
            self.n = np.zeros((1,self.nAA),float)
            self.s = np.zeros((1,self.nAA),float)
            self.zeroIndex = [[] for _ in range(5)]     #n,s,pn,ps,dnds
            self.nd = np.zeros((self.oSize,self.nAA),float)
            self.sd = np.zeros((self.oSize,self.nAA),float)
            self.dn = np.zeros((1,self.nAA),float)
            self.ds = np.zeros((1,self.nAA),float)
            self.dnds = np.zeros((1,self.nAA),float)
            self.mpos = np.zeros((1,self.nAA),float)

    def saveData(self):
        data = {'n':self.n,'s':self.s,
                'dnds':self.dnds,'mpos':self.mpos}
        self.data = data

    def setup(self):
        with open('dictionaries/cdict','rb') as f:
            self.cDict = pickle.load(f)

        with open('dictionaries/pdict','rb') as f:
            self.pDict = pickle.load(f)

        with open('dictionaries/NdSd','rb') as f:
            self.cdDict = pickle.load(f)

    #STEP1: n,s from ref; need to get cdict from somewhere
    def getNS_h(self):
        for i in range(self.nAA):
            self.n += self.cDict.get(self.rseq[i*3:i*3+3])[1]
            self.s += self.cDict.get(self.rseq[i*3:i*3+3])[2]

        self.n/=3
        self.s/=3

    def getNS_v(self):
        for i in range(self.nAA):
            self.n[:,i] = self.cDict.get(self.rseq[i*3:i*3+3])[1]
            self.s[:,i] = self.cDict.get(self.rseq[i*3:i*3+3])[2]
            #self.d[:,i] = self.dDict.get(self.rseq[i*3:i*3+3])
    #STEP2: nd,sd (save position!) ok (csdict aka data)
    def getNdSd_h(self):
        for i in range(self.nAA):
            rcodon = self.rseq[i*3:i*3+3]
            ocodon = self.oseq[i*3:i*3+3]
            if rcodon == ocodon:
                continue
            else:
                dLoc,Nd,Sd = self.codonNdSd(i,rcodon,ocodon)

                self.mpos.append(dLoc)
                self.nd+=Nd
                self.sd+=Sd

    def getNdSd_v(self):
        for i in range(self.nAA):
            rcodon = self.rseq[i*3:i*3+3]
            for j in range(self.oSize):
                seq = self.oseq[j]
                ocodon = seq[i*3:i*3+3]
                if rcodon == ocodon:
                    continue
                else:
                    dLoc,Nd,Sd = self.codonNdSd(i,j,rcodon,ocodon)
                    self.nd[j,i]+=Nd
                    self.sd[j,i]+=Sd
                    self.mpos[0,i]+=1

    #save postiion, stop codon?
    def codonNdSd(self,i,j,c1,c2):
        dLoc = None
        Nd = 0
        Sd = 0

        for codonlist in self.cdDict:
            if codonlist['codon']==c1:
                mycodon = codonlist.get(c2,None)
                if mycodon==None:
                    print('{}{}{}{}'.format('hmm.. ',i,' ',j),end=',')
                    mycodon = ([],0,0)
                Nd,Sd = mycodon[1], mycodon[2]
                dLoc = [i*3+x for x in mycodon[0]]
        return dLoc,Nd,Sd

    def calcF(self,n,s,nd,sd):
        pn = nd/n if n!=0 else np.nan
        ps = sd/s if s!=0 else np.nan
        dn = -3/4*(np.log(1-4*pn/3))
        ds = -3/4*(np.log(1-4*ps/3))
        dnds = dn/ds if ds!=0 else np.nan
        return dn,ds,dnds

    def main(self):
        if self.reload:
            self.setup()

        if self.dirn=='h':
            self.getNS_h()
            self.getNdSd_h()
            self.dn,self.ds,self.dnds = self.calcF(self.n,self.s,self.nd,self.sd)

        elif self.dirn=='v':
            self.getNS_v()
            self.getNdSd_v()

            dn = np.zeros((self.window,self.nAA),float)
            ds = np.zeros((self.window,self.nAA),float)
            dnds = np.zeros((self.window,self.nAA),float)
            for z in range(self.window):
                dn[z],ds[z],dnds[z] = self.slidingWindow(self.window,start=z)

            for i in range(self.nAA):
                self.mpos[0,i]/=self.oSize
                self.dn[0,i] = np.nanmean(dn[:,i])
                self.ds[0,i] = np.nanmean(ds[:,i])
                self.dnds[0,i] = np.nanmean(dnds[:,i])

        self.saveData()

    def slidingWindow(self, window, start=0):
        nBB = int(np.ceil(self.nAA/window) + 1)
        nCC = self.nAA/window

        dn = np.empty((1,self.nAA),float)
        dn[:] = np.nan
        ds = np.empty((1,self.nAA),float)
        ds[:] = np.nan
        dnds = np.empty((1,self.nAA),float)
        dnds[:] = np.nan


        for i in range(nBB):
            ws = 0 if i==0 else None if (i==nBB-1) and (start>=nCC) else (i-1)*window+start
            we = start if i==0 else None if ws+window>self.nAA else ws+window

            n = np.nanmean(self.n[0,ws:we])
            s = np.nanmean(self.s[0,ws:we])
            nd = np.nanmean(self.nd[:,ws:we])
            sd = np.nanmean(self.sd[:,ws:we])
            a,b,c = self.calcF(n,s,nd,sd)

            w=0
            while (w<window) and (ws+w<self.nAA):
                dn[0,ws+w]=a
                ds[0,ws+w]=b
                dnds[0,ws+w]=c
                w+=1

        return dn, ds, dnds

    def manyWindows(self, windList=[3,5,10]):
        windDict = {}
        for windNum in windList:
            windArray = np.empty((windNum,self.nAA),float)
            for z in range(windNum):
                windArray[z] = self.slidingWindow(windNum, start=z)[2]

            dnds = np.empty((1,self.nAA),float)
            for i in range(self.nAA):
                dnds[0,i] = np.nanmean(windArray[:,i])

            windDict.update({str(windNum):dnds})
        return windDict

    def slidingWindow2(self, window, start=0):
        nBB = int(np.ceil(self.nAA/window) + 1)
        nCC = self.nAA/window

        ndsd = np.empty((1,self.nAA),float)
        ndsd[:] = np.nan
        ns = np.empty((1,self.nAA),float)
        ns[:] = np.nan

        for i in range(nBB):
            ws = 0 if i==0 else None if (i==nBB-1) and (start>=nCC) else (i-1)*window+start
            we = start if i==0 else None if ws+window>self.nAA else ws+window

            n = np.nanmean(self.n[:,ws:we])
            s = np.nanmean(self.s[:,ws:we])
            nd = np.nanmean(self.nd[:,ws:we])
            sd = np.nanmean(self.sd[:,ws:we])
            a = nd/sd
            b = n/s

            w=0
            while (w<window) and (ws+w<self.nAA):
                ndsd[0,ws+w]=a
                ns[0,ws+w]=b
                w+=1

        return ndsd,ns

    def manyWindows2(self, windList=[3,5,10]):
        windDict = {}
        for windNum in windList:
            windArray = np.empty((windNum,self.nAA),float)
            for z in range(windNum):
                windArray[z] = self.slidingWindow(windNum, start=z)[2]

            dnds = np.empty((1,self.nAA),float)
            for i in range(self.nAA):
                dnds[0,i] = np.nanmean(windArray[:,i])

            windDict.update({str(windNum):dnds})
        return windDict

#
#ID_step3_sw = dNdS(pa01,sequence1_xAE,dirn='v')
#ID_step3_sw.main()
#dnds3_all[0] = ID_step3_sw.data['dnds']['3'][2]
#dnds3_all[1] = ID_step3_sw.slidingWindow(3,start=1)[2]
#dnds3_all[2] = ID_step3_sw.slidingWindow(3,start=2)[2]

#for i in range(419):
#    dnds3[0,i]=np.nanmean(dnds3_all[:,i])

'''random functions'''
    def slidingWindow3(slide, window, nAA, start=0):
        nBB = int(np.ceil(nAA/window) + 1)
        nCC = nAA/window

        b = np.empty((1,nAA),float)

        for i in range(nBB):
            ws = 0 if i==0 else None if (i==nBB-1) and (start>=nCC) else (i-1)*window+start
            we = start if i==0 else None if ws+window>nAA else ws+window

            a = np.nanmean(slide[ws:we])

            w=0
            while (w<window) and (ws+w<nAA):
                b[0,ws+w]=a
                w+=1

        return b

    def zh(slide, window, nAA):
        zh = np.zeros((window,nAA))
        for z in range(window):
            zh[z]=slidingWindow3(slide,window,nAA,start=z)
        c = np.zeros((1,nAA))
        for i in range(nAA):
            c[0,i] = np.nanmean(zh[:,i])
        return c

    def savePDBList(array,indexList,dataname):
        a = array[:-1].tolist()
        k = 0
        while k<len(indexList):
            a.pop(indexList[k]-k)
            k+=1
        filename='pdb_dnds/'+dataname
        with open(filename,'wb') as f:
            pickle.dump(a,f,protocol=2)
        print(max(a))

    def main(windArray,window,nAA,step3):
        a = np.zeros((4,nAA))
        a[0] = windArray
        for i in range(nAA):
            a[1,i]=dDict.get(step3.rseq[i*3:i*3+3])
        a[3]=zh(a[1],window,nAA)
        for i in range(nAA):
            a[2,i]=a[0,i]/a[3,i]
        return a
