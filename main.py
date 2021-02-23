# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 00:05:18 2018

@author: khaik
"""
import numpy as np
from matplotlib import pyplot as plt
import pickle

'''
CLASS
1. getSeq: output(gene_list) from accession number, gene of interest
2. alignSeq: output(alignedNT) from gene_list
3. dNdS: output(dnds) from alignedNT
'''
#read fasta, output strings
from Bio import SeqIO
filename = 'icd_test.fasta'
seq_data = []
for seq_record in SeqIO.parse(filename,'fasta'):
    seq_data.append(seq_record)

test1 = seq_data[0].seq._data
test2 = seq_data[1].seq._data

#parsing online GenBank records
from Bio import Entrez
from Bio import SeqIO
Entrez.email= 'kks36@cam.ac.uk'
with Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id='ICD[sym]') as handle:
    seq_record = SeqIO.read (handle,'fasta')

def blastsearch(filename, sequence, txid, size=50):
    result_handle = NCBIWWW.qblast('blastn','nt',sequence,
                                   entrez_query='txid'+txid+'[organism]',
                                   hitlist_size=size)

    blastname = filename + '.xml'
    with open(blastname,'w') as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()


#<3
#hsps[0] is always (usually?) the most accurate alignment, but can consider expected < 0.04
#what to do when query_start != 0?
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez


'''
-old blastsearch that worked-
def blastsearch(filename, sequence, txid, size=50):
    result_handle = NCBIWWW.qblast('blastn','nt',sequence,
                                   entrez_query='txid'+txid+'[organism]',
                                   hitlist_size=size)

    blastname = filename + '.xml'
    with open(blastname,'w') as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()
    blastparse = NCBIXML.parse(open(blastname))
    blastlist = list(blastparse)

    return blastlist
'''
class getSeq_blast:
    def __init__(self, filename, geneLength=None):
        self.gene_list = []
        self.filename = filename
        self.geneLength = geneLength

    def setup(self):
        Entrez.email= 'kks36@cam.ac.uk'
        self.blast_records = NCBIXML.read(open(self.filename))
        for alignment in blast_records.alignments:
            orgn = alignment.hit_def
            accs = alignment.accession

            if alignment.hsps[0].align_length==self.geneLength:
                seq = alignment.hsps[0].sbjct

            else:

                geneStart = alignment.hsps[0].sbjct_start
                geneEnd = alignment.hsps[0].sbjct_end
                self.blastfetch()

            data = {'organism':orgn,
                    'accession':accs,
                    'sequence':seq,
                    'length':self.geneLength}
            self.gene_list.append(data)



            data = {'organism':self.currOrgn,
                    'accession':self.currID,
                    'sequence':new_seq,
                    'length':len_seq,
                    'gene':self.currGene}
            self.gene_list.append(data)


            my_dict = {'orgn':alignment.hit_def,
                       'length':alignment.hsps[0].align_length,
                       'query':alignment.hsps[0].query,
                       'match':alignment.hsps[0].match,
                       'sbjct':alignment.hsps[0].sbjct,
                       'ident':alignment.hsps[0].identities/alignment.hsps[0].align_length,
                       'rangeS':[alignment.hsps[0].sbjct_start,alignment.hsps[0].sbjct_end],
                       'rangeQ':[alignment.hsps[0].query_start,alignment.hsps[0].query_end],
                       'score':alignment.hsps[0].score
                       }

            hits.append(my_dict)
        return hits

    def blastfetch(seq_id, n, hsp):
        geneStart = hsp.sbjct_start
        geneEnd = hsp.sbjct_end

        if hsp.query_start != 1:
            geneStart += (1-hsp.query_start)

        if hsp.query_end != self.geneLength:
            geneEnd += (1-geneStart + self.geneLength - 1

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
            self.gene_list.append(data)



from Bio import SeqIO
filename = 'icd_test.fasta'
seq_data = []
for seq_record in SeqIO.parse(filename,'fasta'):
    seq_data.append(seq_record)

test1 = seq_data[0].seq._data
test2 = seq_data[1].seq._data



def getsummary(dictionary, key):
    my_list = []
    for strain in dictionary:
        my_list.append(strain.get(key))
    return my_list

#only the impt bits
pa10 = blastsearch('pseudomonas',test1,'287',size=50) #aka pa10_hits
titlelist = getsummary(pa10,'title')

'''
get sequence from summary - ok
- for ncbi_summary_nogeneLoc: save as Accession list

OUTPUT: gene_list (including wrong genes)

SOURCE ID
0 - ncbi_summary
1 - ncbi_summary_noGeneLoc

Things to catch
1. length
2. start codon?
3. ident

- must include roughGene if name is potentially long?
'''
from Bio import SeqIO, Entrez
import re
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



'''general code'''
#view alignment
def viewAlignment(seq,start,end):
    print(seq.query[start:end])
    print(seq.match[start:end])
    print(seq.sbjct[start:end])

viewAlignment(bp,0,100)

#match codon to aa
def matchAA(codon):
    aa = None
    for key,value in pdict.items():
        if codon in value:
            aa = key
    return aa

#translate nt seq to protein seq
def translate(nseq):
    nAA = int(np.floor(len(nseq)/3))
    pseq = ''
    for i in range(nAA):
        pseq += matchAA(nseq[i*3:i*3+3])
    return pseq
test1_p = translate(test1)

gene_dict = readGR('gene_result')
with open('icd_summary52','wb') as f:
    pickle.dump(gene_dict,f)


#plotting
class plotStuff:

    def four_plots(data):
        plt.figure(1)
        plt.subplot(411)
        plt.plot(data['dn'][0])
        plt.title('dN')

        plt.subplot(412)
        plt.plot(data['ds'][0])
        plt.title('dS')

        plt.subplot(413)
        plt.plot(data['dnds'][0])
        plt.title('dNdS')

        plt.subplot(414)
        plt.plot(data['mpos'][0])
        plt.title('fraction of codon diff')

        plt.tight_layout()
        plt.show()

    def plot_dnds(data,dots):
        s = data['s'][0]
        s0 = np.zeros(len(s))
        for i in range(len(s)):
            if s[i]==0:
                s0[i]=-0.5
            else:
                s0[i]=np.nan

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(data['dnds'][0],'ro')
        ax.plot(dots,'bo')
        ax.plot(s0,'o',mfc='none',mec='r',mew=1)
        ax.set_xlabel('codon posn')
        ax.set_ylabel('dNdS')

    def plot_SW(data):
        for key,value in data.items():
            plt.plot(value[2][0],label=str(key))
        plt.set_xlabel('codon posn')
        plt.set_ylabel('dNdS')

    def plot_circle():
        icd_circle = [113, 115, 119, 129, 153, 160, 230,
                      283, 307,
                      306,
                      153,
                      160,
                      292, 345, 391, 344]
        idh_circle = [132, 135, 139, 145, 547, 420, 255,
                      350, 548,
                      547,
                      411,
                      182,
                      358, 589, 649, 588]
        for i in range
        circle = plt.Circle((),size,color=color)

    def plot_ICDvIDH():
        from cycler import cycler
        size = len(alignments)
        a = plt.get_cmap('gist_rainbow')
        fig, ax = plt.subplots()
        ax.set_prop_cycle(cycler('color',[a(1.*i/size) for i in range(size)]))

        x = [i for i in range(742)]
        y = windArray_idh['10'][0]
        ax.plot(y,c='#F4C2C2')

        for i in range(len(alignments)):
            x0 = alignments[i][0][0]
            x1 = alignments[i][0][1]
            ax.plot(x[x0:x1],y[x0:x1])

        x = [i for i in range(419)]
        y = windArray_icd['10'][0]
        ax.plot(y,c='#F4C2C2')

        for i in range(len(alignments)):
            x0 = alignments[i][1][0]
            x1 = alignments[i][1][1]
            ax.plot(x[x0:x1],y[x0:x1])

        plt.tight_layout

'''
dNdS - ok
0. align (from BLAST)
1. reference seqeunce: total n, total s (rmb to /3) #cDict, getNS
2. compare ref and obs: total nd, total sd          #cdDict
   - if same codon, nd sd = 0
   - if diff by 1 nt: either nd+=1 or sd
   - if diff by >1 nt: trace path and average. total should = nt diff
3. pn = nd/n; ps = sd/s
4. dn = -3/4ln(1-4pn/3); ds = -3/4ln(1-4ps/3)
5. dnds = dn/ds
6. coverage ..?

Qns
- what happens if s=0? return 0 or nan

variables
- seq: h-list(2); v-list

'''

class dNdS:
    def __init__(self,rseq,oseq,dirn='h',
                 window=1,start=0,reload=True):
        self.dirn = dirn
        self.rseq = rseq
        self.oseq = oseq
        self.oSize = len(oseq)
        self.window = window
        self.start = start

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
            self.mpos = np.zeros((1,self.nAA),float)

    def saveData(self):
        data = {'n':self.n,'s':self.s,
                'dn':self.dn,'ds':self.ds,
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
            self.n[:,i] = self.cDict.get(self.rseq[i*3:i*3+3])[1])
            self.s[:,i] = self.cDict.get(self.rseq[i*3:i*3+3])[2])

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

    def vsExpected(self,):



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
            self.dn,self.ds,self.dnds = self.slidingWindow(self.window, self.start)

            for i in range(self.nAA):
                self.mpos[0,i]/=self.oSize

        self.saveData()

    def slidingWindow(self, window, start=0):
        nBB = int(np.ceil(self.nAA/window))

        dn = np.zeros((1,nBB),float)
        ds = np.zeros((1,nBB),float)
        dnds = np.zeros((1,nBB),float)

        for i in range(nBB):
            ws = i+start
            we = k+window if i!=nBB-1 else None

            n = np.nanmean(self.n[0,ws:we])
            s = np.nanmean(self.n[0,:ws:we])
            nd = np.nanmean(self.nd[:,ws:we])
            sd = np.nanmean(self.sd[:,ws:we])
            dn[0,i],ds[0,i],dnds[0,i] = self.calcF(n,s,nd,sd)

        return dn, ds, dnds

#gene_dict
jfr = dNdS(pa01,koba,dirn='v')
jfr.main()
jfr_data0 = jfr.data

if __name__ == "__main__":
    refSeq = test1   #PA01
    obsSeq = test2
    data = []

    #h
    for i in len(gene_dict):
        obsSeq = gene_dict[i][1]
        newGene = dNdS(refSeq,obsSeq,dirn='h')
        newGene.main()
        data.append(newGene)

    #v
    newGene = dNdS(refSeq,seqList,dirn='v')
    newGene.main()


'''align genes
- input nt list (includes rubbish nt, and wrong nt)
- output aligned nt list (should i save the rubbish nt somewhere?
- aligned protein seq is not saved hmmm..
- DO NOT CHANGE REFERENCE SEQ
0. filter
1. to aa seq
2. match aa
3. dNdS of nucleotide seq - how to handle gaps?

OUTPUT: alignedNT
'''
refSeq = test1  #str
obsSeq = list1  #list

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

class alignSeq:
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

    def okr(self,i,oldseq):
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

    def main(self):
        self.setup()

        rlen = len(self.rnseq)
        for i in range(len(self.onseq)):
            seq = self.onseq[i]

            for nt in seq:
                if nt not in ['T','C','A','G']:
                    continue

            if len(seq)%3 != 0:
                continue

            if (len(seq)==rlen):
                self.alignedNT.append({'index':i,'sequence':[self.onseq[i]],
                                       'insertion':[],'score':[]})
            else:
                self.okr(i,seq)

            if i%100==0:
                print(i,end=' ')


'''compare with ecoli'''
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from Bio import pairwise2

refFile = 'icd_RefPA01'
pdbFile = 'protein structure/ACEK-ICDH.pdb'
outputFile = 'pdb_dnds/icd_ecoli'

class zh:
    def __init__(self, refFile, pdbFile):
        with open(refFile, 'rb') as f:
            self.refN = pickle.load(f)

        recordDict = {}
        with open(pdbFile, "rU") as handle:
            for record in SeqIO.parse(handle, "pdb-seqres"):
                recordDict.update({record.id:record.seq})
        self.orgnP = recordDict['4P69:C']._data

        with open('dictionaries/pdict','rb') as f:
            self.pdict = pickle.load(f)

        self.refP = self.translate(self.refN)[:-1]

        with open('icd_dnds','rb') as f:
            self.refD = pickle.load(f)

        self.orgnD = []

    def matchAA(self,codon):
        aa = None
        for key,value in self.pdict.items():
            if codon in value:
                aa = key
        return aa

    def translate(self,nseq):
        nAA = int(np.floor(len(nseq)/3))
        pseq = ''
        for i in range(nAA):
            pseq += self.matchAA(nseq[i*3:i*3+3])
        return pseq

    def an_inspiration(self):
        alignments = pairwise2.align.globalds(self.refP,self.orgnP,blosum62,-10,-0.5)
        aligned = alignments[0][1]
        indexList = []
        for i in range(len(aligned)):
            if aligned[i]=='-':
                indexList.append(i)

        k = 0
        while k<len(indexList):
            self.refD.pop(indexList[k]-k)
            k+=1

    def save_data(self):
        with open(outputFile,'wb') as f:
            pickle.dump(self.refD, outputFile, protocol=2)


'''display sequence alignment'''
from matplotlib import cm, colors
import seaborn as sns
import pandas as pd

class displayHM:
    def __init__(self,nt,pt,indexList,
                 ndict='dictionaries/n2n',pdict='dictionaries/p2n_char'):
        self.nt = nt
        self.pt = pt
        self.indexList = indexList

        with open(ndict,'rb') as f:
            self.ndict = pickle.load(f)
        with open(pdict,'rb') as f:
            self.pdict = pickle.load(f)

        self.ncList=ndict['cList']
        self.pcList=pdict['cList']

    def setup():
        self.Nmat = self.to_df(self.nt,self.ndict)
        self.Pmat = self.to_df(self.pt,self.pdict)

    def to_df(self,sequence,dictionary):
        sSize=len(sequence)
        nSize=len(sequence[0])
        N_mat = np.zeros((sSize, nSize),int)
        for i in range(sSize):
            for j in range(nSize):
                N_mat[i,j] = dictionary.get(sequence[i][j],0)
        df = pd.DataFrame(N_mat)
        return df

    def plotSNS():
        fig,ax = plt.subplots()
        hax = sns.heatmap(self.Nmat,cmap=self.ncList)
        ax.set_xlabel('codon posn')
        ax.set_ylabel('organism')
        ax.set_title('nucleotide')
        plt.tight_layout()

        fig,ax = plt.subplots()
        hax = sns.heatmap(self.Pmat,cmap=self.pcList)
        ax.set_xlabel('amino acid posn')
        ax.set_ylabel('organism')
        ax.set_title('protein')
        plt.tight_layout()


    def main(self):
        self.setup()
        self.plotSNS()

    #acidic (pink), basic (blue), polar (orange), non-polar (green)
    pList0 = ['D','E',
              'H','K','R',
              'S','T','N','Q',
              'A','V','L','I','M','F','Y','W','P','G','C',
              '*']

    sns.color_palette('Greens',20)

'''fasta'''
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio import SeqIO

def to_fasta(seqList, indexList, filename):
    fastafile = []
    for i in range(len(seqList)):
        record = SeqRecord(Seq(seqList[i],generic_nucleotide),
                           id=indexList[i],
                           description='concept')
        fastafile.append(record)
    SeqIO.write(fastafile,filename,'fasta')

'''tree''''
from Bio import Phylo

class drawTree:
    def __init__(self, dndfile, df, orgnList):
        self.tree = Phylo.read(dndfile,'newick')
        self.orgnList = orgnList
        self.df = df

    def sortDF(self, df):
        df = df.assign(species=pd.Series(self.orgnList).values)
        df = df.set_index('species')
        df = df.assign(temp=pd.Series([i for i in range(len(df))]).values)
        df['temp'] = pd.Categorical(df['temp'],self.tree_values)
        df = df.sort_values('temp')
        df = df.drop(columns=['temp'])
        return df

    def setup(self):
        tree_dict = self.tree.get_terminals()
        self.tree_values = []
        for i in range(len(tree_dict)):
            self.tree_values.append(int(tree_dict[i].name))

    def main(self):
        self.setup()
        self.sortDF(self.df)

#from cmd
clustalw2 -infile=C:\Users\khaik\Documents\camb\dNdS\concept.faa

#taxid dict (id, phylum,class,genus)
#phylum key: 0 proteobacteria, 1 chlamydias, 2 spirochetes, 3 cyanobacteria, 4 gram-positive bacteria
#class key: 0 alpha, 1 beta, 2 gamma, 3 delta, 4 epsilon, 5 not proteobacteria
taxid = {'P.aeruginosa':    ['287',   0,2,'Pseudomonas'],
         'P.chlororaphis':  ['587753',0,2,'Pseudomonas'],
         'P.fluorescens':   ['294',   0,2,'Pseudomonas'],
         'P.pertucinogena': ['86175', 0,2,'Pseudomonas'],
         'P.putida':        ['303',   0,2,'Pseudomonas'],
         'P.stutzeri':      ['316',   0,2,'Pseudomonas'],
         'P.syringae':      ['317',   0,2,'Pseudomonas'],
         'E.coli':          ['562',   0,2,'Escherichia'],
         'M.abscessus':     ['36809', 4,5,'Mycobacteroides'],
         'S.sonnei':        ['624',   0,2,'Shigella'],
         'A.baumannii':     ['470',   0,2,'Acinetobacter'],
         'N.meningitidis':  ['487',   0,1,'Neisseria'],
         'B.pseudomallei':  ['28450', 0,1,'Burkholderia'],
         'K.pneumoniae':    ['573',   0,2,'Klebsiella'],
         'S.suis':          ['1307',  4,5,'Streptococcus'],
         'E.cloacae':       ['550',   0,2,'Enterobacter'],
         'L.interrogans':   ['173',   2,5,'Leptospira'],
         'M.xanthus':       ['34',    0,3,'Myxococcus'],
         'S.marcescens':    ['615',   0,2,'Serratia'],
         'S.aureus':        ['1280',  4,5,'Staphylococcus'],
         'L.monocytogenes': ['1639',  4,5,'Listeria'],
         'V.cholerae':      ['666',   0,2,'Vibrio'],
         'B.bronchiseptica':['518',   0,1,'Bordetella'],
         'H.pylori':        ['210',   0,4,'Helicobacter'],
         'S.enterica':      ['28901', 0,2,'Salmonella'],
         'E.albertii':      ['208962',0,2,'Escherichia'],
         'Y.pestis':        ['632',   0,2,'Yersinia']
        }
