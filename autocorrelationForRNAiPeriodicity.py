"""
Marissa Glover, November 05, 2018

Script to perform Spearman and Pearson autocorrelation with ribo-seq data following RNAi for Pule et al.

Inputs: reperiodicity/181104_SJA69MS_ZK617.1d.2_5414.txt
        reperiodicity/181104_SJA70MS_ZK617.1d.2_5414.txt
        reperiodicity/181104_SJA71MS_F11C3.3.2_4082.txt
        reperiodicity/181104_SJA72MS_F11C3.3.2_4082.txt
        reperiodicity/181104_SJA102MS_ZK617.1d.2_5414.txt
        reperiodicity/181104_SJA104MS_ZK617.1d.2_5414.txt
        reperiodicity/181104_SJA105MS_F07A5.7a.4_1147.txt
        reperiodicity/181104_SJA107MS_F07A5.7a.4_1147.txt
        reperiodicity/181104_SJA108MS_F11C3.3.2_2777.txt
        reperiodicity/181104_SJA110MS_F11C3.3.2_2777.txt
        reperiodicity/181104_SJA111MS_F11C3.3.2_556.txt
        reperiodicity/181104_SJA113MS_F11C3.3.2_556.txt

Output: plot of Spearman and Pearson autocorrelation coefficients or -log(p-values) vs. offsets

Need to create dictionary {position, number of reads}, create a list from dictionary values, perform Spearman and Pearson autocorrelation functions, then create plots.

Run as: python AutocorrelationForRNAiPeriodicity.py outPrefix reperiodicity/181104_SJA69MS_ZK617.1d.2_5414.txt reperiodicity/181104_SJA70MS_ZK617.1d.2_5414.txt reperiodicity/181104_SJA71MS_F11C3.3.2_4082.txt reperiodicity/181104_SJA72MS_F11C3.3.2_4082.txt reperiodicity/181104_SJA102MS_ZK617.1d.2_5414.txt reperiodicity/181104_SJA104MS_ZK617.1d.2_5414.txt reperiodicity/181104_SJA105MS_F07A5.7a.4_1147.txt reperiodicity/181104_SJA107MS_F07A5.7a.4_1147.txt reperiodicity/181104_SJA108MS_F11C3.3.2_2777.txt reperiodicity/181104_SJA110MS_F11C3.3.2_2777.txt reperiodicity/181104_SJA111MS_F11C3.3.2_556.txt reperiodicity/181104_SJA113MS_F11C3.3.2_556.txt
"""

import sys, math
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

def autoCorrSpearman(listaa):
    #Will compute Spearman's autocorrelation of x with itself. Will return a list of values, where each value is the Spearman autocorrelation with an offset of that index. Starting with 0 (self correlation, should be 1).
    k=len(listaa)
    a=[]
    for ii in range(0,k-1):
        jj=k-ii
        if len(list(set(listaa[:jj])))==1:
            break
        elif len(list(set(listaa[-jj:])))==1:
            break
        else:
            a.append([spearmanr(listaa[:jj],listaa[-jj:]),ii])
    return a

def autoCorrPearson(listaa):
    #Will compute Pearson's autocorrelation of x with itself. Will return a list of values, where each value is the Pearson autocorrelation with an offset of that index. Starting with 0 (self correlation, should be 1).
    k=len(listaa)
    b=[]
    for ii in range(0,k-1):
        jj=k-ii
        if len(list(set(listaa[:jj])))==1:
            break
        elif len(list(set(listaa[-jj:])))==1:
            break
        else:
            b.append([pearsonr(listaa[:jj],listaa[-jj:]),ii])
    return b

def makeDict(inputFile):
    #Will create dictionary of {position:number of reads} from input file. Can then use dictionary to create list for autocorrelations.
    c={}
    counter=0
    for line in inputFile:
        counter +=1
        if counter >1:
            key,value=line.strip().split()
            c[int(key.strip())]=float(value.strip())
    return c

def listFromDict(dictc):
    #Will make list from dictionary c to use for autocorrelation.
    d=[]
    for jj in dictc:
        d.append(dictc[jj])
    return d

"""
#mock data
keys=range(1,200)
values=[]
for key in keys:
    if key % 14 == 0:
        values.append(10)
    else:
        values.append(0)
c=dict(zip(keys,values))

#test mock data
bb=listFromDict(c)
cc=autoCorrSpearman(bb)
offsetS=[]
coefficientS=[]
pvalueS=[]
for i in range(len(cc)):
    entry=cc[i]
    offsetS.append(entry[1])
    coefficientS.append(entry[0][0])
    if entry[0][1]==0:
        pvalueS.append(-math.log(float(1e-60)))
    else:
        pvalueS.append(-math.log(entry[0][1]))
dd=autoCorrPearson(bb)
offsetP=[]
coefficientP=[]
pvalueP=[]
for i in range(len(cc)):
    entry=dd[i]
    offsetP.append(entry[1])
    coefficientP.append(entry[0][0])
    if entry[0][1]==0:
        pvalueP.append(-math.log(float(1e-60)))
    else:
        pvalueP.append(-math.log(entry[0][1]))
"""

#Will get files to call makeDict, listFromDict, autoCorrPearson, autoCorrSpearman on.
outPrefix=sys.argv[1]
myFiles=sys.argv[2:]
for myFile in myFiles:
    with open(myFile, 'r') as file:
        aa=makeDict(file)
    bb=listFromDict(aa)
    cc=autoCorrSpearman(bb)
    offsetS=[]
    coefficientS=[]
    pvalueS=[]
    for i in range(len(cc)):
        entry=cc[i]
        offsetS.append(entry[1])
        coefficientS.append(entry[0][0])
        if entry[0][1]==0:
            pvalueS.append(-math.log10(float(1e-60)))
        else:
            pvalueS.append(-math.log10(entry[0][1]))
    dd=autoCorrPearson(bb)
    offsetP=[]
    coefficientP=[]
    pvalueP=[]
    for i in range(len(dd)):
        entry=dd[i]
        offsetP.append(entry[1])
        coefficientP.append(entry[0][0])
        if entry[0][1]==0:
            pvalueP.append(-math.log10(float(1e-60)))
        else:
            pvalueP.append(-math.log10(entry[0][1]))
    
    #Will plot 4 line graphs with each file as own line; 2 for Spearman and 2 for Pearson, x-axis as offset, y-axis as coefficient or -log(pvalue).
    ax1=plt.subplot(221)
    ax1.plot(offsetS,coefficientS,label=file)
    ax1.set_ylabel('Correlation coefficient',fontsize=12)
    ax1.set_title('Spearman correlation',fontsize=14)
    ax1.set_xlim([0,50])
    ax1.set_ylim([-0.5,1])
    ax1.spines['left'].set_position('zero')
    ax1.spines['right'].set_position(('data',50))
    ax1.spines['top'].set_bounds(0,50)
    ax1.spines['bottom'].set_bounds(0,50)
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(9))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(3))
    ax2=plt.subplot(223)
    ax2.plot(offsetS,pvalueS,label=file)
    ax2.set_ylabel('-log(p-value)',fontsize=12)
    ax2.set_ylim([0,15])
    ax2.set_xlim([0,50])
    ax2.spines['left'].set_position('zero')
    ax2.spines['right'].set_position(('data',50))
    ax2.spines['top'].set_bounds(0,50)
    ax2.spines['bottom'].set_bounds(0,50)
    ax2.set_xlabel('Offset',fontsize=12)
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(9))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(3))
    plt.plot([0,200],[-math.log10(0.01),-math.log10(0.01)],'k--')
    plt.text(51,2,'pvalue<=\n0.01',fontsize=8)
    
    ax3=plt.subplot(222)
    ax3.plot(offsetP,coefficientP,label=file)
    ax3.set_ylabel('Correlation coefficient',fontsize=12)
    ax3.set_title('Pearson correlation',fontsize=14)
    ax3.set_xlim([0,50])
    ax3.set_ylim([-0.5,1])
    ax3.spines['left'].set_position('zero')
    ax3.spines['right'].set_position(('data',50))
    ax3.spines['top'].set_bounds(0,50)
    ax3.spines['bottom'].set_bounds(0,50)
    ax3.xaxis.set_major_locator(ticker.MultipleLocator(9))
    ax3.xaxis.set_minor_locator(ticker.MultipleLocator(3))
    ax4=plt.subplot(224)
    ax4.plot(offsetP,pvalueP,label=file)
    ax4.set_ylabel('-log(p-value)',fontsize=12)
    ax4.set_ylim([0,15])
    ax4.set_xlim([0,50])
    ax4.set_xlabel('Offset',fontsize=12)
    ax4.spines['left'].set_position('zero')
    ax4.spines['right'].set_position(('data',50))
    ax4.spines['top'].set_bounds(0,50)
    ax4.spines['bottom'].set_bounds(0,50)
    ax4.xaxis.set_major_locator(ticker.MultipleLocator(9))
    ax4.xaxis.set_minor_locator(ticker.MultipleLocator(3))
    plt.plot([0,200],[-math.log10(0.01),-math.log10(0.01)],'k--')
    plt.text(51,2,'pvalue<=\n0.01',fontsize=8)
    plt.legend(loc='upper center',bbox_to_anchor=(0.5,-0.3),ncol=2,fontsize=10)
    plt.tight_layout(w_pad=5)
plt.show()
plt.savefig('test.svg',bbox_inches='tight')
plt.close('all')
