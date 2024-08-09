"""
Marissa Glover, January 11, 2019

Script to identify genes that are hit in more than 1 sample

Input: variantsInORFs files from identifyVariantsinORFs.py

Output: list of genes that are hit in more than 1 sample

Run as python genesHitInMultipleSamples.py outPrefix variantsInORFs_files
"""

import sys
from logJosh import Tee

def parseAndMakeDict(variantsFile):
    #Given the variantsInORFs file, will parse the file and create a dictionary for each file with gene names as keys and number of hits as values
    genes={}
    #open input files
    for fileName in variantsFile:
        with open(fileName,'r') as f:
            for line in f:
                line=line.strip().split('\t')
                geneNames=line[2:]
                for geneName in geneNames:
                    if geneName in genes:
                        genes[geneName]+=1
                    else:
                        genes[geneName]=1
    return genes

def findGenes(dictOfGenes):
    genes=[]
    for entry in dictOfGenes:
        if dictOfGenes[entry]>1:
            #genes.append(entry)
            if dictOfGenes[entry]<7:
                print entry,'\t',dictOfGenes[entry]
        else:
            pass
    print genes

def main(args):
    outPrefix=args[0]
    variantsFile=args[1:]
    #next parse variantsInORFs file and make dictionary for each file with gene names as keys and number of hits as values
    dictOfGenes=parseAndMakeDict(variantsFile)
    genesHitMoreThanOnce=findGenes(dictOfGenes)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
