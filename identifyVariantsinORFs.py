"""
Marissa Glover, January 9, 2019

Script to idenfity variants that are in ORFs

inputs: filtered.vcf filtered variant file from filterHomozygousVariantsFromVCF.py
directory and file prefix /data9/marissa/genomes/171218_historicalGenomeT2A/170612_genomeWithUnc-54ChrAsItAppearsInPD4092 which has annotated text files for each chromosome

output: text file of list of SNPs in ORFs with columns for chromosome, position, and gene name

Run as python identifyVariantsinORFs.py .filtered.vcf_file directory_and_file_prefix outPrefix
"""

import sys, linecache
from logJosh import Tee


def makeListOfVariants(vcfFile):
    #Given vcf File, will make a list of tuples where each tuple is chromosome,position
    #initialize list
    listOfVariants=[]
    #open input file
    with open(vcfFile,'r') as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                line=line.split("\t")
                listOfVariants.append((line[0],line[1]))
    return listOfVariants

def findVariantsInORFs(listOfVariants,annotFilePrefix):
    #Given the listOfVariants as output from makeListOfVariants and the annotated chromosome files, will find variants that are in ORFs
    allVariantsInGenes=[]
    for variant in listOfVariants:
        chromosome=variant[0]
        position=int(variant[1])
        chromosomeFile=annotFilePrefix+"."+chromosome+".txt"
        #Go find line in chromosome file that includes the variant
        line= linecache.getline(chromosomeFile,position)
        if "WBGene" in line:
            variantsInGenes=[]
            splitLine=line.split("\t")
            for entry in splitLine:
                if entry.startswith("WBGene"):
                    splitEntry=entry.split(":")
                    variantsInGenes.append(splitEntry[0])
            allVariantsInGenes.append([chromosome,str(position)]+variantsInGenes)
    return allVariantsInGenes

def makeOutputFile(variantsInORFs,outPrefix):
    #Given variantsInORFs, will create text file for output
    with open(outPrefix,'w') as f:
        for line in variantsInORFs:
            f.write("\t".join(line))
            f.write("\n")

def main(args):
    vcfFile,annotFilePrefix,outPrefix=args[0:]
    #next make list from .vcf file of chromosome, position for all variants
    listOfVariants=makeListOfVariants(vcfFile)
    #next check if variant position falls in an ORF
    variantsInORFs=findVariantsInORFs(listOfVariants,annotFilePrefix)
    #next output file
    makeOutputFile(variantsInORFs,outPrefix)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
