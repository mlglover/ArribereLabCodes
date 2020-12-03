"""
Marissa Glover, April 6, 2020

Script to analyze and compare sequences.

Input: fasta file of genome sequences
Output: nucleotide composition, number ORFs, length of ORFs
"""

import sys
import sequenceAnalysis as seqA
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

class CommandLine():
    def __init__(self, inOpts=None):
        import argparse
        self.parser = argparse.ArgumentParser(description='Analyze nucleotide sequences')
        self.parser.add_argument('inFile', action='store', nargs='+', help='Input is fasta file(s) of coding sequences')
        self.parser.add_argument('outFile', action='store', help='Output file name for results')
        self.parser.add_argument('-s', '--starts', action='append', default=['ATG'], nargs='?', help='Start codons used')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

################################################################################

class DrawCharts:
    def codonHeatMap(codonFreqs, speciesList, nucFreqs, AAFreqs):
        fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, constrained_layout=True)
        #fig, (ax1, ax2, ax3) = plt.subplots(ncols=3)
        
        #heatmap of codon frequencies
        #make dataframe from codon count dictionaries
        codonData = pd.DataFrame(codonFreqs).set_index(species for species in speciesList).transpose()
        #make series from rnaCodonTable to get amino acids associated with codons
        aaData = pd.Series(seqA.NucParams.rnaCodonTable)
        #make multiindexes dataframe indexed by codon and amino acid
        multiFrame = codonData.set_index([[value for value in aaData.values], [index for index in aaData.index]]).sort_index()
        print(multiFrame)
        
        codonHeatmap = sns.heatmap(multiFrame, linewidths=.5, xticklabels=True,
                        yticklabels=True, cbar_kws={'label':'RSCU values', 'shrink':0.7}, ax=ax1,
                        cmap='PiYG', center=1.0)
        #need to figure out how to label ticks with 2 labels (codons and AAs)
        #codonHeatmap.set_yticklabels('%s \n %s'%((multiFrame.index.get_level_values(0)), (multiFrame.index.get_level_values(1))))
        codonHeatmap.tick_params(axis='y', labelsize=8)
        codonHeatmap.tick_params(axis='x', labelsize=12, rotation=90)
        ax1.set_xlabel('Species', fontsize=16)
        ax1.set_ylabel('Codon', fontsize=16)
        ax1.set_title('Relative Synonymous \n Codon Usage', fontsize=18)
        
        #barplot for nucleotide frequencies
        #make dataframe from nucleotide count dictionaries
        nucFreqs2 = pd.DataFrame(nucFreqs, index=speciesList)
        #melt dataframe to get columns of species, nueclotide, and values
        nucFreqs2 = nucFreqs2.reset_index().melt(id_vars='index')
        print(nucFreqs2)
        makeBarPlot = sns.barplot(x='variable', y='value', hue='index', data=nucFreqs2,
                        ax=ax2, palette='Set3', edgecolor=(0,0,0))
        
        ax2.set_xlabel('Nucleotide', fontsize=16)
        ax2.set_ylabel('Frequency', fontsize=16)
        ax2.set_title('Nucleotide Frequencies', fontsize=18)
        
        #heatmap of amino acid frequencies
        AAData = pd.DataFrame(AAFreqs).set_index(species for species in speciesList).transpose().sort_index()
        print(AAData)
        AAHeatmap = sns.heatmap(AAData, linewidths=.5, xticklabels=True,
                        yticklabels=True, cbar_kws={'label':'%', 'shrink':0.7}, ax=ax3)
        AAHeatmap.tick_params(axis='x', labelsize=12, rotation=90)
        AAHeatmap.tick_params(axis='y', labelsize=12, rotation=0)
        ax3.set_xlabel('Species', fontsize=16)
        ax3.set_ylabel('Amino Acid', fontsize=16)
        ax3.set_title('Amino Acid Frequencies', fontsize=18)
        
        #fig.tight_layout()
        plt.savefig(CommandLine().parser.parse_args().outFile, dpi=300)
        plt.show()
################################################################################

def main(inCL=None):
    if inCL is None:
        cmdLine = CommandLine()
    else:
        cmdLine = CommandLine(inCL)
    
    codonDicts = []
    codonFreqs = []
    codonTotals = []
    codonRSCU = []
    relCodonAdapt = []
    
    AADicts = []
    AAFreqs = []
    AATotals = []
    
    nucDicts = []
    nucFreqs = []
    nucTotals = []

    speciesList= []
    
    sa = seqA.NucParams
    for inFile in cmdLine.args.inFile:
        fasta = seqA.FastAreader(inFile)
        keepRunning = 0 #for seq statement to keep running
        runOnce = 0 #for head statement to run once
        
        for head, seq in fasta.readFasta():
            #run only once with first header to get species name
            if runOnce == 0:
                speciesName = seqA.SpeciesFromFasta.getSpeciesFromHeader(head)
                speciesList.append(speciesName)
                runOnce = 1
            
            #run every seq through for analysis
            if keepRunning == 0:
                sa(seq, starts=cmdLine.args.starts)
        
        #save final codon/AA/nuc dictionaries from inFiles
        codonDicts.append(sa.totalCodons)
        codonTotals.append(sum(sa.totalCodons.values()))
        codonFreqDict = sa.totalCodons
        codonFreqDict = {key:codonFreqDict[key]/sum(sa.totalCodons.values())*100 for key in codonFreqDict}
        codonFreqs.append(codonFreqDict)
        
        AADicts.append(sa.totalAAs)
        AATotals.append(sum(sa.totalAAs.values()))
        AAFreqDict = sa.totalAAs
        AAFreqDict = {key:AAFreqDict[key]/sum(sa.totalAAs.values())*100 for key in AAFreqDict}
        AAFreqs.append(AAFreqDict)
        
        nucDicts.append(sa.totalNucs)
        nucTotals.append(sum(sa.totalNucs.values()))
        nucFreqDict = sa.totalNucs
        nucFreqDict = {key:nucFreqDict[key]/sum(sa.totalNucs.values())*100 for key in nucFreqDict}
        nucFreqs.append(nucFreqDict)
        
        sa.rscuValues(sa.totalCodons, sa.totalAAs)
        codonRSCU.append(sa.rscu)
        relCodonAdapt.append(sa.relAdaptOfCodon)
        
        #reset dictionaries so info from following files doesn't add to it
        sa.totalCodons = {codon:0 for codon in sa.rnaCodonTable}
        sa.totalAAs = {aminoAcid:0 for aminoAcid in sa.rnaCodonTable.values()}
        sa.totalNucs = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
        sa.rscu = {codon:0 for codon in sa.rnaCodonTable}
        sa.relAdaptOfCodon = {codon:0 for codon in sa.rnaCodonTable}
    
    inSpecies = []
    for inFile in CommandLine().parser.parse_args().inFile:
        start = 'codingSeqs/'
        end = 'Human'
        inSpecies.append(inFile[inFile.find(start)+len(start):inFile.rfind(end)])
    
    for inFile in cmdLine.args.inFile:
        fasta = seqA.FastAreader(inFile)
        for head, seq in fasta.readFasta():
    
    
    print(codonRSCU, relCodonAdapt)
    DrawCharts.codonHeatMap(codonRSCU, inSpecies, nucFreqs, AAFreqs)

if __name__ == "__main__":
    main()
