#!/usr/bin/env python3

"""
Marissa Glover
Script to perform quality control step for indexing of Illumina libraries.
Contains 3 classes: CommandLine, FastQreader, IndexQC

Inputs: rawReads_I1.fastq, rawReads_I2.fastq
Output: heat map with frequency of reads for each i7/i5 index combination,
        with frequency of invalid index combination frequency printed below

Need to determine the fraction of reads that have each index. Will take files
with raw reads, strip down file to only contain read sequences, and create
dictionary with {keys,values} of {indexes,number of reads with index}

Run as: python illuminaIndexQualityControl.py inFile1 inFile2 outFile
"""

import sys
from itertools import islice
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

class CommandLine():
    """Handle the command line to determine inFiles and outFile.
    Only arguments that are specified in the command line are the file names."""
    def __init__(self, inOpts=None):
        """Implement a parser to interpret command line argv string using argparse."""
        import argparse
        self.parser = argparse.ArgumentParser(description="Process index fastq files.")
        self.parser.add_argument("inFile1", action="store", nargs="+", help="i7 index input file name")
        self.parser.add_argument("inFile2", action="store", nargs="+", help="i5 index input file name")
        self.parser.add_argument("outFile", action="store",
        help="output file name, must be one of the following:.eps, .pdf, .pgf, .png, .ps, .raw, .rgba, .svg, .svgz")

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

##########################################################################################

class FastQreader:
    """Define objects to read two fastq files and only extract lines containing sequence
    reads. Uses a generator to pass lines one at a time to another class.

    Requires that fileNames contain I1 and I2 for the index 1 and index 2 reads.

    Instantiation:
        fastq = FastQreader ("inFile1", "inFile2")
    Usage:
        for indexes in fastq.parseFiles():
            print (indexes)"""
    def __init__ (self, fileName1, fileName2):
        """Constructor: saves attributes fileName1 and fileName2."""
        if 'I1' in fileName1 and 'I2' in fileName2:
            self.fileName1 = fileName1
            self.fileName2 = fileName2
        else:
            self.fileName1 = fileName2
            self.fileName2 = fileName1

    def parseFiles (self):
        """Opens fastq files simultaneously, then compiles each index combo from the
        line-matched index files and passes them one at a time."""
        with open(self.fileName1,'r') as f1, open(self.fileName2, 'r') as f2:
            for line in zip(islice(f1,1,None,4), islice(f2,1,None,4)):
                yield line

##########################################################################################

class IndexQC:
    """Determine the frequency of index hopping from high-throughput sequencing flow cell.
    Plot valid index combination frequencies on a heat map and print other information
    about invalid combinations below the plot.

    Class attributes will allow information for each index combination that is passed
    through to be compiled.

    Instantiation:
        iqc = IndexQC(fastq)
    Usage:
        for indexes in FastQreader.parseFiles(fastq):
            iqc = IndexQC(indexes)
            iqc.countIndexes()
        iqc = IndexQC
        iqc.validOrInvalidCombos(iqc())
        iqc.plotResults(iqc())"""

    #dictionary of all valid index combinations
    i7i5Combos = {
        ('ATTACTCG', 'AGGCTATA'): 0, ('TCCGGAGA', 'AGGCTATA'): 0, ('CGCTCATT', 'AGGCTATA'): 0,
        ('GAGATTCC', 'AGGCTATA'): 0, ('ATTCAGAA', 'AGGCTATA'): 0, ('GAATTCGT', 'AGGCTATA'): 0,
        ('CTGAAGCT', 'AGGCTATA'): 0, ('TAATGCGC', 'AGGCTATA'): 0, ('CGGCTATG', 'AGGCTATA'): 0,
        ('TCCGCGAA', 'AGGCTATA'): 0, ('TCTCGCGC', 'AGGCTATA'): 0, ('AGCGATAG', 'AGGCTATA'): 0,
        ('ATTACTCG', 'GCCTCTAT'): 0, ('TCCGGAGA', 'GCCTCTAT'): 0, ('CGCTCATT', 'GCCTCTAT'): 0,
        ('GAGATTCC', 'GCCTCTAT'): 0, ('ATTCAGAA', 'GCCTCTAT'): 0, ('GAATTCGT', 'GCCTCTAT'): 0,
        ('CTGAAGCT', 'GCCTCTAT'): 0, ('TAATGCGC', 'GCCTCTAT'): 0, ('CGGCTATG', 'GCCTCTAT'): 0,
        ('TCCGCGAA', 'GCCTCTAT'): 0, ('TCTCGCGC', 'GCCTCTAT'): 0, ('AGCGATAG', 'GCCTCTAT'): 0,
        ('ATTACTCG', 'AGGATAGG'): 0, ('TCCGGAGA', 'AGGATAGG'): 0, ('CGCTCATT', 'AGGATAGG'): 0,
        ('GAGATTCC', 'AGGATAGG'): 0, ('ATTCAGAA', 'AGGATAGG'): 0, ('GAATTCGT', 'AGGATAGG'): 0,
        ('CTGAAGCT', 'AGGATAGG'): 0, ('TAATGCGC', 'AGGATAGG'): 0, ('CGGCTATG', 'AGGATAGG'): 0,
        ('TCCGCGAA', 'AGGATAGG'): 0, ('TCTCGCGC', 'AGGATAGG'): 0, ('AGCGATAG', 'AGGATAGG'): 0,
        ('ATTACTCG', 'TCAGAGCC'): 0, ('TCCGGAGA', 'TCAGAGCC'): 0, ('CGCTCATT', 'TCAGAGCC'): 0,
        ('GAGATTCC', 'TCAGAGCC'): 0, ('ATTCAGAA', 'TCAGAGCC'): 0, ('GAATTCGT', 'TCAGAGCC'): 0,
        ('CTGAAGCT', 'TCAGAGCC'): 0, ('TAATGCGC', 'TCAGAGCC'): 0, ('CGGCTATG', 'TCAGAGCC'): 0,
        ('TCCGCGAA', 'TCAGAGCC'): 0, ('TCTCGCGC', 'TCAGAGCC'): 0, ('AGCGATAG', 'TCAGAGCC'): 0,
        ('ATTACTCG', 'CTTCGCCT'): 0, ('TCCGGAGA', 'CTTCGCCT'): 0, ('CGCTCATT', 'CTTCGCCT'): 0,
        ('GAGATTCC', 'CTTCGCCT'): 0, ('ATTCAGAA', 'CTTCGCCT'): 0, ('GAATTCGT', 'CTTCGCCT'): 0,
        ('CTGAAGCT', 'CTTCGCCT'): 0, ('TAATGCGC', 'CTTCGCCT'): 0, ('CGGCTATG', 'CTTCGCCT'): 0,
        ('TCCGCGAA', 'CTTCGCCT'): 0, ('TCTCGCGC', 'CTTCGCCT'): 0, ('AGCGATAG', 'CTTCGCCT'): 0,
        ('ATTACTCG', 'TAAGATTA'): 0, ('TCCGGAGA', 'TAAGATTA'): 0, ('CGCTCATT', 'TAAGATTA'): 0,
        ('GAGATTCC', 'TAAGATTA'): 0, ('ATTCAGAA', 'TAAGATTA'): 0, ('GAATTCGT', 'TAAGATTA'): 0,
        ('CTGAAGCT', 'TAAGATTA'): 0, ('TAATGCGC', 'TAAGATTA'): 0, ('CGGCTATG', 'TAAGATTA'): 0,
        ('TCCGCGAA', 'TAAGATTA'): 0, ('TCTCGCGC', 'TAAGATTA'): 0, ('AGCGATAG', 'TAAGATTA'): 0,
        ('ATTACTCG', 'ACGTCCTG'): 0, ('TCCGGAGA', 'ACGTCCTG'): 0, ('CGCTCATT', 'ACGTCCTG'): 0,
        ('GAGATTCC', 'ACGTCCTG'): 0, ('ATTCAGAA', 'ACGTCCTG'): 0, ('GAATTCGT', 'ACGTCCTG'): 0,
        ('CTGAAGCT', 'ACGTCCTG'): 0, ('TAATGCGC', 'ACGTCCTG'): 0, ('CGGCTATG', 'ACGTCCTG'): 0,
        ('TCCGCGAA', 'ACGTCCTG'): 0, ('TCTCGCGC', 'ACGTCCTG'): 0, ('AGCGATAG', 'ACGTCCTG'): 0,
        ('ATTACTCG', 'GTCAGTAC'): 0, ('TCCGGAGA', 'GTCAGTAC'): 0, ('CGCTCATT', 'GTCAGTAC'): 0,
        ('GAGATTCC', 'GTCAGTAC'): 0, ('ATTCAGAA', 'GTCAGTAC'): 0, ('GAATTCGT', 'GTCAGTAC'): 0,
        ('CTGAAGCT', 'GTCAGTAC'): 0, ('TAATGCGC', 'GTCAGTAC'): 0, ('CGGCTATG', 'GTCAGTAC'): 0,
        ('TCCGCGAA', 'GTCAGTAC'): 0, ('TCTCGCGC', 'GTCAGTAC'): 0, ('AGCGATAG', 'GTCAGTAC'): 0}

    indexFreqDict = {} #count each index combination seen
    invalidCombos = {} #add invalid combinations not in the i7i5Combos dictionary here
    comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'} #for reverse complement of i5 indexes
    totalCounts = 0 #total number of reads
    totalCountsi7i5 = 0 #total number of valid index reads
    totalInvalid = 0 #total number of invalid index reads

    def __init__ (self, indexLine=''):
        """Constructor: save attribute of index combination from fastq parsing."""
        self.indexLine = indexLine

    def countIndexes (self):
        """Determines frequency of each index combination.
        Valid and invalid combos will all be added to one dictionary."""
        indexList = list(self.indexLine) #convert tuple of indexes to list to make mutable
        indexList[0] = indexList[0].strip('\n')

        #reverse complement of i5 indexes
        indexList[1] = ''.join(IndexQC.comp.get(base,base) for base in reversed(indexList[1])).strip('\n')

        #count index combination each time it is seen
        #convert back to tuple to use indexFreqDict
        if tuple(indexList) in IndexQC.indexFreqDict:
            IndexQC.indexFreqDict[tuple(indexList)] += 1
        else:
            IndexQC.indexFreqDict[tuple(indexList)] = 1
        return IndexQC.indexFreqDict

    def validOrInvalidCombos (self):
        """Separates valid index combiantions from invalid.
        Invalid includes any reads with >=1 N or mismatched base.
        Determines total number of reads, total number of valid combos, and total
        number of invalid combos."""
        #only add valid combinations to dictionary
        for key in IndexQC.indexFreqDict:
            if key in IndexQC.i7i5Combos:
                IndexQC.i7i5Combos[key] += IndexQC.indexFreqDict[key]

        #only add invalid combinations to dictionary
        for key in IndexQC.indexFreqDict:
            if key not in IndexQC.i7i5Combos:
                IndexQC.invalidCombos.update({key:IndexQC.indexFreqDict[key]})

        IndexQC.totalCounts = sum(IndexQC.indexFreqDict.values())
        IndexQC.totalCountsi7i5 = sum(IndexQC.i7i5Combos.values())
        IndexQC.totalInvalid = sum(IndexQC.invalidCombos.values())
        return IndexQC.totalCounts, IndexQC.totalCountsi7i5, IndexQC.totalInvalid

    def plotResults (self):
        """Plots results on heat map with i7 indexes on y-axis and i5 indexes on x-axis.
        Each box will contain the percent of total reads containing that index combo and
        will be colored accordingly. The total number of reads, total number of invalid
        combos, and percent of invalid combos will be displayed on a table below the chart."""
        #make series from valid index to use for heat map
        ser = pd.Series(list((line/IndexQC.totalCountsi7i5)*100
                        for line in IndexQC.i7i5Combos.values()),
        index = pd.MultiIndex.from_tuples(IndexQC.i7i5Combos.keys()))

        #create heat map from series annotated with frequency of combinations
        #in the squares
        df = ser.unstack().fillna(0)
        df.shape #get shape of data from data frame
        sns.set(font_scale=1.1)
        sns.heatmap(df, annot=True, linewidths=.5) #use unstacked data frame for heatmap

        #label plot and axes
        plt.xlabel('i5 Indexes', fontsize=16, labelpad=15)
        plt.ylabel('i7 Indexes', fontsize=16, labelpad=15)
        plt.title('Frequency of Valid Indexes', fontsize=20)

        #create array to put extra information into table below chart
        invalidArray = ([['Total # of reads', IndexQC.totalCounts],
                        ['# of reads with invalid combos', IndexQC.totalInvalid],
                        ['% of total reads with invalid combos',
                        '{:.2f}%'.format((IndexQC.totalInvalid/IndexQC.totalCounts)*100)]])

        #make table below heat map
        table = plt.table(cellText=invalidArray, colWidths=[0.3,0.1], cellLoc='center', bbox=[0.0, -0.48, 0.6, 0.15])
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        plt.tight_layout()
        plt.subplots_adjust(bottom = 0.3) #make room for table
        plt.gcf().set_size_inches(11,8.5, forward=True) #adjust so table isn't over heat map

        #save file using outFile from command line and show the output
        plt.savefig(CommandLine().parser.parse_args().outFile, dpi=300)
        plt.show()

##########################################################################################

def main (inCL=None):
    """Determines the frequency of each valid index combination from I1 and I2 reads
    off Illumina machine. Parses the command line for inFiles and outFile, parses the
    fastq files to extract sequence information, determines frequency of each valid
    index combination, determines total frequency of invalid combinations, and finally
    plots a heat map with a table."""
    if inCL is None:
        myCmdLine = CommandLine()
    else:
        myCmdLine = CommandLine(inCL)

    fastq = FastQreader(myCmdLine.args.inFile1[0], myCmdLine.args.inFile2[0])
    for indexes in FastQreader.parseFiles(fastq): #pass each index combo one at a time
        iqc = IndexQC(indexes)
        iqc.countIndexes()
    iqc = IndexQC
    iqc.validOrInvalidCombos(iqc())
    iqc.plotResults(iqc())

if __name__ == '__main__':
    main()
