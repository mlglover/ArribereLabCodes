"""
Marissa Glover, April 9, 2020

Script to make a heat map of codon frequencies.
Can look at 1 or more species to compare.

Input: dictionaries created by sequenceAnalysis.py
Output: heat map with codons on x-axis and species on y-axis.
"""

import sys, sequenceAnalysis as seqA, matplotlib.pyplot as plt, numpy as np, pandas as pd, seaborn as sns

class CommandLine():
    def __init__(self, inOpts=None):
        import argparse
        self.parser = argparse.ArgumentParser(description='Make heat map of codon frequencies')
        self.parser.add_argument('inFile', action='store', nargs='?', help='Input fasta files to compare to each other')
        self.parser.add_argument('outFile', action'store', help='Output file name for results')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

################################################################################

class MakeHeatmap:
    def __init__(self, dictFromSeqA):
        self.codonComp = dictFromSeqA #dictionary of codon frequencies
    
    def heatmap(self):
        pass

def main(inCL=None):
    if inCL is None:
        cmdLine = CommandLine()
    else:
        cmdLine = CommandLine(inCL)

    getDict = seqA.NucParams(seqA.FastAreader)
    plot = MakeHeatmap(getDict)
    plot.heatmap()
