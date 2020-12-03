"""
Marissa Glover, April 6, 2020

Script to analyze and compare sequences.

Input: fasta file of genome sequences
Output: nucleotide composition, number ORFs, length of ORFs
"""

import sys, pyx
from pyx import *
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

class FastAreader:
    def __init__(self, fileName=''):
        self.fileName = fileName
    
    def doOpen(self):
        if self.fileName is '':
            return sys.stdin
        else:
            return open(self.fileName)
    
    def readFasta(self):
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            
            #skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()
            
            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header, sequence

################################################################################

class OrfFinder():

    stopCodons = []
    startCodons = []
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    def __init__(self, seq):
        self.seq = seq
        self.orfs = []

    def findOrfs(self):
        starts = []
        stops = []
        for frame in range(0,3):
            for i in range(frame, len(self.seq), 3):
                codon = self.seq[i:i+3]
                if codon in OrfFinder.startCodons:
                    starts.append(i)
                elif codon in OrfFinder.stopCodons:
                    stops.append(i)

################################################################################

class NucParams:
    """Characterize a genome by nucleotide, codon, and amino acid compositions.
    This class will parse sequences to return information about a genome.
    The nucleotide sequence can be given in DNA or RNA.
    The sequence will be converted to RNA to determine amino acid and codon compositions.
    Non-DNA/RNA letters will not be counted.
    Codons with non-DNA/RNA letters will not be counted.
    """
    rnaCodonTable = {
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    
    synCodons = {
        'C': ['UGU', 'UGC'],
        'D': ['GAU', 'GAC'],
        'S': ['UCU', 'UCG', 'UCA', 'UCC', 'AGC', 'AGU'],
        'Q': ['CAA', 'CAG'],
        'M': ['AUG'],
        'N': ['AAC', 'AAU'],
        'P': ['CCU', 'CCG', 'CCA', 'CCC'],
        'K': ['AAG', 'AAA'],
        '-': ['UAG', 'UGA', 'UAA'],
        'T': ['ACC', 'ACA', 'ACG', 'ACU'],
        'F': ['UUU', 'UUC'],
        'A': ['GCA', 'GCC', 'GCG', 'GCU'],
        'G': ['GGU', 'GGG', 'GGA', 'GGC'],
        'I': ['AUC', 'AUA', 'AUU'],
        'L': ['UUA', 'UUG', 'CUC', 'CUU', 'CUG', 'CUA'],
        'H': ['CAU', 'CAC'],
        'R': ['CGA', 'CGC', 'CGG', 'CGU', 'AGG', 'AGA'],
        'W': ['UGG'],
        'V': ['GUA', 'GUC', 'GUG', 'GUU'],
        'E': ['GAG', 'GAA'],
        'Y': ['UAU', 'UAC']}
    
    #amino acid, codon, and nucleotide frequency for all sequences in fasta file
    totalAAs = {aminoAcid:0 for aminoAcid in rnaCodonTable.values()}
    dnaTotalCodons = {codon:0 for codon in dnaCodonTable}
    totalCodons = {codon:0 for codon in rnaCodonTable}
    totalNucs = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
    
    rscu = {codon:0 for codon in rnaCodonTable}
    relAdaptOfCodon = {codon:0 for codon in rnaCodonTable}
    
    def __init__ (self, inString='', starts=None):
        """Implement dictionaries.
        Pass input to addSequence function to parse sequences and add to the dictionaries."""
        
        """Multiple start codons in command line isn't working yet..."""
        self.starts = starts #list of starts from cmd line or ATG default
        self.aaComp = {aminoAcid:0 for aminoAcid in NucParams.rnaCodonTable.values()}
        self.nucComp = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
        self.dnaCodonComp = {codon:0 for codon in NucParams.dnaCodonTable}
        self.codonComp = {codon:0 for codon in NucParams.rnaCodonTable}
        self.addSequence(inString)
    
    def addSequence (self, inSeq):
        """Parse sequences and update dictionaries for each sequence passed through.
        All work is done here, other functions return dictionaries."""
        for start in self.starts:
            trimmed = inSeq[inSeq.find(start):] #trims off seq before start
            
            for nuc in trimmed:
                if nuc in self.nucComp: #only counts valid nucleotides, including N's
                    self.nucComp[nuc] += 1 #parse and count nucleotides
                    NucParams.totalNucs[nuc] += 1
            
            codons = [trimmed[i:i+3] for i in range(0, len(trimmed), 3)] #codon is every 3 nucleotides
            for codon in codons:
                if codon in NucParams.dnaCodonTable:
                    self.dnaCodonComp[codon] += 1 #parse and count DNA codons
                    NucParams.dnaTotalCodons[codon] += 1
                rna = codon.replace('T','U') #convert DNA to RNA
                if rna in NucParams.rnaCodonTable:
                    self.codonComp[rna] += 1 #parse and count RNA codons
                    NucParams.totalCodons[rna] += 1
                    self.aaComp[NucParams.rnaCodonTable.get(rna)] += 1 #translate and count amino acids
                    NucParams.totalAAs[NucParams.rnaCodonTable.get(rna)] += 1
    
    def aaComposition(self):
        """Get the dictionary for aaComposition that was updated in addSequence."""
        return self.aaComp
    
    def nucComposition(self):
        """Get the dictionary for nucComposition that was updated in addSequence."""
        return self.nucComp
    
    def codonComposition(self):
        """Get the dictionary for codonComposition that was updated in addSequence."""
        return self.codonComp
    
    def nucCount(self):
        """Get the sum of the values in the nucComp dictionary for total count of nuclotides."""
        return sum(self.nucComp.values())
    
    def rscuValues(totalCodonDict, totalAADict):
        """Determine relative synonymous codon values from genome."""
        for aminoAcid in NucParams.synCodons:
            rscusForAA = []
            codons = NucParams.synCodons[aminoAcid]
            for codon in codons:
                codonFreq = totalCodonDict[codon]
                numSynCodons = len(codons)
                synCodonUsage = totalAADict[aminoAcid]
                rscuValue = codonFreq / (synCodonUsage / numSynCodons)
                NucParams.rscu[codon] += rscuValue
                rscusForAA.append(rscuValue)
            
            for codon in codons:
                maxRSCU = max(rscusForAA)
                NucParams.relAdaptOfCodon[codon] += (NucParams.rscu[codon] / maxRSCU)

################################################################################

class CAI:
    """Determines codon adaptation index for each gene.
    Attributes for each genome."""
    inputRscu = {codon:0 for codon in NucParams.rnaCodonTable}
    
    def __init__(self, inString='', inputRscus):
        """Attributes for each gene."""
        self.inputRscus = inputRscus
        self.relAdaptOfCodons = {codon:0 for codon in NucParams.rnaCodonTable}
        self.geneLength
        self.caiOfGene = 0
        self.getCAI = getCAI(inString)
    
    def rscuValues(totalCodonDict, totalAADict):
        """Determine relative synonymous codon values from genome."""
        for aminoAcid in NucParams.synCodons:
            rscusForAA = []
            codons = NucParams.synCodons[aminoAcid]
            for codon in codons:
                codonFreq = totalCodonDict[codon]
                numSynCodons = len(codons)
                synCodonUsage = totalAADict[aminoAcid]
                rscuValue = codonFreq / (synCodonUsage / numSynCodons)
                CAI.rscu[codon] += rscuValue
                rscusForAA.append(rscuValue)
            
            for codon in codons:
                maxRSCU = max(rscusForAA)
                CAI.relAdaptOfCodon[codon] += (CAI.rscu[codon] / maxRSCU)
    
    def getCAI(self):
        """Determine CAI for a gene."""
        

################################################################################

class SpeciesFromFasta:
    def getSpeciesFromHeader (header):
        start = 'Organism:'
        end = '|Strain'
        species = header[header.find(start)+len(start):header.rfind(end)]
        return species
