"""
Marissa Glover, March 27, 2020

Writing a script that can make a protein domain schematic.
Want boxes to scale to size of each domain.
Want line to represent regions that aren't part of a domain.

Input: fasta amino acid sequence of each part of the protein
Output: schematic showing protein domains
"""

import sys, math, pyx
from pyx import *
"""
EXAMPLES
rect = path.path(path.moveto(0,0), path.lineto(0,1),
        path.lineto(1,1), path.lineto(1,0), path.closepath())
c.stroke(rect, [deco.filled([color.grey(0.95)])]) #fill color with decorator
c.writePDFfile('rect')

rect2 = path.rect(0, 0, 1, 1)
"""

class CommandLine():
    """Handle command line to determine inFile and outFile.
    Only arguments specified in command line are file names."""
    def __init__(self, inOpts=None):
        """Implement parser to interpret command line argv string using argparse."""
        import argparse
        self.parser = argparse.ArgumentParser(description="Make protein domain schematic.")
        self.parser.add_argument("inFile", action="store", help="Input is tab delimited text file with each protein domain and total number of amino acids. Example: DomainName/t DomainStart/t DomainStop/n Total/t Total#AminoAcids.")
        self.parser.add_argument("outFile", action="store", help="Output name for protein domain schematic")
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

class ProteinDomains:
    """Will open text file as described in CommandLine class and make protein domain schematic.
    Will make a rectangle for each domain and a line for the whole protein sequence."""
    def __init__(self, inFile):
        self.inFile = inFile
        self.header = '' #header from file
        self.domains = {} #{'Domain Name' : [start, stop]}
        self.total = [] #start and stop of protein seq
    
    def parseFile (self):
        with open(self.inFile, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    self.header = line
                elif 'Total' in line:
                    line = line.strip().split('\t')
                    self.total = [int(line[1]), int(line[2])]
                else:
                    line = line.strip().split('\t')
                    self.domains[line[0]] = [int(line[1]), int(line[2])]
        return self.header, self.total, self.domains
    
    def drawDomains (self):
        c = canvas.canvas()
        #draw line for length of protein
        p = path.line(self.total[0], 0.5, self.total[1], 0.5)
        c.stroke(p, [style.linewidth.THICK, trafo.scale(sx=0.1, sy=1)])
        
        for key, value in self.domains.items():
            domainName = key
            domainStart = value[0]
            domainStop = value[1]
            print(domainStart,domainStop, self.total[0], self.total[1])
            
            #following 2 lines don't draw rect to same length as line above
            #rect = path.rect(self.total[0], 0, self.total[1], 1)
            #rect = path.rect(domainStart, 0, domainStop, 1) #scale length to domain size
            
            rect = path.path(path.moveto(domainStart, 0), path.lineto(domainStart, 1), path.lineto(domainStop, 1), path.lineto(domainStop, 0), path.closepath())
            c.stroke(rect, [style.linewidth.Thick, deco.filled([color.cmyk.Violet]), trafo.scale(sx=0.1, sy=1)])
            c.text(1, 0.5, 'Test')
        c.writePDFfile('rectTest')

def main (inCL=None):
    if inCL is None:
        myCmdLine = CommandLine()
        process = ProteinDomains(myCmdLine.args.inFile)
        process.parseFile()
        process.drawDomains()
    else:
        myCmdLine = CommandLine(inCL)
    print (myCmdLine.args)

if __name__ == "__main__":
    main()
