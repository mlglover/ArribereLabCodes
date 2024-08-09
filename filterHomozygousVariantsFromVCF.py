"""
Joshua Arribere Aug 11, 2016

Script to pull out only the homozygous variants from a vcf file

run as pythong inFile.vcf outPrefix
EDIT: Aug 11, 2016 JOSH added a minimum read count filter as well
"""
import sys, common
from logJosh import Tee

def main(args):
    inFile,outPrefix=args[0:]
    N=5
    print 'Requiring %s reads as well!'%(N)
    #cut=40
    #print 'Requiring RMSQ of at least %s'%cut
    with open(inFile,'r') as f:
        with open(outPrefix,'w') as g:
            for line in f:
                if line.startswith('#'):
                    g.write(line)
                else:
                    #if 'AF=1.00' in line:
                    if '1/1' in line:
                        line2=line.strip().split('\t')
                        cts=sum(map(float,line2[-1].split(':')[1].split(',')))
                        #RMSQ=float(line2[7].split('MQ=')[1].split(';')[0])
                        #if RMSQ>=cut:
                        if cts >=N:
                            g.write(line)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
