"""
Joshua Arribere Aug 16, 2016

Script to call percent hawaiian of a sequencing file. This is
    primarily to automate the read mapping, deduplication steps
    but will not produce an actual plot

Input: inFile.fastq - reads in fastq format, adaptor trimmed

Output: A series of files, see below

run as python pipelineWrapperHawaiian.py inFile.fastq outPrefix
"""
import sys, common, os
from logJosh import Tee

def main(args):
    inFile,outPrefix=args[0:]
    
    #set some parameters
    #bowtieIndex='/data1/genomes/160110_Celegans_rel83/bowtieIndex/160110_Cel_rel83'
    #bowtieIndex='/data4/genomes/171207_Celegans_release90/180402_bowtieIndex/180402_bowtieIndex_rel90'
    #bowtieIndex='/data4/genomes/180402_clone_160110_Celegans_rel83/160110_Celegans_rel83/bowtieIndex/160110_Cel_rel83'
    #bowtieIndex='/data4/genomes/171218_historicalGenomeT2A/bowtieIndex/181110_pd4092Genome'
    #bowtieIndex='/data12/joshua/genomes/171218_historicalGenomeT2A/bowtieIndex/181110_pd4092Genome'
    bowtieIndex='/data9/marissa/genomes/171218_historicalGenomeT2A/bowtieIndex/181110_pd4092Genome'
    ###########################################################
    # Map Reads
    ###########################################################
    print 'Mapping reads %s...'%(inFile)
    os.system('bowtie2 -p 15 -x %s %s -S %s.sam'%(bowtieIndex,inFile,outPrefix)) #-p 15 is 15 cores
    ###sys.exit()
    ###########################################################
    # Convert to BAM and index
    ###########################################################
    print 'Converting to BAM file %s...'%(inFile)
    #190104 added -o to line below due to versioning issues with samtools
    os.system('samtools view -bS %s.sam | samtools sort - -o %s.bam'%(outPrefix,outPrefix))#no need to add .bam extension b/c samtools will anyways
    os.system('samtools index %s.bam'%(outPrefix))
    ###sys.exit()
    #
    ###########################################################
    # Mark Read Duplicates
    ###########################################################
    'Deduplicating reads %s...'%(inFile)
    #picardPath='/home/joshua/programs/picard-2.6.0/picard/build/libs/picard.jar'
    #picardPath='/home/joshua/programs/picard.jar'
    picardPath='/home/marissa/Downloads/picard.jar'
    os.system('java -jar %s MarkDuplicates I=%s.bam O=%s.dedup.bam M=%s.dedup.txt'%(picardPath,outPrefix,outPrefix,outPrefix))
    ###sys.exit()
    #
    ###########################################################
    # Add group names to the file
    ###########################################################
    os.system('java -jar %s AddOrReplaceReadGroups I=%s.dedup.bam O=%s.dedup.good.bam RGSM=dummyName RGLB=dummyName RGPL=illumina RGPU=none VALIDATION_STRINGENCY=LENIENT'%(picardPath,outPrefix,outPrefix))
    ###sys.exit()
    #
    ###########################################################
    # Index the bam file
    ###########################################################
    os.system('samtools index %s.dedup.good.bam'%(outPrefix))
    #
    ###########################################################
    # Annotate variants using CB4856 MMP sequencing
    ###########################################################
    print 'Annotating variants %s...'%(inFile)
    print 'Quitting early not annotating variants'
    sys.exit()
    #hawaiianSet='/data1/160810_Hw/references/CB4856/160811_CB4856.HomoAnd20Ct.vcf'
    #hawaiianSet='/data4/genomes/180402_clone_160110_Celegans_rel83/160811_CB4856.HomoAnd20Ct.vcf'
    #hawaiianSet='/data12/joshua/genomes/180402_clone_160110_Celegans_rel83/160811_CB4856.HomoAnd20Ct.vcf'
    hawaiianSet='/data9/marissa/genomes/180402_clone_160110_Celegans_rel83/160811_CB4856.HomoAnd20Ct.vcf'
    #os.system('java -jar /home/joshua/programs/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data12/joshua/genomes/180402_clone_160110_Celegans_rel83/160110_Celegans_rel83/160110_allChrs.fa -I %s.dedup.good.bam --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles %s -o %s.dedup.good.vcf -allowPotentiallyMisencodedQuals'%(outPrefix,hawaiianSet,outPrefix))
    os.system('java -jar /home/marissa/Downloads/GATK/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar -T HaplotypeCaller -R /data12/joshua/genomes/180402_clone_160110_Celegans_rel83/160110_Celegans_rel83/160110_allChrs.fa -I %s.dedup.good.bam --genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles %s -o %s.dedup.good.vcf -allowPotentiallyMisencodedQuals'%(outPrefix,hawaiianSet,outPrefix))
    #
    ###########################################################
    # Done!
    ###########################################################
    print 'Done with %s.'%(inFile)

if __name__=='__main__':
    Tee()
    main(sys.argv[1:])
