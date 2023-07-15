#!/usr/bin/env python

import sys
import os
import argparse
import re
import pysam
import pybedtools
import subprocess
import array
import shutil

def check(infile, chromsizes, genome, outfile):
    # Check if input files exist
    if not os.path.isfile(infile):
        print("Input file does not exist")
        sys.exit(1)
    if not os.path.isfile(chromsizes):
        print("Chromosome sizes file does not exist")
        sys.exit(1)
    if not os.path.isfile(genome):
        print("Genome file does not exist")
        sys.exit(1)
    # Check if input file is BAM
    if not re.search(r'\.bam$', infile):
        print("Input file is not BAM")
        sys.exit(1)
    # Check if output file is BAM
    if not re.search(r'\.bam$', outfile):
        print("Output file is not BAM")
        sys.exit(1)
    # Check if input file is sorted
    try:
        pysam.view('-H', infile)
    except pysam.SamtoolsError:
        print("Input file is not sorted")
        sys.exit(1)
        
def pairedcheck(infile):
    inputbamfile = pysam.AlignmentFile(infile, "rb")
    PE = 0
    for read in inputbamfile.fetch():
        if read.is_paired:
            PE += 1
            break
    inputbamfile.close()
    return PE
                                            
def pairsplit(PE, infile):
    # Read-pair splitting
    #PE
    if PE != 0:
        inputbamfile = pysam.AlignmentFile(infile, "rb")
        R1reads = pysam.AlignmentFile(infile.replace('.bam', '_R1.bam'), "wb", template=inputbamfile)
        R2reads = pysam.AlignmentFile(infile.replace('.bam', '_R2.bam'), "wb", template=inputbamfile)             
        for read in inputbamfile.fetch():
            if read.is_read1:
                R1reads.write(read)
            if read.is_read2:
                R2reads.write(read)

        R1reads.close()
        R2reads.close()
        inputbamfile.close()
    else: #SE
        shutil.copy(infile, infile.replace('.bam', '_R1.bam'))
    
def blockfind(infile, chromsizes, genome):
    #defining blocks
    pybedtools.BedTool(infile.replace('.bam', '_R1.bam')).bam_to_bed().cut([0, 1, 2, 5]).saveas(infile.replace('.bam', '_R1.bed'))
    cmd = "sed -i 's/\t/\t0\t0\t/3' " + infile.replace('.bam', '_R1.bed') 
    subprocess.call(cmd, shell=True)
    cmd = "sort -V -u " + infile.replace('.bam', '_R1.bed') + " -o" + infile.replace('.bam', '_R1.bed') 
    subprocess.call(cmd, shell=True)
    pybedtools.BedTool(infile.replace('.bam', '_R1.bed')).slop(g=chromsizes, s=True, l=0, r=2).saveas(infile.replace('.bam', '_R1.bed'))
    cmd = "bedtools getfasta " + " -s " " -fi " +  genome +  " -bed " +  infile.replace('.bam', '_R1.bed') + " -bedOut" + " > " + infile.replace('.bam', '_R1_seq.bed')
    subprocess.call(cmd, shell=True)
    pybedtools.BedTool(infile.replace('.bam', '_R1_seq.bed')).filter(lambda x: re.findall(r'CCGG.{0,2}$', x[6], flags=re.IGNORECASE)).saveas(infile.replace('.bam', '_R1_seq2.bed')) 
    pybedtools.BedTool(infile.replace('.bam', '_R1_seq2.bed')).slop(g=chromsizes, s=True, l=0, r=-2).saveas(infile.replace('.bam', '_blocks.bed'))

def msp1split(infile):
    pybedtools.BedTool(infile.replace('.bam', '_R1.bam')).intersect(pybedtools.BedTool(infile.replace('.bam', '_blocks.bed')), s=True, f=1, F=1).saveas(infile.replace('.bam', '_msp1.bam'))
    pybedtools.BedTool(infile.replace('.bam', '_R1.bam')).intersect(pybedtools.BedTool(infile.replace('.bam', '_blocks.bed')), s=True, v=True, f=1, F=1).saveas(infile.replace('.bam', '_msp1neg.bam'))

def logging(infile, outfile):
    #logging
    with open(outfile.replace('.bam', '.log'), 'w') as logs:
        logs.write('Number of unique MspI reads:\n')
        logs.write(str(len(pybedtools.BedTool(infile.replace('.bam', '_blocks.bed')))))
        logs.write('\n')
        logs.write('Number of MspI reads:\n')
        logs.write(str(pysam.view('-c', '-F', '4', infile.replace('.bam', '_msp1.bam'))))
        logs.write('Number of all reads:\n')
        logs.write(str(pysam.view('-c', '-F', '4', infile)))
    #removing temp files  
    deletefiles = [infile.replace('.bam', '_R1.bed'), infile.replace('.bam', '_R1_seq.bed'), infile.replace('.bam', '_R1_seq2.bed'), infile.replace('.bam', '_blocks.bed'), infile.replace('.bam', '_R1.bam')]
    for line in deletefiles:
    	os.remove(line)

def msp1clip(infile, genome):
    #MspI soft clipping
    pysam.index(infile.replace('.bam', '_msp1.bam'))
    msp1bamfile = pysam.AlignmentFile(infile.replace('.bam', '_msp1.bam'), "rb")
    ModReads = pysam.AlignmentFile(infile.replace('.bam', '_msp1_mod.sam'), "wb", template=msp1bamfile)
    
    for read in msp1bamfile.fetch():
        if read.is_forward:
            quals = read.query_qualities
            read.query_sequence = read.query_sequence[:-3] + 'NNN'
            read.query_qualities = quals[:-3] + array.array('B', [0,0,0])
            read.set_tags([('NM', read.tags[0][1]),('MD', read.tags[1][1]),
                           ('XM', read.tags[2][1][:-3] + '...'),
                           ('XR', read.tags[3][1]),('XG', read.tags[4][1])])
        if read.is_reverse:
            quals = read.query_qualities
            read.query_sequence =  'NNN' + read.query_sequence[3:]
            read.query_qualities = array.array('B', [0,0,0]) + quals[3:]
            read.set_tags([('NM', read.tags[0][1]),('MD', read.tags[1][1]),
                           ('XM', '...' + read.tags[2][1][3:]),
                           ('XR', read.tags[3][1]),('XG', read.tags[4][1])])
        ModReads.write(read)
    ModReads.close()
    msp1bamfile.close()
    os.remove(infile.replace('.bam', '_msp1.bam'))
    
def filemerge(PE, infile, outfile):
    if PE != 0:
        bamfiles = [infile.replace('.bam', '_R2.bam'), infile.replace('.bam', '_msp1neg.bam'), infile.replace('.bam', '_msp1_mod.sam')]
    else:
        bamfiles = [infile.replace('.bam', '_msp1neg.bam'), infile.replace('.bam', '_msp1_mod.sam')]
    pysam.merge("-f","-o",outfile,*bamfiles)
    for line in bamfiles:
    	os.remove(line)

class suppress_stdout_stderr(object):
    def __init__(self):
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        self.save_fds = [os.dup(1), os.dup(2)]
    def __enter__(self):
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)
    def __exit__(self, *_):
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        for fd in self.null_fds + self.save_fds:
            os.close(fd)
             	
def main():
    parser = argparse.ArgumentParser(description='MspI-seq pipeline')
    parser.add_argument('-i', '--infile', help='Input BAM file', required=True)
    parser.add_argument('-c', '--chromsizes', help='Chromosome sizes file', required=True)
    parser.add_argument('-g', '--genome', help='Genome file', required=True)
    parser.add_argument('-o', '--outfile', help='Output BAM file', required=True)
    args = parser.parse_args()
    check(args.infile, args.chromsizes, args.genome, args.outfile)
    PE=pairedcheck(args.infile)
    with suppress_stdout_stderr():
    	pairsplit(PE, args.infile)
    	blockfind(args.infile, args.chromsizes, args.genome)
    	msp1split(args.infile)
    	logging(args.infile, args.outfile)
    	msp1clip(args.infile, args.genome)
    	filemerge(PE, args.infile, args.outfile)
if __name__ == '__main__':
    main()
