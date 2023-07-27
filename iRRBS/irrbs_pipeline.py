import os
import re
import pysam
import pybedtools
import subprocess
import array
import shutil
import argparse

class IPipeline:
    
    def __init__(self, infile, chromsizes, genome, outfile):
        self.infile = infile
        self.chromsizes = chromsizes
        self.genome = genome
        self.outfile = outfile
        self.PE = 0


    def check_dependencies(self):
        try:
            result_samtools = subprocess.run(['samtools', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            result_bedtools = subprocess.run(['bedtools', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        except FileNotFoundError:
            print("Error: 'samtools' or 'bedtools' not found. Please make sure they are installed and available in the PATH.")
            exit(1)
            
    def check_files_exist(self):
        if not os.path.isfile(self.infile):
            print("Input file does not exist")
            exit(1)
        if not os.path.isfile(self.chromsizes):
            print("Chromosome sizes file does not exist")
            exit(1)
        if not os.path.isfile(self.genome):
            print("Genome file does not exist")
            exit(1)

    def check_input_formats(self):
        if not re.search(r'\.bam$', self.infile):
            print("Input file is not BAM")
            exit(1)
        if not re.search(r'\.bam$', self.outfile):
            print("Output file is not BAM")
            exit(1)

    def check_if_sorted(self):
        try:
            pysam.view('-H', self.infile)
        except pysam.SamtoolsError:
            self.log_error("Input file is not sorted")
            exit(1)


    def paired_check(self):
        input_bamfile = pysam.AlignmentFile(self.infile, "rb")
        for read in input_bamfile.fetch():
            if read.is_paired:
                self.PE = 1
                break
        input_bamfile.close()

    def pair_split(self):
        if self.PE != 0:
            input_bamfile = pysam.AlignmentFile(self.infile, "rb")
            R1_reads = pysam.AlignmentFile(self.infile.replace('.bam', '_R1.bam'), "wb", template=input_bamfile)
            R2_reads = pysam.AlignmentFile(self.infile.replace('.bam', '_R2.bam'), "wb", template=input_bamfile)

            for read in input_bamfile.fetch():
                if read.is_read1:
                    R1_reads.write(read)
                if read.is_read2:
                    R2_reads.write(read)

            R1_reads.close()
            R2_reads.close()
            input_bamfile.close()
        else:  # SE
            shutil.copy(self.infile, self.infile.replace('.bam', '_R1.bam'))

    def block_find(self): 
        # Defining blocks
        pybedtools.BedTool(self.infile.replace('.bam', '_R1.bam')).bam_to_bed().cut([0, 1, 2, 5]).saveas(
            self.infile.replace('.bam', '_R1.bed'))
        cmd = "sed -i 's/\t/\t0\t0\t/3' " + self.infile.replace('.bam', '_R1.bed')
        subprocess.call(cmd, shell=True)
        cmd = "sort -V -u " + self.infile.replace('.bam', '_R1.bed') + " -o" + self.infile.replace('.bam', '_R1.bed')
        subprocess.call(cmd, shell=True)
        pybedtools.BedTool(self.infile.replace('.bam', '_R1.bed')).slop(g=self.chromsizes, s=True, l=0, r=2).saveas(
            self.infile.replace('.bam', '_R1.bed'))
        cmd = "bedtools getfasta " + " -s " " -fi " + self.genome + " -bed " + self.infile.replace(
            '.bam', '_R1.bed') + " -bedOut" + " > " + self.infile.replace('.bam', '_R1_seq.bed')
        subprocess.call(cmd, shell=True)
        pybedtools.BedTool(self.infile.replace('.bam', '_R1_seq.bed')).filter(
            lambda x: re.findall(r'CCGG.{0,2}$', x[6], flags=re.IGNORECASE)).saveas(
            self.infile.replace('.bam', '_R1_seq2.bed'))
        pybedtools.BedTool(self.infile.replace('.bam', '_R1_seq2.bed')).slop(g=self.chromsizes, s=True, l=0,
                                                                            r=-2).saveas(
            self.infile.replace('.bam', '_blocks.bed'))

    def msp1_split(self):
        pybedtools.BedTool(self.infile.replace('.bam', '_R1.bam')).intersect(pybedtools.BedTool(self.infile.replace('.bam', '_blocks.bed')), s=True, f=1, F=1).saveas(self.infile.replace('.bam', '_msp1.bam'))
        pybedtools.BedTool(self.infile.replace('.bam', '_R1.bam')).intersect(pybedtools.BedTool(self.infile.replace('.bam', '_blocks.bed')), s=True, v=True, f=1, F=1).saveas(self.infile.replace('.bam', '_msp1neg.bam'))

    def logging(self):
        # Logging
        with open(self.outfile.replace('.bam', '.log'), 'w') as logs:
            logs.write('Number of unique MspI reads:\n')
            logs.write(str(len(pybedtools.BedTool(self.infile.replace('.bam', '_blocks.bed')))))
            logs.write('\n')
            logs.write('Number of MspI reads:\n')
            logs.write(str(pysam.view('-c', '-F', '4', self.infile.replace('.bam', '_msp1.bam'))))
            logs.write('Number of all reads:\n')
            logs.write(str(pysam.view('-c', '-F', '4', self.infile)))
        # Removing temp files  
        deletefiles = [self.infile.replace('.bam', '_R1.bed'), self.infile.replace('.bam', '_R1_seq.bed'), self.infile.replace('.bam', '_R1_seq2.bed'), self.infile.replace('.bam', '_blocks.bed'), self.infile.replace('.bam', '_R1.bam')]
        for line in deletefiles:
            os.remove(line)

    def msp1_clip(self):
        #MspI soft clipping
        pysam.index(self.infile.replace('.bam', '_msp1.bam'))
        msp1bamfile = pysam.AlignmentFile(self.infile.replace('.bam', '_msp1.bam'), "rb")
        ModReads = pysam.AlignmentFile(self.infile.replace('.bam', '_msp1_mod.sam'), "wb", template=msp1bamfile)
        
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
        os.remove(self.infile.replace('.bam', '_msp1.bam'))
        
    def file_merge(self):
        if self.PE != 0:
            bamfiles = [self.infile.replace('.bam', '_R2.bam'), self.infile.replace('.bam', '_msp1neg.bam'), self.infile.replace('.bam', '_msp1_mod.sam')]
        else:
            bamfiles = [self.infile.replace('.bam', '_msp1neg.bam'), self.infile.replace('.bam', '_msp1_mod.sam')]
        pysam.merge("-f", "-o", self.outfile, *bamfiles)
        for line in bamfiles:
            os.remove(line)


def main():
    pass
    
if __name__ == '__main__':
    main()

