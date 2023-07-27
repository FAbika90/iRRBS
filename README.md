# iRRBS
# iRRBS tool
## Overview
3’ ends of RRBS reads overlapping with genomic MspI sites include non-methylated cytosines introduced through end-repair. These cytosines are not recognized by Trim Galore and are therefore not trimmed but considered during methylation calling. To avoid methylation bias we developed iRRBS, which identifies and hides end-repaired cytosines from methylation calling.

## Features
- Detecting whether the input file is single-read or paired-end
- Logging the "Number of unique MspI reads", the "Number of MspI reads" and the "Number of all reads"
- Outputting a BAM file without the biased cytosines


## Usage
To run iRRBS the following input parameters are required in this order:
- infile: path to input sorted BAM file
- chromsizes: path to chrom.sizes file to define the chromosome lengths for a given genome
- genome: path to genome file
-outfile: name for the output file

## Dependencies
- samtools
- bedtools
