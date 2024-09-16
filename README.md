# improve-RRBS tool
## Overview
3’ ends of RRBS reads overlapping with genomic MspI sites include non-methylated cytosines introduced through end-repair. These cytosines are not recognized by Trim Galore and are therefore not trimmed but considered during methylation calling. To avoid methylation bias we developed improve-RRBS, which identifies and hides end-repaired cytosines from methylation calling.

## Features
- Detecting whether the input file is single-read or paired-end
- Logging the "Unique MspI site-aligned reads", the "Total MspI site-aligned reads" and the "Total aligned reads in dataset"
- Outputting a BAM file without the biased cytosines

## Installation
improve-RRBS can be installed by pip:
`pip install iRRBS`

## Usage

`python -m iRRBS.run_irrbs -i <infile> -c <chromsizes> -g <genome> -o <outputfile>`

To run improve-RRBS the following input parameters are required in this order:
- infile (-i): path to input sorted BAM file with an associated index file
- chromsizes (-c): path to chrom.sizes file to define the chromosome lengths for a given genome
- genome (-g): path to genome file
- outfile (-o): name for the output file


## Setup and Test Instructions
This guide walks you through the process of testing the iRRBS (IMPROVE-RRBS) analysis tool using a specific genomic data file. 
Follow these steps to ensure everything is set up and running correctly.
Prerequisites

Before you begin, ensure you have the following installed:

- Python
- samtools
- bedtools
    
### Test Data Preparation

Unzip the rn6 chromosome 20 fasta file

The rn6 (Rattus norvegicus) chromosome 20 fasta file comes compressed. Use the following command to decompress it:

`gunzip rn6_chr20.fa.gz`

### Run iRRBS

With the data prepared, run the iRRBS tool using the command below. Test data can be found in the "test_data" folder. 
This command specifies the input BAM file, chromosome sizes file, genome fasta file, and the output file names.

`python -m iRRBS.run_irrbs -i test_data.bam -c rn6_chr20.chrom.sizes -g rn6_chr20.fa -o test_data_out_test.bam`

### Log File

After running the test, a log file named test_data_out_test.log will be generated. This file contains information about the summary statistics of the analysis.
    
Please ensure you check the log file for any potential issues or to verify the successful completion of the test.
Compare it with the test_data_out.log file from the test_data folder.
    
