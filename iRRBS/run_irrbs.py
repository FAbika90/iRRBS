import argparse
from iRRBS.irrbs_pipeline import IPipeline 
from iRRBS.suppress import Suppress


def main():
    parser = argparse.ArgumentParser(description='iRRBS pipeline')
    parser.add_argument('-i', '--infile', help='Input BAM file', required=True)
    parser.add_argument('-c', '--chromsizes', help='Chromosome sizes file', required=True)
    parser.add_argument('-g', '--genome', help='Genome file', required=True)
    parser.add_argument('-o', '--outfile', help='Output BAM file', required=True)
    args = parser.parse_args()

    pipeline = IPipeline(args.infile, args.chromsizes, args.genome, args.outfile)
    
    pipeline.check_dependencies()
    pipeline.check_files_exist()
    pipeline.check_input_formats()
    pipeline.check_if_sorted()
    # Wrap pipeline execution with Suppress to hide output
    with Suppress():
        pipeline.paired_check()
        pipeline.pair_split()
        pipeline.block_find()
        pipeline.msp1_split()
        pipeline.logging()
        pipeline.msp1_clip()
        pipeline.file_merge()


if __name__ == '__main__':
    main()

