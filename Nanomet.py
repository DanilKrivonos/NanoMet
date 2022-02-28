from asyncio import FastChildWatcher
import os 
import sys
import logging
import argparse 
from Bio.SeqIO import parse
from subprocess import call
from pandas import read_csv
from os import listdir 
from tqdm import tqdm
from Technical_function import cut_tails

def main():
    parser = argparse.ArgumentParser(description='BioCAT is a tool, which estimates the' + 
                                    'likelihood that a given orgnism is capable of producing of a given NRP')
    group1 = parser.add_argument_group("Main")
    group2 = parser.add_argument_group("Tech nical")
    
    group1.add_argument('-barcode_directory',
                        type=str,
                        help='barcode dir',
                        default=None)

    group1.add_argument('-fastq',
                        type=str,
                        help='path to fastq/fastq.gz',
                        default=None)
                        
    group1.add_argument('-barcode_seq',
                        type=str,
                        help='path to fastq',
                        default=None)

    group1.add_argument('--trimm', 
                        help='path to fastq',
                        default=None,
                        action='store_true')

    group1.add_argument('-primer_hangF',
                        type=str,
                        help='path to fastq/fastq.gz',
                        default=None)
                        
    group1.add_argument('-primer_hangR',
                        type=str,
                        help='path to fastq',
                        default=None)

    group2.add_argument('--disable_filter', 
                        help='path to fastq',
                        default=None,
                        action='store_true')

    group2.add_argument('-out',
                        type=str,
                        help='output path',
                        default=None)

    if len(sys.argv)==1:

        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    barcode_directory = args.barcode_directory
    fastq = args.fastq
    barcode_seq = args.barcode_seq
    trimm = args.trimm
    hangF = args.primer_hangF
    hangR = args.primer_hangR
    disable_filter = args.disable_filter
    out_dir = args.out

    if fastq is None and barcode_directory is None:

        print("Give an fastq/fastq.gz or directiry with one barcodes' fastq/fastq.gz")
        sys.exit()

    if not fastq is None:
        if 'gz' in fastq:

            call('gzip -dk {}'.format(fastq), shell=True)
            barcode_name = fastq.split('.fastq')[0]

    else:
        if 'gz' in listdir(barcode_directory)[0]:
            
            [call('gzip -dk {}'.format(fastq), shell=True) for fastq in listdir(barcode_directory)]
            barcode_name = barcode_directory.split('/')[-1]

    if trimm != False:

        os.mkdir('{}/TRIMMED_{}/'.format(out_dir, barcode_name))
        barcode_directory = cut_tails(barcode_directory, out_dir, hangF, hangR)

    if disable_filter != False:

        os.mkdir('{}/FILTER_{}/'.format(out_dir, barcode_name))
    




if __name__ == "__main__":
    main()