from asyncio import FastChildWatcher
import os 
import sys
import shutil
import logging
import argparse 
from pickle import dump
from subprocess import call
from pandas import read_csv
from os import listdir 
from tqdm import tqdm
from src.Technical_function import cut_tails, get_read_filter, get_cluster_sequences
from src.Clustering_stage import get_clustering
from src.Get_consensus import get_closer_seq, get_consensus
from src.Find_OTU import align_consensus, get_presentative

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

    group1.add_argument('-path_to_SILVA',
                        type=str,
                        help='path to fastq',
                        default='./data/SILVA_138.1_SSURef_NR99_tax_silva.fasta')

    group2.add_argument('--trimm_adapter', 
                        help='path to fastq',
                        default=False,
                        action='store_true')

    group2.add_argument('--trimm_primer', 
                        help='',
                        default=False,
                        action='store_true')

    group2.add_argument('--disable_filter', 
                        help='',
                        default=False,
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
    trimm = args.trimm
    silva_db = args.path_to_SILVA
    hangF = args.primer_hangF
    hangR = args.primer_hangR
    disable_filter = args.disable_filter
    trimm_adapter = args.trimm_adapter
    trimm_primer = args.trimm_primer
    out_dir = args.out

    if out_dir is None:
        try:
            out_dir = './NanoMet_output'

            if os.path.exists(out_dir):

                shutil.rmtree(out_dir)
            os.makedirs(out_dir)

        except:
            print('Dorectory is exist!')

    a_logger = logging.getLogger()
    a_logger.setLevel(logging.DEBUG)

    output_file_handler = logging.FileHandler("{}/NanoMet.log".format(out_dir))
    stdout_handler = logging.StreamHandler(sys.stdout)

    a_logger.addHandler(output_file_handler)
    a_logger.addHandler(stdout_handler)

    if fastq is None and barcode_directory is None:

        a_logger.debug("Give an fastq/fastq.gz or directiry with one barcodes' fastq/fastq.gz")
        sys.exit()
    
    if not fastq is None:
        if '.gz' in fastq:

            call('gzip -dk {}/{}'.format(barcode_directory, fastq), shell=True)
            barcode_name = fastq.split('.fastq')[0]
    
    else:
        
        barcode_name = barcode_directory.split('/')
        barcode_name.remove('')
        barcode_name = barcode_name[-1]
        
        if '.gz' in listdir(barcode_directory)[0]:
            
            remove_before = []
            
            for fastq in listdir(barcode_directory):

                call('gzip -dk {}/{}'.format(barcode_directory, fastq), shell=True)
                remove_before.append('{}/{}'.format(barcode_directory, fastq.split('.gz')[0]))

            call('cat {} > {}/{}.fastq'.format(" ".join(remove_before), out_dir, barcode_name.split('/')[-1]), shell=True)
            [call('rm {}'.format(i), shell=True) for i in remove_before]
            fastq = '{}/{}'.format(out_dir,  barcode_name.split('/')[-1] + '.fastq')

        else:
            
            call('cat {} > {}/{}'.format(barcode_directory, out_dir, barcode_name.split('/')[-1]), shell=True)

    fastq = cut_tails(fastq, out_dir, trimm_adapter, trimm_primer, hangF, hangR)
    fastq = get_read_filter(fastq, disable_filter, out_dir)
    RESULT_DF = get_clustering(fastq, out_dir)
    cluster_out = get_cluster_sequences(RESULT_DF, fastq, out_dir)
    get_closer_seq(cluster_out, out_dir)
    get_consensus(cluster_out, out_dir)
    align_consensus(out_dir, silva_db)
    abs_presentative = get_presentative(out_dir, silva_db)
    abs_presentative.to_csv('{}/Result_absolute_abundance.txt'.format(out_dir), sep='\t', index=False)
    #if trimm != False:

     #   os.mkdir('{}/TRIMMED_{}/'.format(out_dir, barcode_name))
      #  barcode_directory = cut_tails(barcode_directory, out_dir, hangF, hangR)


    #if disable_filter != False:

     #   os.mkdir('{}/FILTER_{}/'.format(out_dir, barcode_name))
    
    print('nononionionio')
    a_logger.debug('NanoMet processing is done!')
if __name__ == "__main__":
    main()