import os
import logging
import shutil
import numpy as np
from tqdm import tqdm
from subprocess import call
from Bio.SeqIO import parse

a_logger = logging.getLogger()
a_logger.setLevel(logging.DEBUG)

def create_dir(out_dir):
    if out_dir is None:
        
        out_dir = './NanoMet_output'
        
    if os.path.exists(out_dir):

        shutil.rmtree(out_dir)
    
    os.mkdir(out_dir)
    
    return out_dir

def cut_tails(fastq, out_dir, trimm_adapter, trimm_primer, hangF, hangR):
    """
    The functuion ...
    Parameters
    ----------
    reads : str
        path to ...
    out_dir : str
        path to ...
    hang1 : str
        Sequence ...
    hang2 : str
        Sequence ...
    Returns
    -------
    -
    """
    output = fastq
    # cut barcodes
    if trimm_adapter == True:

        call('porechop -i {} --verbosity 0 -t 100 --require_two_barcodes --extra_end_trim 0 -o {}/trimmed_barcode.fastq'.format(fastq, out_dir), shell=True)
        fastq = out_dir + "/trimmed_barcode.fastq"
        output = out_dir + "/trimmed_barcode.fastq"

    # cut primers
    if trimm_primer == True:

        opn_fastq = parse(fastq, 'fastq')

        # cut primers
        with open('{}/trimmed_primer.fastq'.format(out_dir), 'w') as trimmed_fasta:
            for record in opn_fastq:
                for idx in range(4):
                    if idx != 1 or idx != 3:

                        trimmed_fasta.write(record.format('fastq').split('\n')[idx] + '\n')

                    else:

                        trimmed_fasta.write(record.format('fastq').split('\n')[idx][len(hangF): -len(hangR)] + '\n')

        output = '{}/trimmed_primer.fastq'.format(out_dir)

    return output

def get_read_filter(fastq, disable_filter, out_dir):
    """
    The functuion ...
    Parameters
    ----------
    reads : str
        path to ...
    out_dir : str
        path to ...
    hang1 : str
        Sequence ...
    hang2 : str
        Sequence ...
    Returns
    -------
    -
    """
    if disable_filter != False:

        return fastq

    else:

        output = '{}/filter_reads.fastq'.format(out_dir)
        call('filtlong --min_length 350 --max_length 600 {} > {}'.format(fastq, output), shell=True)
       # call('filtlong --min_length 350 --max_length 600 --mean_q_weight 12 {} > {}'.format(fastq, output), shell=True)

        return output

def get_cluster_sequences(RESULT_DF, fastq, out_dir):
    """
    The functuion ...
    Parameters
    ----------
    reads : str
        path to ...
    out_dir : str
        path to ...
    hang1 : str
        Sequence ...
    hang2 : str
        Sequence ...
    Returns
    -------
    -
    """
    classes = list(set(RESULT_DF['Class'].values))
    read_ids = list(RESULT_DF['Read id'].values)
    os.mkdir('{}/Clusters'.format(out_dir))

    for cs in classes:
        
        os.mkdir('{}/Clusters/{}'.format(out_dir, cs))

    opn_fastq = parse(fastq, 'fastq')

    for sequence in opn_fastq:
        try:

            cs = RESULT_DF[RESULT_DF['Read id'] == sequence.id]['Class'].values[0]

            with open('{}/Clusters/{}/{}.fasta'.format(out_dir, cs, sequence.id), 'w') as open_read_fasta:
            
                open_read_fasta.write('>{}\n{}\n'.format(sequence.id, sequence.seq))

        except:
            pass
    
    cluster_out = '{}/Clusters/'.format(out_dir)

    return cluster_out
#########################################################################################################
def get_SILVA_taxon():
    """
    The functuion ...
    Parameters
    ----------
    reads : str
        path to ...
    out_dir : str
        path to ...
    hang1 : str
        Sequence ...
    hang2 : str
        Sequence ...
    Returns
    -------
    -
    """
    silva_fasta = parse('../data/SILVA_DB.fasta', 'fasta')
    full_tax = {}
    logging.debug('Getting SILVA taxons ...')

    for line in tqdm(silva_fasta):
        
        ID = line.description.split()[0]
        taxonomy = ' '.join(line.description.split()[1: ])
        full_tax[ID] = taxonomy     
################################################################################################################