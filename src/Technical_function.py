import os
import logging
import numpy as np
from tqdm import tqdm
from subprocess import call
from Bio.SeqIO import parse

a_logger = logging.getLogger()
a_logger.setLevel(logging.DEBUG)

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
        
        call('porechop -i {} --verbosity 1 -t 100 --require_two_barcodes --extra_end_trim 0 -o {}/trimmed_barcode.fastq'.format(fastq, out_dir), shell=True)
        fastq = out_dir + "/trimmed_barcode.fastq"
        output = out_dir + "/trimmed_barcode.fastq"

    # cut primers
    if trimm_primer == True:

        opn_fastq = open(fastq)
        switch = 0
        # cut primers
        with open('{}/trimmed_primer.fastq'.format(out_dir), 'w') as trimmed_fasta:
            for line in opn_fastq:
                if switch == 0:

                    trimmed_fasta.write(line)
                    switch = 1
                    continue

                if switch == 1:

                    trimmed_fasta.write(line[len(hangF): -len(hangR)] + "\n")
                    switch = 0
                    continue
        
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
        call('filtlong --min_length 350 --max_length 600 --mean_q_weight 12 {} > {}'.format(fastq, output), shell=True)

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
        