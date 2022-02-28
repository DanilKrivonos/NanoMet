import os
from subprocess import call
from Bio.SeqIO import parse
from os import listdir
from tqdm import tqdm
import numpy as np
import logging

def cut_tails(barcode_directory, out_dir, hang1, hang2):
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
    
    out_trimmed = os.mkdir(out_dir + '/TRIMMED_BARCODES/')

    for fastq in tqdm(barcode_directory):
        if 'gz' in fastq:
            continue
        # cut barcodes
        call('porechop -i {}/{} --verbosity 2 -t 100 --require_two_barcodes --extra_end_trim 0 -o {}/{}'.format(barcode_directory, fastq, out_dir, fastq), shell=True)
        opn_fastq = open('{}/{}'.format(barcode_directory, fastq))
        switch = 0
        # cut primers
        with open('{}/{}'.format(out_trimmed, fastq), 'w') as trimmed_fasta:
            for line in opn_fastq:
                if switch == 0:

                    trimmed_fasta.write(line)
                    switch = 1
                    continue

                if switch == 1:

                    trimmed_fasta.write(line[len(hang1): -len(hang2)] + "\n")
                    switch = 0
                    continue

    call('cat {}/* > {}/CONCATENATE_BARCODE/TRIMMED_{}.fastq'.format(out_trimmed, out_dir, barcode_directory.split('/')[-1]))

    return '{}/CONCATENATE_BARCODE/TRIMMED_{}.fastq'.format(out_dir, barcode_directory.split('/')[-1])

def get_read_filter(out_dir):


    lengths = []
    fastq_opn = parse(, 'fastq')
    try:
        os.mkdir('{}/'.format(out_dir))
    for line in fastq_opn:
        
        lengths.append(len(line.seq))
    
    sdand_div = np.std(lengths)
    plus_2std = int(443  + 2*sdand_div)
    minus_2std = int(443  - 2*sdand_div)
    
    if minus_2std > 0:
        
        call('filtlong --min_length {} --max_length {} --mean_q_weight 12 {}/{} > {}/{}'.format(minus_2std,
                                                                                                plus_2std,
                                                                                                out_dir,
                                                                                                barcode,
                                                                                                out_dir,
                                                                                                barcode), shell=True)
    else:
        
        call('filtlong --min_length 400 --max_length {} --mean_q_weight 12 {}/{} > {}/{}'.format(plus_2std,
                                                                                                out_dir,
                                                                                                barcode,
                                                                                                out_dir,
                                                                                                barcode), shell=True)
    



def get_SILVA_taxon():


    

    silva_fasta = parse('../data/SILVA_DB.fasta', 'fasta')
    full_tax = {}
    logging.debug('Getting SILVA taxons ...')

    for line in tqdm(silva_fasta):
        
        ID = line.description.split()[0]
        taxonomy = ' '.join(line.description.split()[1: ])
        full_tax[ID] = taxonomy     
        