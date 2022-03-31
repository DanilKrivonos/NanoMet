import os
from os import listdir
from subprocess import call
from pandas import read_csv

def align_consensus(out_dir, path_to_SILVA):

    clusters = listdir('{}/Consensus/'.format(out_dir))
    os.mkdir('{}/OTU'.format(out_dir))

    for cluster in listdir(clusters):

        os.mkdir('{}/Consensus/{}'.format(out_dir, cluster))
        call('minimap2 -ax map-ont {}/Consensus/medaka_result/consensus.fasta {} > {}/OTU/{}/{}_otu.sam'.format(out_dir, path_to_SILVA, out_dir, cluster, cluster), shell=True) 
        call('samtools view -F 260 {}/OTU/{}_otu.sam| awk -F "{} "{}}print {} > {}/OTU/{}_filter_otu.txt'.format(out_dir, cluster, out_dir, '\t"', '{', '$3}', cluster), shell=True)

def get_otu(out_dir):

    clusters = listdir('{}/Consensus/'.format(out_dir))

    for cluster in clusters:

        m = read_csv('{}/OTU/{}_filter_otu.txt'.format(out_dir, cluster), sep='\t', header=None)
