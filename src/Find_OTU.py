import os
from os import listdir
from subprocess import call
from pandas import read_csv
from Bio.SeqIO import parse

def align_consensus(out_dir, path_to_SILVA):

    clusters = listdir('{}/Consensus/'.format(out_dir))
    os.mkdir('{}/OTU'.format(out_dir))

    for cluster in clusters:

        os.mkdir('{}/OTU/{}'.format(out_dir, cluster))
        call('minimap2 -ax map-ont {} {}/Consensus/{}/medaka_result/consensus.fasta > {}/OTU/{}/{}_otu.sam'.format(path_to_SILVA, out_dir, cluster, out_dir, cluster, cluster), shell=True) 
        awk_print = "'{" + 'print $1"\t"$3' + "}'"
        call('samtools view -F 260 {}/OTU/{}/{}_otu.sam| awk -F "{}" {} > {}/OTU/{}/{}_filter_otu.txt'.format(out_dir, cluster, cluster, '\t', awk_print, out_dir, cluster, cluster), shell=True)
