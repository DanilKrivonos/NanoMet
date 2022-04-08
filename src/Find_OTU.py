import os
import numpy as np
from os import listdir
from subprocess import call
from pandas import read_csv, DataFrame
from Bio.SeqIO import parse

def align_consensus(out_dir, path_to_SILVA):

    clusters = listdir('{}/Consensus/'.format(out_dir))
    os.mkdir('{}/OTU'.format(out_dir))

    for cluster in clusters:

        os.mkdir('{}/OTU/{}'.format(out_dir, cluster))
        call('minimap2 -ax map-ont {} {}/Consensus/{}/medaka_result/consensus.fasta > {}/OTU/{}/{}_otu.sam'.format(path_to_SILVA, out_dir, cluster, out_dir, cluster, cluster), shell=True) 
        awk_print = "'{" + 'print $1"\t"$3' + "}'"
        call('samtools view -F 260 {}/OTU/{}/{}_otu.sam| awk -F "{}" {} > {}/OTU/{}/{}_filter_otu.txt'.format(out_dir, cluster, cluster, '\t', awk_print, out_dir, cluster, cluster), shell=True)

def get_cluster_presentstive(out_dir, cluster, abs_presentative):
    
    ani_table = read_csv('{}/proOTU/{}/{}.txt'.format(out_dir, cluster, cluster), sep='\t', header=None)

    for sample in set(ani_table[0]):

        subset = ani_table[ani_table[0] == sample]
        max_val = np.max(subset[2].values)
        abs_presentative[subset[subset[2] == max_val][1].values[0].split('/')[-1][: -6]] += 1
        
    return abs_presentative

def get_taxonomy(abs_presentative, taxonomy):
    
    abs_presentative = DataFrame({'Taxonomy' : list(abs_presentative.keys()), 'Presentative' : list(abs_presentative.values())})
    
    for tax_id in abs_presentative.copy().index:

        ID = abs_presentative['Taxonomy'][tax_id]
        abs_presentative['Taxonomy'][tax_id] =  taxonomy[ID]

    return abs_presentative

def get_presentative(out_dir, path_to_SILVA, threads):

    clusters = listdir('{}/Consensus/'.format(out_dir))
    ALL_POTENTIAL_OTU = {} 
    output = '{}/potential_otus'.format(out_dir)
    os.mkdir(output)
    silva_db = parse(path_to_SILVA, 'fasta') 
    all_otus  = []

    for cluster in clusters:

        potential_OTU = list(read_csv('{}/OTU/{}/{}_filter_otu.txt'.format(out_dir, cluster, cluster), sep='\t', header=None)[1].values)
        ALL_POTENTIAL_OTU[cluster] = potential_OTU
        all_otus.extend(potential_OTU)
        
    all_otus = list(set(all_otus))
    abs_presentative = {i: 0 for i in all_otus}
    taxonomy = {i: '' for i in all_otus}

    for record in silva_db:
        
        if record.id in all_otus:

            taxonomy[record.id] = record.description.split()[1].split(';')[-2]
            
            with open('{}/potential_otus/{}.fasta'.format(out_dir, record.id), 'w') as opne_fasta:
                
                opne_fasta.write('>{}\n{}'.format(record.id, record.seq))

    os.mkdir('{}/proOTU/'.format(out_dir))

    for cluster in clusters:

        os.mkdir('{}/proOTU/{}'.format(out_dir, cluster))

        for pro_otu in ALL_POTENTIAL_OTU[cluster]:
        #    print('fastANI --fragLen 200 --ql {}/FastANI_result/{}/path.txt  -r {}/potential_otus/{}.fasta -o {}/proOTU/{}/{}.txt'.format(out_dir, cluster, out_dir, pro_otu, out_dir, cluster, cluster))
            call('fastANI -t {} --fragLen 350 --ql {}/FastANI_result/{}/path.txt  -r {}/potential_otus/{}.fasta -o {}/proOTU/{}/{}.txt'.format(threads, out_dir, cluster, out_dir, pro_otu, out_dir, cluster, cluster), shell=True)
            abs_presentative = get_cluster_presentstive(out_dir, cluster, abs_presentative)

    abs_presentative = get_taxonomy(abs_presentative, taxonomy)

    return abs_presentative