import os
import numpy as np
from os import listdir
from tqdm import tqdm
from subprocess import call
from pandas import read_csv

def get_closer_seq(cluster_out, out_dir):

    output = '{}/FastANI_result/'.format(out_dir)
    os.mkdir(output)

    for cluster in tqdm(listdir(cluster_out)):
        
        os.mkdir(output + cluster)

        with open('{}/path.txt'.format(output + cluster), 'w') as path:
            for fasta in listdir('{}/{}'.format(cluster_out, cluster)):

                path.write('{}/Clusters/{}/{}\n'.format(out_dir, cluster, fasta))

        call('fastANI --fragLen 200 --ql {}/path.txt --rl {}/path.txt -o {}/fastANI_result.txt'.format(output + cluster, output + cluster, output + cluster), shell=True)
    


def get_consensus(cluster_out, out_dir):
    print('biobnonoin')
    os.mkdir('{}/Consensus'.format(out_dir))
    print('uibuiguiguiburm -r')
    for cluster in tqdm(listdir(cluster_out)):

        ani = read_csv('{}/FastANI_result/{}/fastANI_result.txt'.format(out_dir, cluster), sep='\t', header=None)
        consensus_prototipe = ''
        consensus_prototipe_value = 0
        polish_data = []

        for sample in set(ani[0].values):

            polish_data.append('{}/Clusters/{}/{}'.format(out_dir, cluster, sample.split('/')[-1]))
            mean_val = np.mean(ani[ani[1] == sample][2].values)
            
            if consensus_prototipe_value < mean_val:
                
                consensus_prototipe = sample.split('/')[-1]
                consensus_prototipe_value = mean_val

        polish_data.remove('{}/Clusters/{}/{}'.format(out_dir, cluster, consensus_prototipe))
        os.mkdir('{}/Consensus/{}'.format(out_dir, cluster))
        call('cp {out}/Clusters/{clstr}/{cons_prot} {out}/Consensus/{clstr}/consensus.fasta'.format(out=out_dir, 
                                                                                                    clstr=cluster, 
                                                                                                    cons_prot=consensus_prototipe), 
                                                                                                    shell=True)
        call('cat {} > {}/Consensus/{}/polishing_data.fasta'.format(' '.join(polish_data), out_dir, cluster), shell=True)

        call('minimap2 -ax map-ont {out}/Consensus/{clstr}/consensus.fasta {out}/Consensus/{clstr}/polishing_data.fasta > {out}/Consensus/{clstr}/polishing_data.sam'.format(out=out_dir, 
                                                                                                                                                                             clstr=cluster), 
                                                                                                                                                                             shell=True)
        
        
        call('racon {out}/Consensus/{clstr}/polishing_data.fasta {out}/Consensus/{clstr}/polishing_data.sam {out}/Consensus/{clstr}/consensus.fasta > {out}/Consensus/{clstr}/consensus_polished_racon.fasta'.format(out=out_dir,
                                                                                                                                                                                                                     clstr=cluster), 
                                                                                                                                                                                                                     shell=True)

        call('medaka_consensus -i {out}/Consensus/{clstr}/polishing_data.fasta -d {out}/Consensus/{clstr}/consensus_polished_racon.fasta -o {out}/Consensus/{clstr}/medaka_result'.format(out=out_dir,
                                                                                                                                                                                          clstr=cluster), 
                                                                                                                                                                                          shell=True)
 