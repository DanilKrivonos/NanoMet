import logging
import numpy as np
from tqdm import tqdm
from umap import UMAP
from subprocess import call
from src.Kmer_abundance import get_kmers_fereq
from scipy.spatial.distance import pdist, squareform
from Bio.SeqIO import parse
from hdbscan import HDBSCAN
from pandas import DataFrame, read_csv

a_logger = logging.getLogger()
a_logger.setLevel(logging.DEBUG)

def select_read(fastq, silva, out_dir):

    awk_command = "'{" + 'print $1"\t"$3' + "}'"
    call('minimap2 -ax map-ont {} {} |samtools view -F 260 |awk -F "{}" {} > {}/selected_reads.txt'.format(silva, fastq, '\t', awk_command, out_dir), shell=True)
    print('minimap2 -ax map-ont {} {} |samtools view -F 260 |awk -F "{}" {} > {}/selected_reads.txt'.format(silva, fastq, '\t', awk_command, out_dir))
    selected_reds = list(read_csv('{}/selected_reads.txt'.format(out_dir), sep='\t', header=None)[0].values)

    return selected_reds

def get_clustering(fastq, silva, out_dir):
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
    selected_reds = select_read(fastq, silva, out_dir)
    umap_model = UMAP(n_neighbors=200,
                    min_dist=0.4,
                    metric='manhattan')

    hdbscan_model = HDBSCAN(min_samples=1, 
                            cluster_selection_epsilon=0.1, 
                            metric='manhattan')

    READ_IDS = []
    DATASET = []
    open_fastq = parse(fastq, 'fastq')

    a_logger.debug('UMAP stage ...')

    for sequence in tqdm(open_fastq):
        if sequence.id not in selected_reds:
            continue
        READ_IDS.append(sequence.id)
        DATASET.append(list(get_kmers_fereq(sequence.seq, 5).values()))
    
    DATASET = np.array(DATASET)
    READ_IDS = np.array(READ_IDS)
    DISTANCE_MATRIX = squareform(pdist(DATASET, 'cosine'))  
    indexex_to_drop = []
    a_logger.debug('Denosing ...')

    for idx in tqdm(range(len(DISTANCE_MATRIX))):
        
        dist_vector = list(DISTANCE_MATRIX[idx])
        dist_vector.remove(0.0)
        
        if min(dist_vector) < 0.1:
            continue
        
        indexex_to_drop.append(idx)

    a_logger.debug('HDBSCAN clustering ...')

    DATASET = np.delete(DATASET, (indexex_to_drop), axis=0)
    READ_IDS = np.delete(READ_IDS, (indexex_to_drop), axis=0)
    RESULT = umap_model.fit_transform(DATASET)    
    hdbscan_pedicted = hdbscan_model.fit(RESULT)
    hdbscan_labels = hdbscan_pedicted.labels_
    RESULT_DICT = {'Read id' : [], 
                   '1 UMAP COMPONENT' : [], 
                   '2 UMAP COMPONENT' : [], 
                   'Class' : []}
    
    for idx in range(len(hdbscan_labels)):
        
        RESULT_DICT['Read id'].append(READ_IDS[idx])
        RESULT_DICT['1 UMAP COMPONENT'].append(RESULT[:, 0][idx])
        RESULT_DICT['2 UMAP COMPONENT'].append(RESULT[:, 1][idx])
        RESULT_DICT['Class'].append(hdbscan_labels[idx])

    RESULT_DF = DataFrame(RESULT_DICT)
    RESULT_DF = RESULT_DF[RESULT_DF['Class'] != -1]
    RESULT_DF.to_csv('./{}/Read_clusters.tsv'.format(out_dir), sep='\t')

    return RESULT_DF
