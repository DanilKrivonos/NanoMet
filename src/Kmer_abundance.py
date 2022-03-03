from itertools import product 

def normalize(kmers):
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
    norm = sum(list(kmers.values()))
    
    for kmer in kmers.keys():
        
        kmers[kmer] = kmers[kmer]/ norm
    
    return kmers


def get_kmers_fereq(seq, k=2):
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
    kmers = {"".join(kmer) : 0 for kmer in list(product("AGTC", repeat=k))}
    
    step = 1
    start = 0
    end = k
    cken = []
    while end != len(seq) - 1:
        #try:
        kmers[str(seq[start: end])] += 1
    #    except:
        
        start, end = start + step, end + step
        
    step = 1
    start = 0
    end = k

    while end != len(seq.reverse_complement()) - 1:
        
        #try:
        kmers[str(seq.reverse_complement()[start: end])] += 1
   #     except:
        
        start, end = start + step, end + step
        
    kmers = normalize(kmers)
    
    return kmers


