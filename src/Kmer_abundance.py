import random 
from os import listdir
from tqdm import tqdm

DATASET = []
TAXOMONY = []
read_ids = []
COLOR_DICT = {}
COLOR_DICT['UNMAPPED ON SILVIA DB'] = 'black'

def get_color(tax_list):
    
    color = ''
    
    while color not in tax_list.values() and color == '':
        
        color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
        
        if color == '#000000':
            
            color = ''
            
    return color

#for barcode in listdir('./DATA/FILTER_BARCODES/'):
    
c = 0
#    sample = parse('./DATA/FILTER_BARCODES/' + barcode, 'fastq')

for sequene in tqdm(sample):


    read_ids.append(sequene.id)

    try:
       # if _16s_card[_16s_card[0] == sequene.id][1].values[0].split(';')[-2] in top_10:

        TAXOMONY.append(_16s_card[_16s_card[0] == sequene.id][1].values[0].split(';')[-2])#[-1].split()[1])

        if _16s_card[_16s_card[0] == sequene.id][1].values[0].split(';')[-2] not in COLOR_DICT:

            COLOR_DICT[_16s_card[_16s_card[0] == sequene.id][1].values[0].split(';')[-2]] = get_color(COLOR_DICT)
        #else:

         #   TAXOMONY.append('UNMAPPED ON SILVIA DB')

    except IndexError:

        TAXOMONY.append('UNMAPPED ON SILVIA DB')

    #DATASET.append(np.concatenate((np.array(list(get_kmers_fereq(sequene.seq, 2).values())), np.array(list(get_kmers_fereq(sequene.seq.reverse_complement(), 2).values())))))
    DATASET.append(list(get_kmers_fereq(sequene.seq, 5).values()))
c += 1
        #for i in range(2, 7):
print()
print('Done!')