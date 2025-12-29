from chaining import chain_driver_np
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from collections import defaultdict
from itertools import combinations
from multiprocessing import shared_memory
import numpy as np
import time

'''
Performs pairwise co-linear chaining of genomic regions, allowing varied scores and starting from motif lists of
genomic regions.

Uses preprocessed data (list of motifs in each region) and the motif pair data 
to output the anchors for each pair across all regions

Note that this program parallelizes chaining but NOT finding anchors.
'''

##########################################################
# Settings

# Adjust score mean to this value
# 'None' for no adjustment
SCORE_CENTER = 5

###########################################################
# Global memory variables

global_anchor_array = None
global_anchor_shm = None
global_anchor_name = f"anchors_{int(time.time())}"
BATCH_SIZE = 100000


###########################################################
# Helpers for parallelization

def chunkify(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

#######################################################
# Motif dict creation
   
# Puts motif order into dictionary
# Takes in motif list file
def get_motif_loc_dict(data) :
    # Holds where each motif is located {MOTIF: {region: [loc, ..., loc] region:[loc, ..., loc]}}
    motif_loc_dict = defaultdict(lambda: defaultdict(list))

    print(f"Filling Motif Dict")

    # Open just the file, format: region \n MOTIF \n MOTIF ... region \n ...
    # Duplicates should have been removed already
    with open(data, 'r') as reg_motifs :
        curr_reg = ""
        curr_idx = 0
        for line in reg_motifs :
            if 'ACR: ' in line:
                curr_reg = line.rstrip().split('ACR: ')[1]
                curr_idx = 0

            else :
                motif_loc_dict[line.rstrip()][curr_reg].append(curr_idx)
                curr_idx += 1

    return motif_loc_dict

# Gets score dictionary (if custom score file was passed in)
# Adjusts score if SCORE_CENTER is not 'None'
def get_scores(CUSTOM_SCORE) :
    score_dict = None
    if CUSTOM_SCORE != None :
        # Score dict: {motif : score}
        score_dict = {}
        with open(CUSTOM_SCORE, 'r') as cs :
           next(cs)
           for line in cs :
               score_dict[line.split('\t')[0]] = float(line.rstrip().split('\t')[1])

        # Adjust the score
        if SCORE_CENTER != None :
            current_mean = sum(score_dict.values()) / len(score_dict)
            scale = SCORE_CENTER / current_mean
            score_dict = {k: v * scale for k, v in score_dict.items()}

    return score_dict

#####################################################################
# Calculates the total number anchors for all pairs of regions and then allocates space


# Calculate the total number of anchors (motif matches between pairs of regions)
# Returns an int: total number of anchors, and a dict: maps each region to its start and stop points in the anchor space
def calculate_no_anchors(motif_loc_dict) :
    print("Calculating Anchor Dict Allocation Space")

    # Running total of space needed
    total = 0

    # To track individual space needed
    # {(reg1, reg2) : space}
    reg_space_dict = defaultdict(int)

    # To give start and stop indices
    # {(reg1, reg2) : (start, stop)}
    reg_start_stop_dict = defaultdict(lambda: (0, 0))

    # single_motif_dict = {region: [loc, loc, ...], region: [loc, loc, ...], ...}
    for motif_name, single_motif_dict in tqdm(motif_loc_dict.items()) :
        # Number of motifs in each region
        sizes = {reg: len(single_motif_dict[reg]) for reg in list(single_motif_dict.keys())}

        # For every pair of regions, add needed space to reg_space_dict and total
        for reg1, reg2 in combinations(single_motif_dict.keys(), 2) :
            key = (reg1, reg2) if reg1 < reg2 else (reg2, reg1)
            pair_count = sizes[reg1] * sizes[reg2]
            reg_space_dict[key] += pair_count
            total += pair_count
    
    curr_start = 0
    for key, space in reg_space_dict.items() :
        reg_start_stop_dict[key] = (curr_start, curr_start + space)
        curr_start += space
    
    return total, reg_start_stop_dict


# Allocate space for anchors
def init_anchor_array(total_anchors, CUSTOM_SCORE):
    global global_anchor_array, global_anchor_shm
    if CUSTOM_SCORE == None :
        global_anchor_shm = shared_memory.SharedMemory(create=True, size=total_anchors * 2 * 4, name = global_anchor_name)  # 2 ints per anchor, 4 bytes each
        print(f"Allocating {total_anchors * 8} Bytes")
        global_anchor_array = np.ndarray((total_anchors, 2), dtype=np.int32, buffer=global_anchor_shm.buf)
    else :
        global_anchor_shm = shared_memory.SharedMemory(create=True, size=total_anchors * 3 * 4, name = global_anchor_name)  # 3 floats per anchor, 8 bytes each
        print(f"Allocating {total_anchors * 12} Bytes")
        global_anchor_array = np.ndarray((total_anchors, 3), dtype='float32', buffer=global_anchor_shm.buf)


def cleanup_anchor_array():
    global global_anchor_shm
    global_anchor_shm.close()
    global_anchor_shm.unlink()



######################################################################
# Anchors

# Find anchors given location dicts
def find_anchors(anchor_loc_dict, motif_loc_dict, score_dict, CUSTOM_SCORE) :
    print("Finding Anchors")

    # {(reg1, reg2) : curr_filled}, where curr_filled is the number of anchors currently
    # in global_anchor_array for the given pair
    curr_filled_per_pair = defaultdict(int)

    # single_motif_dict = {region: [loc, loc, ...], region: [loc, loc, ...], ...}
    for motif_name, single_motif_dict in tqdm(motif_loc_dict.items()) :

        for reg1, reg2 in combinations(single_motif_dict.keys(), 2):
            key = (reg1, reg2) if reg1 < reg2 else (reg2, reg1)

            # locations of motif in smaller (alphabetically) region
            anchors1 = single_motif_dict[reg1 if reg1 < reg2 else reg2]
            # locations of motif in larger (alphabetically) region
            anchors2 = single_motif_dict[reg2 if reg1 < reg2 else reg1]

            # Write to the correct slice of the global array
            n1, n2 = len(anchors1), len(anchors2)
            n_total = n1 * n2

            start_idx = anchor_loc_dict[key][0] + curr_filled_per_pair[key]
            
            # Add anchors to global_anchor_array with score if custom score is given
            if CUSTOM_SCORE == None :
                i = 0
                for a1 in anchors1:
                    for a2 in anchors2:
                        global_anchor_array[start_idx + i, 0] = a1
                        global_anchor_array[start_idx + i, 1] = a2
                        i += 1
            else :
                i = 0
                for a1 in anchors1:
                    for a2 in anchors2:
                        global_anchor_array[start_idx + i, 0] = a1
                        global_anchor_array[start_idx + i, 1] = a2
                        global_anchor_array[start_idx + i, 2] = score_dict[motif_name]
                        i += 1

            curr_filled_per_pair[key] += n_total

####################################################################
# Global chaining


# Parallelized chaining worker
# Returns a list where each item is f"{reg1}\t{reg2}\t{chain_score}\t{no_anchors}"
def batch_pair_chaining(args):
    pairs, total_anchors, custom_score, global_anchor_name = args
    shm = shared_memory.SharedMemory(name=global_anchor_name)

    # Create NumPy array wrapper (no copy)
    if custom_score :
        shared_array = np.ndarray((total_anchors, 3), dtype='float32', buffer=shm.buf)
    else :
        shared_array = np.ndarray((total_anchors, 2), dtype=np.int32, buffer=shm.buf)
    results = []
    for reg1, reg2, start, stop in pairs:  
        chain_len = chain_driver_np(shared_array[start:stop], custom_score)
        results.append(f"{reg1}\t{reg2}\t{chain_len}\t{stop - start}\n")

    return results

# Finds pairwise chain length (global)
# Uses global anchor dict dict from earlier and the output directory     
# Make sure to clear the output folder first since we are appending to files 
# Printed file in out_file in the format reg1\treg2\tchain_score\tno_anchors
def chain_all_pairs(anchor_loc_dict, total_anchors, out_file, CUSTOM_SCORE) :

    print("Chaining")

    # Parallelized chaining
    anchor_items = [(a, b, start, stop) for (a, b), (start, stop) in anchor_loc_dict.items()]

    batches = list(chunkify(anchor_items, BATCH_SIZE))
    parallelized_inputs =  [(batch, total_anchors, CUSTOM_SCORE != None, global_anchor_name) for batch in batches]

    with Pool(cpu_count()) as pool, open(out_file, 'w') as out:
        try :
            for lines in tqdm(pool.imap_unordered(batch_pair_chaining, parallelized_inputs), total=len(batches)):
                out.writelines(lines)
        finally:
            pool.close()
            pool.join()           


##############################################################
# Driver

def chain_global_driver(input_path, out_file, custom_score) :
    # Set to None if not weighted
    score_dict = get_scores(custom_score)
    motif_l_dict = get_motif_loc_dict(input_path)
    space, anchor_loc_dict = calculate_no_anchors(motif_l_dict)
    init_anchor_array(space, custom_score)
    find_anchors(anchor_loc_dict, motif_l_dict, score_dict, custom_score)
    chain_all_pairs(anchor_loc_dict, space, out_file, custom_score)
    cleanup_anchor_array()


if __name__ == "__main__":
    # Path to file containing the genomic region motif lists
    MOTIFS = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_random_v2_50-50/all_motifs.txt'

    # Set to None if not weighted
    CUSTOM_SCORE = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_upstream_50-50/motif_scoring.tsv'

    # Path to output file
    # Note we are appending to this file (so make sure it is empty or does not exist)
    OUTPUT = "/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/no_cluster_random_v2_50-50/all_chain.tsv"

    start_time = time.time()
    chain_global_driver(MOTIFS, OUTPUT, CUSTOM_SCORE)
    print(f"Finished in {time.time() - start_time} s", flush=True)




