
from collections import defaultdict
from pathlib import Path
import pandas as pd

# Set to 'None' to use all ACRs
# Otherwise, specify ACRs to use (file with one ACR per line, formatted ChrX_<start>to<end>)
ACR_SET = None

# Specify which chromosome files to use
CHRS_TO_USE = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']

###########################################################################################

def dict_with_keys_zero():
    return {key: 0 for key in CHRS_TO_USE}

def z_score(x1, n1, x2, n2):
    p1 = x1 / n1
    p2 = x2 / n2
    p = (x1 + x2) / (n1 + n2)
    numerator = p1 - p2
    denominator = (p * (1 - p) * (1/n1 + 1/n2)) ** 0.5
    return numerator / denominator if denominator != 0 else 0

# Driver
def motif_relevance(FULL_SEQUENCE_FOLDER, ACR_FILE, OUT_FILE):
    if ACR_SET != None :
        acrs_to_use = set()
        with open(ACR_SET, 'r') as f :
            for line in f :
                acrs_to_use.add(line.rstrip())


    # {motif: {chr1 : num, chr2 : num...}}
    full_seq_chr_counts = defaultdict(dict_with_keys_zero)
    # {motif: num}
    full_seq_total_counts = defaultdict(int)

    total_full = 0

    for file_path in Path(FULL_SEQUENCE_FOLDER).iterdir():
        chromosome = str(file_path).rstrip().split('/')[-1]
        with open(file_path, 'r') as motif_file :
            for line in motif_file :
                if chromosome in CHRS_TO_USE :
                    if 'BREAK' not in line.split('\t')[0] :

                        full_seq_chr_counts[line.split('\t')[0]][chromosome] += 1
                        full_seq_total_counts[line.split('\t')[0]] += 1

                        total_full += 1

    # {motif: {chr1 : num, chr2 : num...}}
    acr_chr_counts = defaultdict(dict_with_keys_zero)
    # {motif: num}
    acr_total_counts = defaultdict(int)

    total_acr = 0

    with open(ACR_FILE, 'r') as acr_file : 

        curr_acr = ""
        for line in acr_file :
            # End of one
            if 'ACR: ' in line:
                # RESET
                curr_acr = line.rstrip().split('ACR: ')[1]

            elif ACR_SET != None :
                if curr_acr in acrs_to_use :
                    if curr_acr.split('_')[0] in CHRS_TO_USE :
                        acr_chr_counts[line.rstrip()][curr_acr.split('_')[0]] += 1
                        acr_total_counts[line.rstrip()] += 1
                        total_acr += 1
            else :
                if curr_acr.split('_')[0] in CHRS_TO_USE :

                    acr_chr_counts[line.rstrip()][curr_acr.split('_')[0]] += 1
                    acr_total_counts[line.rstrip()] += 1

                    total_acr += 1

    ################################################################################
    # Output

    with open(OUT_FILE, 'w') as out :
        out.write('Motif\t')
        out.write('Total_ACR\tTotal_Full\tACR_percent\tFull_percent\tZScore\n')


        sorted_acrs = sorted(
            full_seq_chr_counts.keys(),
            key = lambda acr: z_score(acr_total_counts[acr], total_acr, full_seq_total_counts[acr], total_full),
            reverse=True  # Optional: sort descending by difference
        )
        for acr in sorted_acrs :
            out.write(f"{acr}\t")
            out.write(f"{acr_total_counts[acr]}\t{full_seq_total_counts[acr]}\t")
            out.write(f"{acr_total_counts[acr]/total_acr}\t{full_seq_total_counts[acr]/total_full}\t")
            out.write(f"{z_score(acr_total_counts[acr], total_acr, full_seq_total_counts[acr], total_full)}\n")

if __name__ == "__main__":

    # Folder with files from get_motif_sequence.py 
    FULL_SEQUENCE_FOLDER = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/all_acrs_ml_exp1/fimo_seq_full_DAPv1_clustered'
    
    # List of motifs in all known ACRs
    ACR_FILE = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/all_acrs_ml_exp1/ref_set_motifs.txt'

    OUT_FILE = '/home/mwarr/Data/arabidopsis_one_genome/all_ACRs/all_acrs_ml_exp1/motif_relevance_ref_set.txt'

    motif_relevance(FULL_SEQUENCE_FOLDER, ACR_FILE, OUT_FILE)