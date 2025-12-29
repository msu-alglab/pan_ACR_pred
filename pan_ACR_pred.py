import argparse
import os.path
from os import mkdir
import shutil
from get_motifs_from_full import get_motifs_from_full
from get_motif_sequence import get_motif_sequence
from motif_relevance import motif_relevance
from create_motif_scoring import create_motif_scoring
from anchor_chain_driver import chain_global_driver
from chain_to_csv import chain_to_csv

if __name__ == "__main__":
    # Parse input --------------------------------------------------------

    parser = argparse.ArgumentParser(prog='pan_ACR_pred', description='Pangenomic Chromatin Accessibility Prediction')
    parser.add_argument("-o", "--output", help="Output Directory", required=True)
    parser.add_argument("-k", "--known", help="List of known ACRs formatted as ChrX_<start>to<end>", required=True)
    parser.add_argument("-u", "--unknown", help="List of unknown regions formatted as ChrX_<start>to<end>. _pos or _neg can also be appended to the" \
    " end of the region name for supervised machine learning.", required=True)
    parser.add_argument("-kf", "--known_fimo", help="A directory containing a fimo_out_*/ folder with a fimo.tsv file detailing all candidate CRE motifs in the genome from which" \
    " the known ACRs originate", required=True)
    parser.add_argument("-uf", "--unknown_fimo", help="A directory containing a fimo_out_*/ folder with a fimo.tsv file detailing all candidate CRE motifs in the genome from which" \
    " the unknown regions originate (if different from fimo_known)")
    parser.add_argument("-s", "--scoring", help="'None' for unweighted scoring, 'Default' for scoring weighted" \
    " by motif relevance, filepath for custom scoring")
    parser.add_argument("-r", "--reset", action="store_true", help="Clear all current temp files when this flag is present (rerun all processes).")

    args = parser.parse_args()

    # Validate input --------------------------------------------------------

    # Optional
    if args.scoring != None and args.scoring != "None" and args.scoring != "Default" and not os.path.isfile(args.scoring):
        raise FileNotFoundError(f"File not found: {args.scoring}")
    if args.unknown_fimo != None and args.fimo != "None" and not os.path.exists(args.unknown_fimo):
        raise FileNotFoundError(f"File not found: {args.unknown_fimo}")
    
    # Required
    if not os.path.isfile(args.known):
        raise FileNotFoundError(f"File not found: {args.known}")
    if not os.path.isfile(args.unknown):
        raise FileNotFoundError(f"File not found: {args.unknown}")
    if not os.path.exists(args.known_fimo):
        raise FileNotFoundError(f"File not found: {args.known_fimo}")
    if not os.path.exists(args.output):
        raise FileNotFoundError(f"Directory not found: {args.output}")
    
    # Set up -----------------------------------------------------------------

    if args.reset:
        try:
            shutil.rmtree("./.temp")
        except FileNotFoundError:
            pass
    
    try:
        mkdir("./.temp")
    except FileExistsError:
        pass


    if args.scoring == None:
        args.scoring = "Default"
    if args.unknown_fimo == None:
        args.unknown_fimo = "None"

    # Extract motif sequence for known ACRs ---------------------------------
    if not os.path.isfile("./.temp/acr_motifs.txt"):
        get_motifs_from_full(args.known, args.known_fimo, "./.temp/acr_motifs.txt")

    # Scoring ----------------------------------------------------------------

    if args.scoring == "Default":
        try:
            mkdir("./.temp/motif_seq_full_known")
            get_motif_sequence(args.known_fimo, "./.temp/motif_seq_full_known")
        except FileExistsError:
            pass
        if not os.path.isfile("./.temp/motif_relevance.txt"):
            motif_relevance("./.temp/motif_seq_full_known", "./.temp/acr_motifs.txt", "./.temp/motif_relevance.txt")
        if not os.path.isfile("./.temp/motif_scoring.txt"):
            create_motif_scoring("./.temp/motif_relevance.txt", "./.temp/motif_scoring.txt")

    # Extract motif sequence for unknown ACRs ---------------------------------
    
    if not os.path.isfile("./.temp/unknown_motifs.txt"):
        if args.unknown_fimo == "None":
            get_motifs_from_full(args.unknown, args.known_fimo, "./.temp/unknown_motifs.txt")
        else:
            get_motifs_from_full(args.unknown, args.unknown_fimo, "./.temp/unknown_motifs.txt")

    # Create a single motif list -----------------------------------------------

    if not os.path.isfile("./.temp/all_motifs.txt"):
        with open("./.temp/all_motifs.txt", 'w') as destination_file:
            with open("./.temp/acr_motifs.txt", 'r') as file1:
                destination_file.write(file1.read())
            with open("./.temp/unknown_motifs.txt", 'r') as file2:
                destination_file.write(file2.read())
    
    # Co-linear chaining ---------------------------------------------------------
    if not os.path.isfile("./.temp/chain_all.tsv"):
        if args.scoring == "None":
            chain_global_driver("./.temp/all_motifs.txt", "./.temp/chain_all.tsv", None)
        elif args.scoring == "Default":
            chain_global_driver("./.temp/all_motifs.txt", "./.temp/chain_all.tsv", "./.temp/motif_scoring.txt")
        else:
            chain_global_driver("./.temp/all_motifs.txt", "./.temp/chain_all.tsv", args.scoring)

    # Output to csv file ----------------------------------------------------------

    if args.output[-1] == "/":
        outpath = f"{args.output}"
    else:
        outpath = f"{args.output}/"
    chain_to_csv("./.temp/chain_all.tsv", args.unknown, args.known, f"{outpath}acr_pred.csv")

    

