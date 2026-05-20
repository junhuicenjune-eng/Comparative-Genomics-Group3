# 1. Transform proteasoms to directories
# 2. read ORFs from best_cluster list
# 3. Use ORFs as key to grab AA sequences
# 4. Write fasta files for each cluster; limited amount

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import random

random.seed(9) #for reproducibility of selected clusters

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract top BLAST hit sequence for each genome"
    )
    parser.add_argument("-i", "--input", required=True, help="Cluster list")
    parser.add_argument("-r", "--reference", required=True, help="refrence file e.g. proteom")
    parser.add_argument("-n", "--sample",type=int, required=True, help="Number of selected clusters")
    return parser.parse_args()

def main():
    args = parse_args()

    proteoms = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta")) #prepare a dictionary of all ORFs

    with open(args.input, "r") as infile:
        clusters = [line.strip() for line in infile] #read clusters a list of lines
        clusters = random.sample(clusters, args.sample) #select a random sample of the clusters

    for cluster in clusters:
        ORFs = cluster.strip().split("\t") #split cluster lines into ORFs
        cluster_file = ORFs[0] + "_cluster.fa" #used to name generated file after ref_ORF

        with open(cluster_file, "w") as outfile: 
            for ORF in ORFs: #iterate through selected ORFs
                outfile.write(">"+ORF+"\n") 
                outfile.write(str(proteoms[ORF].seq)+"\n") #get protein seqeunce from ORF dictionary


if __name__ == "__main__":
    main()
