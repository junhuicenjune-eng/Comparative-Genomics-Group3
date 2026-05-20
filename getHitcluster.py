#
# Scan combined basepair list; use ref_ORF as key, store query_ORFs as list; check list length
#

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract top BLAST hit sequence for each genome"
    )
    parser.add_argument("-i", "--input", required=True, help="Combined list of paired ref- and query ORFs")
    parser.add_argument("-o", "--output", required=True, help="Output list of hit clusters")
    return parser.parse_args()


def main():
    args = parse_args()

    clusters = {}

    with open(args.input, "r") as infile:
        for line in infile: #iterate through list of paired ref_ORF - query_ORF
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t") #collects ref_ORF amd query_ORF

            ref_ORF = fields[0]
            query_ORF = fields[1]

            #Collect clusters in a dictionary usign the ref_ORF as key
            if ref_ORF not in clusters: #if ref_ORF is not included initialise new entry
                clusters[ref_ORF] = {
                    "ref_ORF":ref_ORF,
                    "query_ORFs": [query_ORF]
                }
            else: #if ref_ORF is included extend the list of matched query_ORFs
                clusters[ref_ORF]["query_ORFs"].append(query_ORF)

    with open(args.output, "w") as outfile:
        m = 0 #total ref_ORFs
        n = 0 #complete clusters

        for ref_ORFs, cluster in clusters.items():
            m += 1
            ref_ORF = cluster["ref_ORF"]
            query_ORFs = cluster["query_ORFs"]

            if len(query_ORFs) == 3:
                outfile.write(ref_ORF + "\t" + query_ORFs[0] + "\t" + query_ORFs[1] + "\t" + query_ORFs[2] + "\n")
                n += 1

    print("Compiled %d complete clusters out of %d ref_ORFs"%(n,m))


if __name__ == "__main__":
    main()