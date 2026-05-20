#
# Read Balstp result, extract the best matched query ORF for each reference ORF
#


import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract top BLAST hit sequence for each genome"
    )
    parser.add_argument("-b", "--blast", required=True, help="BLAST output file, outfmt 6")
    parser.add_argument("-r", "--reference", required=True, help="refrence file e.g. proteom")
    parser.add_argument("-o", "--output", required=True, help="Output fasta file")
    parser.add_argument("-f", "--filter", required=True, help="Determines category from which the best hit is selected for each value")
    return parser.parse_args()


def main():
    args = parse_args()

    # Load all reference sequences into a dictionary
    reference = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    best_hits = {}

    with open(args.blast) as blast_file:
        for line in blast_file: #iterate through hits in blast result
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t")

            # Standard BLAST outfmt 6:
            # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            qseqid = fields[0] #query ID
            sseqid = fields[1] #reference ID
            sstart = int(fields[8]) #start postion in reference
            send = int(fields[9]) #end position in refence
            evalue = float(fields[10]) #evalue of hit

            #filter_cat = args.filter #choose category to filter for best result

            if args.filter == "qseqid":
                filter_cat = qseqid
            elif args.filter == "sseqid":
                filter_cat = sseqid
            else:
                raise Exception("Filter must be qseqid or sseqid")

            #uses values in filter category as key in dictionary; keeps onyl the hit with the lowest e-value per genome
            #Keep first hit if e-values are equal
            if filter_cat not in best_hits or evalue < best_hits[filter_cat]["evalue"]:
                best_hits[filter_cat] = {
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "sstart": sstart,
                    "send": send,
                    "evalue": evalue
                }

    with open(args.output, "w") as out:
        for filter_cat, hit in best_hits.items(): #iterate through optimal hits in each genome
            qseqid = hit["qseqid"]
            sseqid = hit["sseqid"]

            out.write(sseqid + "\t" + qseqid + "\n")


if __name__ == "__main__":
    main()