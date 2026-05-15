import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract top BLAST hit sequence for each genome"
    )
    parser.add_argument("-b", "--blast", required=True, help="BLAST output file, outfmt 6")
    parser.add_argument("-g", "--genomes", required=True, help="Concatenated genomes fasta file")
    parser.add_argument("-o", "--output", required=True, help="Output fasta file")
    return parser.parse_args()


def main():
    args = parse_args()

    # Load all genome sequences into a dictionary
    genomes = SeqIO.to_dict(SeqIO.parse(args.genomes, "fasta"))

    best_hits = {}

    with open(args.blast) as blast_file:
        for line in blast_file:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t")

            # Standard BLAST outfmt 6:
            # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            qseqid = fields[0]
            sseqid = fields[1]
            sstart = int(fields[8])
            send = int(fields[9])
            evalue = float(fields[10])

            genome_id = sseqid

            # Keep first hit if e-values are equal
            if genome_id not in best_hits or evalue < best_hits[genome_id]["evalue"]:
                best_hits[genome_id] = {
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "sstart": sstart,
                    "send": send,
                    "evalue": evalue
                }

    with open(args.output, "w") as out:
        for genome_id, hit in best_hits.items():
            sseqid = hit["sseqid"]
            sstart = hit["sstart"]
            send = hit["send"]

            genome_seq = genomes[sseqid].seq

            if sstart <= send:
                subseq = genome_seq[sstart - 1:send]
                strand = "+"
            else:
                subseq = genome_seq[send - 1:sstart].reverse_complement()
                strand = "-"

            header = (
                f">{genome_id}_top_hit "
                f"coords={sstart}-{send} "
                f"strand={strand} "
                f"evalue={hit['evalue']}"
            )

            out.write(header + "\n")
            out.write(str(subseq) + "\n")


if __name__ == "__main__":
    main()