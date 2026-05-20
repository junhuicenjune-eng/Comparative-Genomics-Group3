import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Input FASTA file with 4 proteomes")
    parser.add_argument("--output_prefix", default="Metaproteome", help="Prefix for output files")
    args = parser.parse_args()

    # Store sequences by proteome, keeping order
    proteome_names = ["9.fa", "11.fa", "22.fa", "27.fa"]

    proteomes = {
        proteome_names[0]: [],
        proteome_names[1]: [],
        proteome_names[2]: [],
        proteome_names[3]: [],
    }

    # Read and route each record
    for record in SeqIO.parse(args.fasta, "fasta"):
        header = record.id

        if header.startswith(proteome_names[0]):
            proteomes[proteome_names[0]].append(record)
        elif header.startswith(proteome_names[1]):
            proteomes[proteome_names[1]].append(record)
        elif header.startswith(proteome_names[2]):
            proteomes[proteome_names[2]].append(record)
        elif header.startswith(proteome_names[3]):
            proteomes[proteome_names[3]].append(record)

    # For each proteome, concatenate all sequences and write one entry
    for proteome_name, record_list in proteomes.items(): #iterates through the proteomes dictionary storing the SeqI0 entries  
        concat_seq = Seq("".join(str(rec.seq) for rec in record_list)) # Concatenate all sequences as one string
        # Create one SeqRecord with the requested header
        single_record = SeqRecord(
            concat_seq,
            id= proteome_name + ".Metaproteome",
            name="",
            description="",
        )

        # Write to a single‑entry FASTA for this proteome
        outfile = f"{args.output_prefix}_{proteome_name}.fasta"
        SeqIO.write(single_record, outfile, "fasta")
        print(f"Wrote concatenated sequence (n_records={len(record_list)}) to {outfile}")

if __name__ == "__main__":
    main()