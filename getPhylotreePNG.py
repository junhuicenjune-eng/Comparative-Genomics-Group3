import argparse
from Bio import Phylo
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract top BLAST hit sequence for each genome"
    )
    parser.add_argument("-i","--infile", required=True, type=str, help="Treefile in .nwk format gerated by belvu")
    parser.add_argument("-o","--outfile", required=True, type=str, help="name of outputfile .png")
    parser.add_argument("-l","--label", required=False, type=int, default=0, help="determines displayed branch label; 0 none, 1 branchlenght, 2 branch confidence")
    
    return parser.parse_args()

def main():
    args = parse_args()

    tree = Phylo.read(args.infile, "newick")

    print(set(clade.branch_length for clade in tree.get_terminals()))
    print(set(clade.branch_length for clade in tree.get_nonterminals()))

    for clade in tree.get_nonterminals():
        print(clade, getattr(clade, "confidence", None))

    fig = plt.figure(figsize=(16, 10))
    axes = fig.add_subplot(1, 1, 1)

    if args.label == 1:    
        Phylo.draw(tree, do_show=False, axes=axes, branch_labels=lambda clade: clade.branch_length)
    elif args.label == 2:
        Phylo.draw(tree, do_show=False, axes=axes, branch_labels=lambda clade: getattr(clade, "confidence", None))
    else:
        Phylo.draw(tree, do_show=False, axes=axes)
    axes.set_axis_off()
    plt.savefig(args.outfile, dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()