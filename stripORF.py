#!/usr/bin/env python3
"""
strip_orf_labels.py

Reads all Newick files matching the pattern '9.fa*.fa.nwk' in a given
directory and writes modified copies named '9.fa*-ORF.fa.nwk'.

For every line, any leaf label of the form
    <PREFIX>_<anything>/...
where PREFIX is one of {9.fa, 11.fa, 22.fa, 27.fa},
is rewritten to
    <PREFIX>/...
i.e. the ORF identifier between the prefix and the first '/' is removed.
"""

import argparse
import glob
import os
import re

# Prefixes we want to clean. The order does not matter because we anchor
# them with a regex word-boundary-like lookbehind below.
PREFIXES = ["9.fa", "11.fa", "22.fa", "27.fa"]

# Regex 1: clean leaf labels.
#
# We match a whole leaf token and replace it with just the genome prefix:
#   27.fa_orf03035_rev/1-297:0.267  ->  27.fa
#   9.fa/1-263:0.590                ->  9.fa
#   9.fa:0.590                      ->  9.fa
# so that (a) names are identical across input trees and (b) branch
# lengths are stripped, leaving a pure-topology cladogram.
#
# Breakdown:
#   (?<![\w.])         negative lookbehind: prefix is not part of a longer
#                      token (prevents '19.fa' matching '9.fa').
#   (9\.fa|11\.fa|...) capture group 1: the prefix itself.
#   [^\s,():;]*        optional ORF/strand suffix and '/start-end' range.
#                      Stops at any Newick structural character so we can
#                      never eat across nodes.
#   (?::[^\s,():;]+)?  optional ':<branchlength>'. The (?:...) is a
#                      non-capturing group so it doesn't disturb \1.
LEAF_PATTERN = re.compile(
    r"(?<![\w.])("
    + "|".join(re.escape(p) for p in PREFIXES)
    + r")[^\s,():;]*(?::[^\s,():;]+)?"
)

# Regex 2: strip the internal-node annotation between ')' and the next
# Newick structural character. This removes both placeholder labels like
# ')0:' and any branch length that follows, e.g. ')0:0.042,' -> '),'.
# Whitespace/newlines between ')' and the annotation are preserved via
# capture group 1 so the file keeps its line breaks.
INTERNAL_ANNOT_PATTERN = re.compile(r"\)(\s*)[^,();\s]+")


def transform_text(text: str) -> str:
    """Clean leaf labels and strip internal-node labels/branch lengths."""
    text = LEAF_PATTERN.sub(r"\1", text)
    text = INTERNAL_ANNOT_PATTERN.sub(r")\1", text)
    return text


def output_name(input_path: str) -> str:
    """
    Map '9.faXYZ.fa.nwk' -> '9.faXYZ-ORF.fa.nwk', preserving the directory.
    We replace the trailing '.fa.nwk' with '-ORF.fa.nwk'.
    """
    directory, base = os.path.split(input_path)
    if not base.endswith(".fa.nwk"):
        raise ValueError(f"Unexpected filename (no .fa.nwk suffix): {base}")
    new_base = base[: -len(".fa.nwk")] + "-ORF.fa.nwk"
    return os.path.join(directory, new_base)


def process_file(input_path: str) -> str:
    """Read input_path, transform the whole text, write the new file."""
    out_path = output_name(input_path)
    with open(input_path, "r") as fin:
        text = fin.read()
    with open(out_path, "w") as fout:
        fout.write(transform_text(text))
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-d", "--directory", default=".",
        help="Directory to search for input files (default: current directory).",
    )
    parser.add_argument(
        "-p", "--pattern", default="9.fa*.fa.nwk",
        help="Glob pattern for input files (default: '9.fa*.fa.nwk').",
    )
    args = parser.parse_args()

    search = os.path.join(args.directory, args.pattern)
    # Exclude already-processed files so re-runs are idempotent.
    inputs = [f for f in sorted(glob.glob(search)) if "-ORF.fa.nwk" not in f]

    if not inputs:
        print(f"No files matched: {search}")
        return

    for path in inputs:
        out = process_file(path)
        print(f"{path}  ->  {out}")


if __name__ == "__main__":
    main()