#
# Execute tree generation and picture display for all cluster msa
#

#!/usr/bin/env python3

import glob
import os
import subprocess

pattern = "9.fa*.fa.msa"

for msa_file in glob.glob(pattern):
    # 1) Make tree with belvu
    # Run belvu -Tn -o tree infile > outfile
    root, ext = os.path.splitext(msa_file)      # .msa -> .nwk
    nwk_file = root + ".nwk"

    print(f"Belvu: {msa_file} -> {nwk_file}")
    with open(nwk_file, "w") as out_f:
        subprocess.run(
            ["belvu", "-Tn", "-o", "tree", msa_file],
            stdout=out_f,
            check=True
        )

    # 2) Convert tree to PNG with getPhylotreePNG.py
    # Run python getPhylotreePNG.py -i infile -o outfile -l 1
    png_file = root + ".png"
    print(f"getPhylotreePNG: {nwk_file} -> {png_file}")
    subprocess.run(
        ["python", "getPhylotreePNG.py", "-i", nwk_file, "-o", png_file, "-l", "1"],
        check=True
    )