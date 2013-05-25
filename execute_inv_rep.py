# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import argparse

# <codecell>

try:
    reload(ir)
except:
    import inverted_repeats as ir

# <codecell>

desc = ("""Strips off the 3' copy of any inverted repeat.\n
Input: one fastq file\n
Output:
1) a new fastq file with cleaned sequences
'infile.fastq' generates 'infile.clean.fastq'
2) a file called infile.inv_reps.txt with the stripped sequences and their counts.
3) a file called infile.inv_rep_lengths.txt with the length distribution of the stripped sequences.
""")

# <codecell>

parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-i','--input', help='Input file name',required=True)

if __name__ == "__main__":
    args = parser.parse_args()
    ## show values ##
    print ("Input file: %s" % args.input )
    ir.process(args.input)

# <codecell>


