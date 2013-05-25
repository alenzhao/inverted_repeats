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

parser = argparse.ArgumentParser(description='Help text to be added.')
parser.add_argument('-i','--input', help='Input file name',required=True)

if __name__ == "__main__":
    args = parser.parse_args()
    ## show values ##
    print ("Input file: %s" % args.input )
    print ("Output file: %s" % args.output )

# <codecell>

infile = '/Users/alexajo/Dropbox/current/inverted_repeats/data/Subset_Lex.fastq'
ir.process(infile)

# <codecell>


