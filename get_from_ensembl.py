# adjusted from 
# https://raw.github.com/nickloman/benchtop-sequencing-comparison/master/scripts/get_from_sra.py

import os
import sys

fullfn = sys.argv[1]
[species, asm, release] = fullfn.split('.')[0:3]

ftp_path = 'ftp://ftp.ensembl.org/pub/release-%s/fasta/%s/dna/%s' % (release, species.lower(), fullfn)

cmd = "curl --retry 3 -s -O %s" % ftp_path
print cmd
os.system(cmd)
