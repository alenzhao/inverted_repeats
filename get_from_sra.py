# adjusted from 
# https://raw.github.com/nickloman/benchtop-sequencing-comparison/master/scripts/get_from_sra.py

import os
import sys

fullfn = sys.argv[1]
dir = fullfn[0:6]
fn = fullfn[0:9]

cmd = "curl -s -S -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s.fastq.gz" % (dir, fn, fullfn)
print cmd
os.system(cmd)
