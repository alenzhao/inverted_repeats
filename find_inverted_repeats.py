# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from clip_inverted_repeats import find_inv_repeat
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# <codecell>

revcomp = {}
seed = 3 # number of bp to start searching for inv_rep
a = 'ACGT'
for i in a:
    for j in a:
        for k in a:
            revcomp[i+j+k] = str(Seq(i+j+k).reverse_complement())

# <codecell>

#def reverse_complement(seq):
#    return str(Seq(seq).reverse_complement())

# <codecell>

seq = 'TGTGTGTGTGTATCGTGTGGGGTTGTGTTTGGCTGCGGTCGGTCCACACACACA' + 'CTAATATTTGCCAGCATGCAAGCGAAAAATATTAG' + 'AAATGGGCCCACATGCGTGCCTGTGTGTGTGTTTGTGGGCCCATTT'
#seq = 'AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGAAGCTACGACT' # first test case from clip_inverted_repeats
len(seq)

# <codecell>

def test_for_inv_repeat(subseq):
    inv_rep_len = find_inv_repeat(subseq,str(Seq(subseq).reverse_complement()))
    if inv_rep_len > 0:
        return str(inv_rep_len) + ": " + subseq[:inv_rep_len] + "-" + subseq[inv_rep_len:-inv_rep_len]  + "-" + subseq[-inv_rep_len:]
    else:
        return "No hit"

# <codecell>

for pos in range(len(seq) - seed + 1):
    mer = seq[pos:pos + seed]
    hit = seq.rfind(revcomp[mer]) + len(mer)
    subseq = seq[pos:hit]
    print str(pos) + '--> ' + test_for_inv_repeat(subseq)

