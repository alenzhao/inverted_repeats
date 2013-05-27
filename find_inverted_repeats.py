# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from clip_inverted_repeats import find_inv_repeat
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# <codecell>

revcomp = {}
a = 'ACGT'
for i in a:
    for j in a:
        for k in a:
            revcomp[i+j+k] = str(Seq(i+j+k).reverse_complement())

# <codecell>

seq = 'TGTGTGTGTGTATCGTGTGGGGTTGTGTTTGGCTGCGGTCGGTCCACACACACA' + 'CTAATATTTGCCAGCATGCAAGCGAAAAATATTAG' + 'AAATGGGCCCACATGCGTGCCTGTGTGTGTGTTTGTGGGCCCATTT'
#seq = 'AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGAAGCTACGACT' # first test case from clip_inverted_repeats
len(seq)

# <codecell>

startpos = 54 # 89 # 54 0
trimer = seq[startpos:startpos + 3]
trimer

# <codecell>

hit = seq.rfind(revcomp[trimer]) + len(trimer)
len(seq), hit

# <codecell>

subseq = seq[startpos:hit]
test_for_inv_repeat(subseq)

# <codecell>

def test_for_inv_repeat(subseq):
    inv_rep_len = find_inv_repeat(subseq,str(Seq(subseq).reverse_complement()))
    if inv_rep_len > 0:
        print subseq, inv_rep_len
        print subseq[:inv_rep_len] + "-" + subseq[inv_rep_len:-inv_rep_len]  + "-" + subseq[-inv_rep_len:]
    else:
        print "No hit"

# <codecell>

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

# <codecell>

seq[0:4], seq[0:48][-4:]

# <codecell>

seq2 = 'TGTGTGTGTGTATCGTGTGGGGTTGTGTTTGGCTGCGGTCGGTCCACACACACA'
inv_rep_len = find_inv_repeat(seq2, reverse_complement(seq2))
print seq2[-inv_rep_len:]

# <codecell>


# <codecell>

revcomp

# <codecell>


