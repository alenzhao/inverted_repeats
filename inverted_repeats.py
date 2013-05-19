# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

# <codecell>

# global parameters
shortest_length_to_check = 3
longest_length_to_check = 200

# <codecell>

def isSequence(seq):
    """Check input is a sequence"""
    test = True if type(seq) is Seq else False
    test = True if len(seq) > 0 else False
    return test

# <codecell>

def check(seq):
    """Check whether first x bases are the reverse complement of the last x bases
       Returns sequence with the palindromic part removed, and the palindromic part"""
    # expects a Bio.SeqRecord.SeqRecord
    assert type(seq) == SeqRecord, "Not a sequence record: " + str(seq)
    assert len(seq) > 0, "Sequence length is zero: " + str(seq)
    palindrome_length = 0
    # need to test whether seq is a reverse complement of itself
    pass
    # check if first x bases are a reverse complement of the last x bases
    for i in range(len(seq)):
        subseq = seq[0:i]
        #print str(rec[-10:].reverse_complement().seq) == str(rec[0:10].seq)
        if str(subseq.reverse_complement().seq) == str(seq[-i:].seq):
            #print i, subseq.reverse_complement(), seq[-i:]
            palindrome_length = len(subseq)
        # check when no longer adding bases to palindrome
        if palindrome_length > 0 and i > palindrome_length:
            break
        # check for maximum length to check reached
        if i == longest_length_to_check - 1:
            break
        if palindrome_length > 0:
            new_seq = seq[:-palindrome_length]
            palindrome = str(seq[-palindrome_length:].seq)
        else:
            new_seq = seq
            palindrome = ''
    return [new_seq, palindrome, palindrome_length]

# <codecell>

def test(seq):
    [new_seq, palindrome, palindrome_length] = check(SeqRecord(Seq(seq)))
    return [str(new_seq.seq), palindrome, palindrome_length]

# <codecell>

def get_outfnames(infile):
    """Adds '.clean' to input filename (before the extension).
    Example: 'infile.fastq' becomes 'infile.clean.fastq'"""
    out_fnames = []
    in_fname = os.path.basename(infile)
    [in_fname, in_ext] = os.path.splitext(in_fname)
    # fastq outfile
    out_fnames.append(in_fname + '.clean' + in_ext)
    # palindrome sequences + counts
    out_fnames.append(in_fname + '.palindromes.txt')
    # palindrome lengths + counts
    out_fnames.append(in_fname + '.palindrome_lengths.txt')
    return out_fnames

# <codecell>

# first/last 10 RC of eachother
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGAAGCTACGACT') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGA','AGCTACGACT', 10]
# one base
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAT') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTA', 'T', 1]
# no palindrome
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA', '', 0]
#test('')

# <codecell>

def process(rec):
    if isSequence(rec):
		new_rec, palindrome, length = check(rec)
		# don't chenge when length of palindrom is below minimum
		if length < shortest_length_to_check:
			new_rec = rec
		else:
			# change record ID
			new_rec.description = '_cleaned_off_' + palindrome
			#print new_rec
			# log palindrome data
			palindromes[palindrome] = palindromes.get(palindrome, 0) +1
		lengths[length] = lengths.get(length, 0) +1
    return new_rec

# <codecell>

from matplotlib import pylab as plt
ks = []
vs = []
for k, v in lengths.iteritems():
    if k >0:
        ks.append(k)
        vs.append(v)

plt.plot(ks, vs)
#plt.xlim(1,max(lengths.keys()))

# <codecell>

if __name__ == "__main__":
    path = '/Users/alexajo/Dropbox/current/inverted_repeats/data/'
#    path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/'
    infile = 'Subset_Lex.fastq'
    #infile = 'Determ.Underterm.collapsed.fastq'
    #infile = 'All_trimmed_forward.fastq'
    #infile = 'temp.fastq'
    #infile = 'ancient_dna_terminal_palindrome_test.fastq'
    #infile = 'ancient_dna_terminal_palindrome_test.fastq'

    [out_fname, out_pname, out_lname] = get_outfnames(path + infile)
    
    palindromes = {}
    lengths = {}
    total_trimmed = 0
    processed = 0

    max_rec_to_process = 1000
    with open(path + out_fname, 'w') as out_fh:
        for rec in SeqIO.parse(path+infile, "fastq"):
            processed += 1
            new_rec = process(rec)
            out_fh.write(new_rec.format("fastq"))
            if processed == max_rec_to_process:
                break
    print "Processed %i records" % processed

# <codecell>


