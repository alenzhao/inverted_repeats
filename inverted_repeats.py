# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os

# <codecell>

# global parameters
shortest_length_to_check = 4
longest_length_to_check = 200

# <codecell>

def isSequence(seq):
    """Check input is a sequence"""
    test = True if type(seq) is Seq else False
    test = True if len(seq) > 0 else False
    return test

# <codecell>

def test(seq):
    result = extract_inv_repeat(SeqRecord(Seq(seq, IUPAC.unambiguous_dna)))
    return [str(result[0].seq), result[1], result[2]]

# <codecell>

def get_outfnames(infile):
    """Adds '.clean' to input filename (before the extension).
    Example: 'infile.fastq' becomes 'infile.clean.fastq'"""
    out_fnames = []
    in_fname = os.path.basename(infile)
    [in_fname, in_ext] = os.path.splitext(in_fname)
    # fastq outfile
    out_fnames.append(in_fname + '.clean' + in_ext)
    # inv_rep sequences + counts
    out_fnames.append(in_fname + '.inv_reps.txt')
    # inv_rep lengths + counts
    out_fnames.append(in_fname + '.inv_rep_lengths.txt')
    return out_fnames

# <codecell>

def plot(lengths):
    from matplotlib import pylab as plt
    ks = []
    vs = []
    for k, v in lengths.iteritems():
        if k >0:
            ks.append(k)
            vs.append(v)
    plt.plot(ks, vs)
    plt.xlim(1,max(lengths.keys()))

# <markdowncell>

# original: 10 loops, best of 3: 2.15 s per loop  
# only once reverse_complement: 10 loops, best of 3: 1.89 s per loop  
# only one time str(): 10 loops, best of 3: 591 ms per loop  
# after refactoring: 10 loops, best of 3: 193 ms per loop

# <codecell>

def find_inv_repeat(seq, seq_rc):
    """Check for inverted repeat:
       whether first x bases are the reverse complement of the last x bases
       Returns starting position of the inverted repeat (or 0)"""

    inv_rep_length = 0

    # need to test whether seq is a reverse complement of itself
    if seq == seq_rc:
        inv_rep_length = len(seq)
    else:
        # check if first x bases are a reverse complement of the last x bases
        for i in range(1, len(seq)):
            if seq_rc[-i:] == seq[-i:]:
                inv_rep_length = i
            else:
                break
    return inv_rep_length

# <codecell>

def extract_inv_repeat(seq):
    """After finding posisiton of inverted repeat - if any - 
       returns Bio.SeqRecord with the palindromic part removed, 
       and the sequence of the palindromic part"""

    assert shortest_length_to_check > 0, "Shortest length to remove needs to be larger than zero, not %s" % shortest_length_to_check
    # expects a Bio.SeqRecord.SeqRecord
    assert type(seq) == SeqRecord, "Not a sequence record: '%s'" % str(seq)
    assert len(seq) > 0, "Sequence length is zero:\n" + str(seq)

    # get sequence and reverse complement as text format
    seq_txt = str(seq.seq)
    rc_seq_txt = str(seq.reverse_complement().seq)
    
    inv_rep_length = find_inv_repeat(seq_txt, rc_seq_txt)
    if inv_rep_length >= shortest_length_to_check:
        new_seq = seq[:-inv_rep_length]
        inv_rep = str(seq[-inv_rep_length:].seq)
        new_seq.description += ' cleaned_off_' + inv_rep
    else:
        new_seq = seq
        inv_rep = ''

    return [new_seq, inv_rep, inv_rep_length]

# <codecell>

temp_test = shortest_length_to_check
shortest_length_to_check = 1
# first/last 10 RC of eachother
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGAAGCTACGACT') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGA', 'AGCTACGACT', 10]
# one base
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAT') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTA', 'T', 1]
# no inv_rep
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA', '', 0]
# entire sequence it's own reverse complement
assert test('ACACAGGCCTGTGT') == ['', 'ACACAGGCCTGTGT', 14]
shortest_length_to_check = temp_test

# <codecell>

#%%timeit -n 10
if __name__ == "__main__":
    path = '/Users/alexajo/Dropbox/current/inverted_repeats/data/'
    #path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/'
    #path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/test_Lex/'
    infile = 'Subset_Lex.fastq'
    #infile = 'Determ.Underterm.collapsed.fastq'
    #infile = 'All_trimmed_forward.fastq'
    #infile = 'temp.fastq'
    #infile = 'ancient_dna_terminal_palindrome_test.fastq'

    [out_fname, out_rname, out_lname] = get_outfnames(path + infile)
    
    inv_reps = {}
    lengths = {}
    total_trimmed = 0
    processed = 0
    
    max_rec_to_process = 1000
    print "Processing sequences..."
    
    with open(path + out_fname, 'w') as out_fh:
        for rec in SeqIO.parse(path+infile, "fastq"):
            processed += 1
            new_rec, inv_rep, inv_rep_length = extract_inv_repeat(rec)
            out_fh.write(new_rec.format("fastq"))
            if inv_rep_length > shortest_length_to_check:
                inv_reps[inv_rep] = inv_reps.get(inv_rep, 0) +1
            lengths[inv_rep_length] = lengths.get(inv_rep_length, 0) +1
            if inv_rep_length > shortest_length_to_check:
                total_trimmed += 1
            if processed == max_rec_to_process:
                break
        print "Writing summary files..."
        
	with open(path + out_rname, "w") as p_out_fh:
		p_out_fh.write("inverted_repeat\tcount\n")
		for p in inv_reps.keys():
			p_out_fh.write("%s\t%s\n" %(p, inv_reps[p]))
                           
	with open(path + out_lname, "w") as l_out_fh:
		l_out_fh.write("repeat_length\tcount\n")
		for l in sorted(lengths.iterkeys()):
			l_out_fh.write("%s\t%s\n" %(l, lengths[l]))
            
    print "Processed %i records, found %i inverted repeats" % (processed, total_trimmed)

# <codecell>

test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA') # == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA', '', 0]

# <codecell>

test('ACACAGGCCTGTGT')

# <codecell>


