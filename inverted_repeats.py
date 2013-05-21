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
    inv_rep_length = 0
    rc_seq = seq.reverse_complement()
    # convert to text format
    seq_txt = str(seq.seq)
    rc_seq_txt = str(rc_seq.seq)
    
    # need to test whether seq is a reverse complement of itself
    pass

    # check if first x bases are a reverse complement of the last x bases
    for i in range(1, len(seq_txt)):
        # slower: if str(seq_txt[0:i].reverse_complement().seq) == str(seq[-i:].seq):
        #if str(rc_seq[-i:].seq) == str(seq[-i:].seq):
        if rc_seq_txt[-i:] == seq_txt[-i:]:
            inv_rep_length = i
        # check when no longer adding bases to inv_rep
        if inv_rep_length > 0 and i > inv_rep_length:
            break
        # check for maximum length to check reached
        if i == longest_length_to_check - 1:
            break
        if inv_rep_length > 0:
            new_seq = seq[:-inv_rep_length]
            inv_rep = str(seq[-inv_rep_length:].seq)
        else:
            new_seq = seq
            inv_rep = ''
    return [new_seq, inv_rep, inv_rep_length]

# <codecell>

def test(seq):
#    [new_seq, inv_rep, inv_rep_length] = check(SeqRecord(Seq(seq)))
#    return [str(new_seq.seq), inv_rep, inv_rep_length]
    print extract_inv_repeat(SeqRecord(Seq(seq)))

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

# first/last 10 RC of eachother
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGAAGCTACGACT') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGA','AGCTACGACT', 10]
# one base
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAT') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTA', 'T', 1]
# no inv_rep
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA') == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA', '', 0]
#test('')

# <codecell>

def process(rec):
    if isSequence(rec):
		new_rec, inv_rep, length = check(rec)
		# don't chenge when length of palindrom is below minimum
		if length < shortest_length_to_check:
			new_rec = rec
		else:
			# change record ID
			new_rec.description = 'cleaned_off_' + inv_rep
    return new_rec, inv_rep, length

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

# <codecell>

#%%timeit  -n 10
# OLD
if __name__ == "__main__":
    path = '/Users/alexajo/Dropbox/current/inverted_repeats/data/'
    #path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/'
    #path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/test_Lex/'
    infile = 'Subset_Lex.fastq'
    #infile = 'Determ.Underterm.collapsed.fastq'
    #infile = 'All_trimmed_forward.fastq'
    #infile = 'temp.fastq'
    #infile = 'ancient_dna_terminal_inv_rep_test.fastq'
    #infile = 'ancient_dna_terminal_inv_rep_test.fastq'

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
            new_rec, inv_rep, inv_rep_length = process(rec)
            out_fh.write(new_rec.format("fastq"))
            out_fh.flush()
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

# <markdowncell>

# original: 10 loops, best of 3: 2.15 s per loop  
# only once reverse_complement: 10 loops, best of 3: 1.89 s per loop  
# only one time str(): 10 loops, best of 3: 591 ms per loop

# <codecell>

def find_inv_repeat(seq, seq_rc):
    """Check for inverted repeat:
       whether first x bases are the reverse complement of the last x bases
       Returns starting position of the inverted repeat (or 0)"""

    # need to test whether seq is a reverse complement of itself
    pass

    inv_rep_length = 0
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

    # expects a Bio.SeqRecord.SeqRecord
    assert type(seq) == SeqRecord, "Not a sequence record: '%s'" % str(seq)
    assert len(seq) > 0, "Sequence length is zero:\n" + str(seq)

    # get sequence and reverse complement as text format
    seq_txt = str(seq.seq)
    rc_seq_txt = str(seq.reverse_complement().seq)
    
    inv_rep_length = find_inv_repeat(seq_txt, rc_seq_txt)
    if inv_rep_length > 0:
        new_seq = seq[:-inv_rep_length]
        inv_rep = str(seq[-inv_rep_length:].seq)
    else:
        new_seq = seq
        inv_rep = ''

    return [str(new_seq.seq), inv_rep, inv_rep_length]

# <codecell>

# 10 bases
test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGAAGCTACGACT')# == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGA', 'AGCTACGACT', 10]
#one_base
test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAT')# == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTA', 'T', 1]
# no inv_rep
test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA')# == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA', '', 0]

# <codecell>

# NEW
#%%timeit -n 10
if __name__ == "__main__":
    path = '/Users/alexajo/Dropbox/current/inverted_repeats/data/'
    #path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/'
    #path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/test_Lex/'
    infile = 'Subset_Lex.fastq'
    #infile = 'Determ.Underterm.collapsed.fastq'
    #infile = 'All_trimmed_forward.fastq'
    #infile = 'temp.fastq'
    #infile = 'ancient_dna_terminal_inv_rep_test.fastq'
    #infile = 'ancient_dna_terminal_inv_rep_test.fastq'

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
            new_rec, inv_rep = extract_inv_repeat(rec)
            inv_rep_length = len(inv_rep)
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


