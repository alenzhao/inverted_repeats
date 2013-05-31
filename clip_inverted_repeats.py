# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Some Illumina reads from ancient DNA samples were found to consist of an inverted repeat
# with a different seqeunce inbetween
# in other words, the first x bases of a read are the reverse complement of the last x bases
# This script is meant to clip the 3' part of an inverted repeat when present
# A new fastq file is generated mentioning in the sequences ID which sequence was clipped, if any
# Two metrics files on the repeats found (and clipped) are proeduced as well
# When a sequence is its own reverse complement, this does not get clipped
# but a mention is made sequence identifier and
# these are also reported in the metrics file
#
# Written by Lex Nederbragt, with input from Bastiaan Star
# Version 1.0 release candidate 1, May 2013
#
# on Abel, needs 'module load python2'
# run as 'python script_name.py -h' for instructions

# <codecell>

# Not implemented yet:
# longest_length_to_check
# maximum number of reads to process as command-line argument

# <codecell>

from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
import argparse

# <codecell>

# help text and argument parser
desc = '\n'.join(["Strips off the 3' copy of any inverted repeat.",
                 "Input: one fastq file.",
                 "Output:",
                 "1) a new fastq file with cleaned sequences: 'infile.fastq' gives 'infile.clean.fastq'",
                 "2) a file called 'infile.inv_reps.txt' with the stripped sequences and their counts",
                 "3) a file called 'infile.inv_rep_lengths.txt' with the length distribution of the stripped sequences.",
                 "An optional argument -s/--shortest_length can be used to set the minumum length of repeat to clip (default: 4 bp)"
                  ])
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-i','--input', help='Input file name',required=True)
parser.add_argument('-s', '--shortest_length', help='Shortest repeat length to clip', type=int, default=4, required = False)

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

def extract_inv_repeat(seq, shortest_length_to_clip):
    """After finding position of inverted repeat - if any - 
       returns Bio.SeqRecord with the inverted repeated part removed, 
       and the sequence of the removed part"""

    assert shortest_length_to_clip > 0, "Shortest length to remove needs to be larger than zero, not %s" % shortest_length_to_clip
    # expects a Bio.SeqRecord.SeqRecord
    assert type(seq) == SeqRecord, "Not a sequence record: '%s'" % str(seq)

    # get sequence and reverse complement as text format
    seq_txt = str(seq.seq)
    rc_seq_txt = str(seq.reverse_complement().seq)
    
    # locate the inverted repeat - if any
    inv_rep_length = find_inv_repeat(seq_txt, rc_seq_txt)
    
    # process results
    if inv_rep_length == len(seq_txt):
        # sequence is its own reverse complement
        new_seq = seq
        inv_rep = seq_txt
        new_seq.description += ' self_reverse_complement'
    elif inv_rep_length >= shortest_length_to_clip:
        # hit
        new_seq = seq[:-inv_rep_length]
        inv_rep = str(seq[-inv_rep_length:].seq)
        new_seq.description += ' cleaned_off_' + inv_rep
    else:
        # either no hit, or hit shorter than minimum to report
        new_seq = seq
        inv_rep = ''

    return [new_seq, inv_rep, inv_rep_length]

# <codecell>

def test(seq, shortest_length_to_clip):
    """Performs the 'unit tests'"""
    result = extract_inv_repeat(SeqRecord(Seq(seq, IUPAC.unambiguous_dna)), shortest_length_to_clip)
    return [str(result[0].seq), result[1], result[2]]

# <codecell>

# set of 'unit tests'
# first/last 10 RC of eachother
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGAAGCTACGACT', 1) == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGA', 'AGCTACGACT', 10]
# one base
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAT', 1) == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTA', 'T', 1]
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAT', 4) == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAT', '', 1]
# no inv_rep
assert test('AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA', 1) == ['AGTCGTAGCTGATGCTTAGGGGCTTACTAGGCTTGATGAGGATTAA', '', 0]
# entire sequence it's own reverse complement
assert test('ACACAGGCCTGTGT', 1) == ['ACACAGGCCTGTGT', 'ACACAGGCCTGTGT', 14]
# empty sequence
assert test('', 4) == ['', '', 0]

# <codecell>

def process(infile, shortest_length_to_clip):
    """ Does the actual work:
        Goes through the input file and streams the content through the inverted_repeat locator
        Collects the new sequences and repeats found and reports them """
    
    # test for existing inut file
    assert os.path.exists(infile), "Input file '%s' appears not to exist." %infile
    [out_fname, out_rname, out_lname] = get_outfnames(infile)
    
    inv_reps = {}
    lengths = {}
    total_trimmed = 0
    total_skipped = 0
    processed = 0
    
    max_rec_to_process = 1e30
    print "Processing sequences..."
    
    with open(out_fname, 'w') as out_fh:
        for rec in SeqIO.parse(infile, "fastq"):
            processed += 1
            if len(rec) < 1:
                # skip zero-length sequences
                total_skipped += 1
                continue
            new_rec, inv_rep, inv_rep_length = extract_inv_repeat(rec, shortest_length_to_clip)
            out_fh.write(new_rec.format("fastq"))
            if inv_rep_length >= shortest_length_to_clip:
                inv_reps[inv_rep] = inv_reps.get(inv_rep, 0) +1
                total_trimmed += 1
            lengths[inv_rep_length] = lengths.get(inv_rep_length, 0) +1
            if processed == max_rec_to_process:
                break
        print "Writing summary files..."
        
	with open(out_rname, "w") as p_out_fh:
		p_out_fh.write("inverted_repeat\tcount\n")
		for p in inv_reps.keys():
			p_out_fh.write("%s\t%s\n" %(p, inv_reps[p]))
                           
	with open(out_lname, "w") as l_out_fh:
		l_out_fh.write("repeat_length\tcount\n")
		for l in sorted(lengths.iterkeys()):
			l_out_fh.write("%s\t%s\n" %(l, lengths[l]))
            
    print "\nProcessed %i records:\n- skipped %i (zero-length)\n- found %i inverted repeat(s)\n" % (processed, total_skipped, total_trimmed)

# <codecell>

if __name__ == "__main__":
    args = parser.parse_args()
    print ("Input file: %s" % args.input )
    print ("Shortest length to clip: %s" % args.shortest_length)
    process(args.input, args.shortest_length)

