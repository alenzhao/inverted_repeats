# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

try:
    reload(ir)
except:
    import inverted_repeats as ir

# <codecell>

if __name__ == "__main__":
    #path = '/Users/alexajo/Dropbox/current/inverted_repeats/data/'
    #path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/'
    #path = '/projects/454data/in_progress/bastiaan/GM_historic_DNA/Hiseq_pilot/data/test_Lex/'
    infile = '/Users/alexajo/Dropbox/current/inverted_repeats/data/Subset_Lex.fastq'
    #infile = 'Determ.Underterm.collapsed.fastq'
    #infile = 'All_trimmed_forward.fastq'
    #infile = 'temp.fastq'
    #infile = 'ancient_dna_terminal_palindrome_test.fastq'
    ir.process(infile)

# <codecell>


