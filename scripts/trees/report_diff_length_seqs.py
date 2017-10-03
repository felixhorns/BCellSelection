""" Remove columns containing all gaps from a FASTA alignment """

import sys
from Bio import AlignIO

infile = sys.argv[1]

# aln = AlignIO.read(infile, "fasta")
#quit()

L = 0
with open(infile) as f:
    for line in f:
        L += 1

i = 0

len_prev_seq = 0

with open(infile) as f:
    while i < L:

        name = f.readline().rstrip()
        seq = f.readline().rstrip()

        if len(seq) != len_prev_seq:
            print name
            print seq
            print

        len_prev_seq = len(seq)
        
        i += 2
