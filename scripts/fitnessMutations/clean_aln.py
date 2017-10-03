# Removes positions in the alignment where all observed sequences have a gap

import sys
import time
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO, Align, AlignIO, Phylo

def load_aln(infile):

    aln = Align.MultipleSeqAlignment([])
    aln_dict = {}

    with open(infile, 'r') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            aln.append(seq_record)
            aln_dict[seq_record.id] = str(seq_record.seq)
            
    return aln, aln_dict

if __name__ == "__main__":

    infile = sys.argv[1]
    outfile = sys.argv[2]

    print "Infile is", infile

    start_time = time.time()
    
    # Load alignment
    aln, aln_dict = load_aln(infile)

    # Find positions to remove (all gaps in observed sequences)

    names_observed = [x for x in aln_dict.keys() if "_" not in x and "germline" not in x] # names of observed sequences
    L = len(aln[0]) # length of alignment (bases)

    positions_to_drop = []
    for i in range(L):
        bases_observed = list(set([aln_dict[name][i] for name in names_observed]))
        if bases_observed == ["-"]:
            positions_to_drop.append(i)

    positions_to_keep = list(set(range(L)) - set(positions_to_drop))

    # Remove positions from sequences
    aln_dict_clean = {}
    for name, seq in aln_dict.items():
        seq_clean = "".join([seq[i] for i in positions_to_keep])
        aln_dict_clean[name] = seq_clean

    # Write to file
    with open(outfile, 'w') as out:
        for name in sorted(aln_dict_clean.keys()):
            seq = aln_dict_clean[name]
            out.write(">" + name + "\n")
            out.write(seq + "\n")

    print "Done!!"
    print "Elapsed time (s)", time.time() - start_time
