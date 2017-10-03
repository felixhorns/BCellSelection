# Get germline V and J sequence corresponding to a lineage

import sys
from Bio import SeqIO

infile = sys.argv[1]
lineage_uid = int(sys.argv[2])
infile_V_dict = sys.argv[3]
infile_J_dict = sys.argv[4]

# Load V and J germline sequence dictionaries

germline_V_seqs = {}
records = list(SeqIO.parse(infile_V_dict, "fasta"))
for record in records:
    germline_V_seqs[record.id] = str(record.seq).upper()

germline_J_seqs = {}
records = list(SeqIO.parse(infile_J_dict, "fasta"))
for record in records:
    germline_J_seqs[record.id] = str(record.seq).upper()

# Get V and J germline gene names and look up sequences

with open(infile) as f:
    for line in f:
        
        vals = line.rstrip().split()
        # seqid = vals[0]
        my_lineage_uid = int(vals[10])
        # sequence_string_uid = int(vals[8])

        if my_lineage_uid == lineage_uid:

            V_gene = vals[25]
            J_gene = vals[27]

            germline_V_seq = germline_V_seqs[V_gene]
            germline_J_seq = germline_J_seqs[J_gene]
            germline_VJ_seq = germline_V_seq + germline_J_seq

            print ">germline"
            print germline_VJ_seq
            quit()
