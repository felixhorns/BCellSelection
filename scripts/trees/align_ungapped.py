# Align sequences by concatenating V and J around D

import sys

infile_seq_strings = sys.argv[1]
infile_seqids = sys.argv[2]

# Load list of sequences to align
S = []
with open(infile_seqids) as f:
    for line in f:
        seqid = int(line.rstrip())
        S.append(seqid)
        
# Load sequences

seqids = []
V_seqs = []
D_seqs = []
J_seqs = []

with open(infile_seq_strings) as f:
    for line in f:
        
        vals = line.rstrip().split()
        seqid = int(vals[0])

        if seqid in S:
            
            V_seq = vals[4]
            D_seq = vals[5]
            J_seq = vals[6]

            seqids.append(seqid)
            V_seqs.append(V_seq)
            D_seqs.append(D_seq)
            J_seqs.append(J_seq)

# Find maximum length of V and J
L_max_V = max([len(x) for x in V_seqs])
L_max_J = max([len(x) for x in J_seqs])

# Concatenate V, D, and J sequences with gaps at ends
VDJ_seqs = []

for V_seq, D_seq, J_seq in zip(V_seqs, D_seqs, J_seqs):
    V_seq_padded = "".join(["-"] * (L_max_V - len(V_seq)) + [V_seq])
    J_seq_padded = "".join([J_seq] + ["-"] * (L_max_J - len(J_seq)))
    VDJ_seq = V_seq_padded + D_seq + J_seq_padded
    VDJ_seqs.append(VDJ_seq)

# Print result as fasta
for seqid, VDJ_seq in zip(seqids, VDJ_seqs):
    print ">" + str(seqid)
    print VDJ_seq
