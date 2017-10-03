# Get all sequence names for a lineage

import sys

infile = sys.argv[1]
lineage_uid = int(sys.argv[2])

with open(infile) as f:
    for line in f:
        
        vals = line.rstrip().split()
        # seqid = vals[0]
        my_lineage_uid = int(vals[10])
        sequence_string_uid = int(vals[8])

        if my_lineage_uid == lineage_uid:
            
            print sequence_string_uid
