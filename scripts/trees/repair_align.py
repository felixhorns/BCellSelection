import sys
import os
import shutil
import subprocess
import uuid
import time
from Bio import Phylo, Align, SeqIO, AlignIO

def load_tree(f):
    t = Phylo.read(f, 'newick')
    t.root_with_outgroup("germline")
    t.get_nonterminals()[0].branch_length = 0.0
    t.ladderize(reverse=True)
    return t

def load_aln(infile):
                        
    aln = Align.MultipleSeqAlignment([])
    aln_dict = {}

    with open(infile, 'r') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            aln.append(seq_record)
            aln_dict[seq_record.id] = str(seq_record.seq)

    return aln, aln_dict

def find_long_branches(T, depth):
    long_branches = []    
    depths = T.depths()
    for node in T.get_terminals():
        if depths[node] > depth:
            long_branches.append(node.name)
    return long_branches

def load_aln_to_repair(infile, omit):

    aln = Align.MultipleSeqAlignment([])
    aln_dict = {}

    with open(infile, 'r') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            
            aln_dict[seq_record.id] = str(seq_record.seq)

            if seq_record.name not in omit:
                aln.append(seq_record)
                
    return aln, aln_dict

def muscle_profile_profile_repair1(infile, seq, name):

    path_to_muscle_bin = "/local10G/rfhorns/resources/muscle3.8.31/muscle"

    infile_temp = os.path.join(os.path.dirname(infile), str(uuid.uuid4()) + ".fa~")
    with open(infile_temp, 'w') as out:
        out.write(">" + name + "\n")
        out.write(seq)

    outfile = os.path.join(os.path.dirname(infile), str(uuid.uuid4()) + ".fa~")        
    cmd = path_to_muscle_bin + " -profile -in1 " + infile + " -in2 " + infile_temp + " -out " + outfile + " -quiet -maxiters 2 -diags"
    result = subprocess.call(cmd, shell=True)
    os.remove(infile_temp)
    return outfile

def muscle_profile_profile_repair(infile, long_branches, aln_dict, outfile):
    for name in long_branches:
        next_infile = muscle_profile_profile_repair1(infile, aln_dict[name], name)
        os.remove(infile)
        infile = next_infile
    shutil.move(infile, outfile)
    return outfile

def fasttree(infile, outfile):
    fasttree_bin = "/local10G/rfhorns/resources/FastTree/fasttree"
    cmd = fasttree_bin + " -nt -gtr -out " + outfile + " " + infile
    subprocess.call(cmd, shell=True)
    return outfile

if __name__ == "__main__":

    infile_aln = sys.argv[1]
    infile_tree = sys.argv[2]
    outfile_aln = sys.argv[3]
    outfile_tree = sys.argv[4]
    depth = float(sys.argv[5])

    verbose = True

    # Load alignment and tree
    aln, aln_dict = load_aln(infile_aln)
    T = load_tree(infile_tree)

    print len(aln), "sequences"
    print len(aln[0].seq), "bp"
    print "Depth cutoff", depth

    # Find branches longer than desired depth
    long_branches = find_long_branches(T, depth)

    print "Number of long branches:", len(long_branches)

    # If there is nothing to repair, then copy inputs to outputs
    if len(long_branches) == 0:
        print "Nothing to do."
        shutil.copyfile(infile_aln, outfile_aln)
        shutil.copyfile(infile_tree, outfile_tree)
        print "Done!!"
        quit()

    # Load alignment without long branches
    aln_without_long_branches, _ = load_aln_to_repair(infile_aln, long_branches)

    # Write alignment without long branches to temporary fasta file
    temp_fasta_file = str(uuid.uuid4()) + ".fa~"
    with open(temp_fasta_file, 'w') as out:
        AlignIO.write(aln_without_long_branches, out, 'fasta')

    # Perform realignment of each long branch individually

    if verbose:
        print "Running muscle repair..."

    start_time = time.time()        
    
    muscle_profile_profile_repair(temp_fasta_file, long_branches, aln_dict, outfile_aln)

    if verbose:
        print "Elapsed time (wall clock):", time.time() - start_time, "s"    

    # Infer tree based on repaired alignment

    if verbose:
        print "Running fasttree..."

    start_time = time.time()
    
    fasttree(outfile_aln, outfile_tree)

    if verbose:
        print "Elapsed time (wall clock):", time.time() - start_time, "s"

    print "Done!!"
    
