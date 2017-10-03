import sys
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO, Align, AlignIO, Phylo
from itertools import izip

def load_tree(f):
    t = Phylo.read(f, 'newick')
    t.root_with_outgroup("germline")
    t.get_nonterminals()[0].branch_length = 0.0
    # t.ladderize(reverse=True)
    return t

def load_aln(infile):

    aln = Align.MultipleSeqAlignment([])
    aln_dict = {}

    with open(infile, 'r') as f:
        for seq_record in SeqIO.parse(f, 'fasta'):
            aln.append(seq_record)
            aln_dict[seq_record.id] = str(seq_record.seq)

    return aln, aln_dict

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def str_diffs(X, Y):
    diffs = []
    for i, (x, y) in enumerate(izip(X,Y)):
        if x != y:
            d = [i, x, y]
            diffs.append(d)
    return diffs

def get_mutations(T, aln_dict):
    """ Get mutations on each branch of the tree """
    header = ["name", "parent_name", "position", "base_before", "base_after"]
    df = pd.DataFrame(columns=header)
    i = 0
    for clade in T.find_clades():
        if clade.name in [None, "germline", "2_"]: continue
        parent = get_parent(T, clade)
        seq_parent = aln_dict[parent.name]
        seq_clade = aln_dict[clade.name]
        diffs = str_diffs(seq_parent, seq_clade)
        for diff in diffs:
            position, base_before, base_after = tuple(diff)
            features = [clade.name, parent.name, position, base_before, base_after]
            df.loc[i] = features
            i += 1
    return df

def find_frame(s):
    # Finds longest ORF
    s = s.replace("-", "")
    seq1 = Seq(s, generic_dna).translate() # translate in every frame
    seq2 = Seq(s[1:], generic_dna).translate()
    seq3 = Seq(s[2:], generic_dna).translate()
    L_seq1 = max([len(x) for x in seq1.split("*")]) # find longest ORF in each frame
    L_seq2 = max([len(x) for x in seq2.split("*")])
    L_seq3 = max([len(x) for x in seq3.split("*")])
    Ls = [L_seq1, L_seq2, L_seq3]
    L_max = max(Ls) # get longest ORF among all frames
    frames_max = [i for i, x in enumerate(Ls) if x == L_max] # get frame of longest ORF
    if len(frames_max) > 1:
        print "Warning: more than one reading frame had max length ORF"
    return frames_max[0]

def annotate_coding(df_mutations, aln_dict):
    """ Annotate each mutation as either nonsynonymous or synonymous """
    coding_status = []
    for i, row in df_mutations.iterrows():
        seq_parent = aln_dict[row["parent_name"]] # get parent sequence
        seq_mutated = list(seq_parent)
        seq_mutated[int(row["position"])] = row["base_after"] # introduce mutation
        seq_mutated = "".join(seq_mutated)
        seq_parent = seq_parent.replace("-", "") # collapse gaps
        seq_mutated = seq_mutated.replace("-", "")
        AA_parent = Seq(seq_parent, generic_dna).translate() # translate
        AA_mutated = Seq(seq_mutated, generic_dna).translate()
        if AA_parent != AA_mutated: # compare AA before and after mutation
            coding_status.append("N")
        else:
            coding_status.append("S")
    df_mutations["coding_status"] = coding_status
    return df_mutations

def map_positions(s, positions):
    """ Maps positions in an ungapped sequence to corresponding positions in a gapped sequence """
    # count number of gaps before each position
    counter = 0
    gaps = []
    for i, x in enumerate(s):
        if x == "-":
            counter += 1
        gaps.append(counter)

    # transform boundaries to corresponding positions in new sequence
    positions_transformed = []
    for x in positions:
        my_gaps = gaps[x]
        x_transformed = x + my_gaps
        positions_transformed.append(x_transformed)
    return positions_transformed

def annotate_regions(df_mutations, aln_dict, df_seqs):
    """ Annotate region of each mutation (CDR/FWR) """
    # get one sequence
    sequence_uids = [x for x in aln_dict.keys() if "_" not in x]
    my_sequence_uid = int(sequence_uids[0])
    s = aln_dict[str(my_sequence_uid)]
    # transform positions of region boundaries to corresponding positions in gapped alignment
    fields = ["FWR1_start", "CDR1_start", "FWR2_start", "CDR2_start", "FWR3_start", "CDR3_start", "FWR4_start", "C_start"]
    boundaries_ungapped = df_seqs.loc[my_sequence_uid][fields]
    boundaries_ungapped = np.array(boundaries_ungapped) - 1 # transform to zero-indexed positions
    boundaries_ungapped[-1] -= 1 # decrement C region boundary (end of sequence) to fit within array
    boundaries_gapped = map_positions(s, boundaries_ungapped)
    # boundaries_gapped = np.array(boundaries_gapped) - 1 # not used anymore (we do transform earlier)
    boundaries_gapped[0] = 0
    boundaries_gapped[-1] += 1
    # map mutations to regions using boundaries
    labels = ["FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3", "FWR4"]
    regions = pd.cut(df_mutations["position"], boundaries_gapped, include_lowest=True, right=False, labels=labels)
    df_mutations["region"] = regions
    return df_mutations

if __name__ == "__main__":

    infile_aln = sys.argv[1]
    infile_fitness_tree = sys.argv[2]
    outfile = sys.argv[3]

    print infile_fitness_tree
    print infile_aln
    
    infile_df_seqs = "/local10G/rfhorns/Bcell/flu_highres/figures/v5/data/FitnessMutations.df_seqs_raw.csv"
    df_seqs = pd.read_csv(infile_df_seqs, header=0, index_col=0)

    aln, aln_dict = load_aln(infile_aln)
    fitness_tree = load_tree(infile_fitness_tree)

    df_mutations = get_mutations(fitness_tree, aln_dict)
    df_mutations = annotate_coding(df_mutations, aln_dict)
    df_mutations = annotate_regions(df_mutations, aln_dict, df_seqs)

    df_mutations.to_csv(outfile)

    print "Done!!"
