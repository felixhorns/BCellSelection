import sys
import pandas as pd
from Bio import SeqIO, Align, AlignIO, Phylo

import matplotlib
matplotlib.use('Agg')
import pylab as plt

sys.path.append("/local10G/rfhorns/resources/FitnessInference/prediction_src/")
import sequence_ranking
import tree_utils

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

def get_fitness_changes(fit):
    """ Get fitness changes on each branch of the tree """
    header = ["name", "depth", "mean_fitness", "var_fitness",
              "parent_name", "parent_depth", "parent_mean_fitness", "parent_var_fitness",
              "delta_mean_fitness", "length"]
    df = pd.DataFrame(columns=header)
    i = 0
    depths = fit.T.depths()
    for clade in fit.T.find_clades():
        if clade.name in [None, "germline", "2_"]: continue
        parent = get_parent(fit.T, clade)
        delta = clade.mean_fitness - parent.mean_fitness
        features = [clade.name, depths[clade], clade.mean_fitness, clade.var_fitness,
                    parent.name, depths[parent], parent.mean_fitness, parent.var_fitness,
                    delta, depths[clade] - depths[parent]]
        df.loc[i] = features
        i += 1
    df.set_index("name", inplace=True)
    return df

if __name__ == "__main__":

    infile_tree = sys.argv[1]
    infile_aln = sys.argv[2]
    outfile_fitness_tree = sys.argv[3]
    outfile_df_fitness = sys.argv[4]
    outfile_tree_pdf = sys.argv[5]
    outfile_tree_pdf_labeled = sys.argv[6]

    print infile_tree
    print infile_aln

    # Load alignment
    aln = Align.MultipleSeqAlignment([])
    outgroup = None
    with open(infile_aln) as f:
        for record in SeqIO.parse(f, 'fasta'):
            if record.name != "germline":
                aln.append(record)
            else:
                outgroup = record

    if outgroup is None:
        print "outgroup not in alignment -- FATAL"
        quit()
        
    # Load tree
    T = load_tree(infile_tree)

    # Create sequence data object
    seq_data = sequence_ranking.alignment(aln, outgroup, build_tree=False)
    seq_data.T = T # use predetermined tree

    # Infer fitness
    eps_branch_length = 1e-4
    diffusion = 0.5
    distance_scale = 2.0
    samp_frac = 0.1

    prediction = sequence_ranking.sequence_ranking(seq_data, eps_branch_length=eps_branch_length, pseudo_count = 5,
                                methods = ['mean_fitness'], D=diffusion,
                                distance_scale=distance_scale, samp_frac=samp_frac)

    best_node = prediction.predict()

    # Write fitness tree to file (for mutation annotation)
    Phylo.write(prediction.T, outfile_fitness_tree, "newick")

    # Get fitness changes on each branch
    df_fitness = get_fitness_changes(prediction)
    df_fitness.sort_values("delta_mean_fitness", ascending=False, inplace=True)

    # Write fitness changes to file
    df_fitness.to_csv(outfile_df_fitness)

    # Plot tree colored by fitness with and without node labels
    tree_utils.plot_prediction_tree(prediction)
    plt.savefig(outfile_tree_pdf)

    tree_utils.plot_prediction_tree(prediction, node_label_func=lambda x: x.name)
    plt.savefig(outfile_tree_pdf_labeled)

    print "Done!!"