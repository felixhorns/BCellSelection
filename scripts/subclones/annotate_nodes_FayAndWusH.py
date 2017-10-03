from __future__ import division
import sys
import pickle

import numpy as np
from scipy.stats import johnsonsu
import pandas as pd

import ete3

def name_internal_nodes(T):
    # Name internal nodes of a tree
    i = 1
    for node in T.traverse():
        if not node.is_leaf():
            node.name = str(i) + "_"
            i += 1
    return None

def flatten(L):
    """ Convert list of lists into flat list"""
    return [item for sublist in L for item in sublist]

def unpack(X):
    """ Unpack a comma separated list of values into a flat list """
    return flatten([x.split(",") for x in list(X)])

def get_muts_counts_subtree(df_seqs, lineage_uid, T, subtree_parent, seq_string_uid_to_uid):
    """ Get unique somatic mutations and occurrence of each for a subtree of a lineage """

    leaf_names = [node.name for node in T.search_nodes(name=subtree_parent)[0].get_leaves()]
    leaf_names = [x for x in leaf_names if x != "germline"] # remove germline
    leaf_seq_string_uids = [int(x.split("~")[0]) for x in leaf_names] # get sequence string uids
    leaf_uids = [seq_string_uid_to_uid[x] for x in leaf_seq_string_uids] # convert to sequence uids

    df = df_seqs.loc[leaf_uids] # subset to only leaves in the subtree
    N = df.shape[0] # number of sequences in subtree
    
    df = df[df["mut_germline_positions"] != ""] # filter out sequences without germline mutations
    df = df.loc[-pd.isnull(df["mut_germline_positions"])] # filter out sequences without germline mutations

    # Count occurrences of each mutation
    muts = zip(unpack(df.mut_germline_positions), unpack(df.mut_germline_before), unpack(df.mut_germline_after))
    counts = {}
    for mut in muts:
        try:
            counts[mut] += 1
        except:
            counts[mut] = 1
            
    return counts, N

def calc_H(mut_counts, n):
    # Calculate Fay and Wu's H based on counts of mutations
    counts = pd.Series(mut_counts).value_counts()
    theta_H = sum(2 * np.array(counts.index)**2 * counts) / (n * (n-1))
    theta_pi = sum(2 * counts * np.array(counts.index) * (n - np.array(counts.index))) / (n * (n-1))
    H = theta_pi - theta_H
    return H

def find_nearest(L,value):
    array = np.array(L)
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def calc_pvalue_matchedSimulations(H_focal, N_focal, params, model=johnsonsu):
    N = find_nearest(params.keys(), N_focal) # Find nearest N in ensemble
    myParams = params[N]
    p_low = model.cdf(H_focal, myParams[0], myParams[1], myParams[2], myParams[3]) # Calculate p of H under nearest N
    p = p_low
    return p

def calc_FayAndWusH_subtree(T, node_name, lineage_uid, df_seqs,
                            fit_params_kingman, fit_params_BSC,
                            seq_string_uid_to_uid):
    # Calculate Fay and Wu's H and P values for H under neutral models for a subtree
    
    # Get counts of derived mutations in subtree
    mut_counts, N = get_muts_counts_subtree(df_seqs, lineage_uid, T, node_name, seq_string_uid_to_uid)

    # Remove mutations at fixation within subtree
    for mut, count in mut_counts.items():
        if count == N:
            mut_counts.pop(mut) # mutation is fixed in subtree, drop it

    # Calculate Fay and Wu's H based on mutations within subtree
    H = calc_H(mut_counts.values(), N)

    # Calculate P values based on comparison with simulations of neutrality and strong selection
    pvalue_kingman = calc_pvalue_matchedSimulations(H, N, fit_params_kingman, model=johnsonsu)
    pvalue_BSC = calc_pvalue_matchedSimulations(H, N, fit_params_BSC, model=johnsonsu)
    
    return H, pvalue_kingman, pvalue_BSC

def annotate_FayAndWusH(T, lineage_uid, df_seqs, fit_params_kingman, fit_params_BSC, seq_string_uid_to_uid):
    """ Traverse tree, calculate Fay and Wu's H for each node, and calculate significance """

    annotations = []

    # Condition for stopping traversal    
    def stop(node):
        if node.name == "germline" or node.name == "1_":
            return False
        if len(node) < 100:
            # print "Stopping branch (too small)", node.name
            return True
        else:
            return False

    # Traverse tree and calculate Fay and Wu's H at each node
    for node in T.traverse(is_leaf_fn=stop):
        if node.is_leaf() or node.name == "1_": continue
        H, pvalue_kingman, pvalue_BSC = calc_FayAndWusH_subtree(T, node.name, lineage_uid, df_seqs,
                                                                fit_params_kingman, fit_params_BSC,
                                                                seq_string_uid_to_uid)
        myAnnotation = [node.name, len(node), node.dist, H, pvalue_kingman, pvalue_BSC]
        annotations.append(myAnnotation)

    return annotations

if __name__ == "__main__":

    lineage_uid = int(sys.argv[1])

    print "Lineage uid", lineage_uid

    # Set paths to data
    infile_basename = "/local10G/rfhorns/Bcell/flu_highres/figures/v5/data"
    infile_df_seqs = infile_basename+"/FindSubclonesSelected.df_seqs_raw.csv"
    infile_seq_string_uid_to_uid = infile_basename+"/FindSubclonesSelected.seq_string_uid_to_uid.pickle"
    infile_fit_params_kingman = infile_basename+"/fit_params_kingman.pickle"
    infile_fit_params_BSC = infile_basename+"/fit_params_BSC.pickle"    

    # Load sequence data
    df_seqs = pd.read_csv(infile_df_seqs, index_col=0, header=0)
    seq_string_uid_to_uid = pickle.load(open(infile_seq_string_uid_to_uid))

    # Load parameters of fits to null distributions
    fit_params_kingman = pickle.load(open(infile_fit_params_kingman))
    fit_params_BSC = pickle.load(open(infile_fit_params_BSC))

    # Set name of tree file
    infile_tree = "/local10G/rfhorns/Bcell/flu_highres/trees/" + str(lineage_uid) + "/fasttree.rep.nwk"

    print "infile_tree", infile_tree

    # Set name of output file
    outfile = "/local10G/rfhorns/Bcell/flu_highres/trees/" + str(lineage_uid) + "/annotation_FayAndWusH.csv"

    print "outfile", outfile

    # Load tree
    T = ete3.Tree(infile_tree, format=1)

    print "Leaves", len(T)

    # Set names of internal nodes (e.g., 1_)
    name_internal_nodes(T)

    # Annotate nodes with Fay and Wu's H and P value 
    annotations = annotate_FayAndWusH(T, lineage_uid, df_seqs,
                                      fit_params_kingman, fit_params_BSC,
                                      seq_string_uid_to_uid)

    # Convert output to dataframe
    cols = ["name", "N_leaves", "dist", "H", "pvalue_kingman", "pvalue_BSC"]
    df_annotations = pd.DataFrame(columns=cols)
    for i in range(len(annotations)):
        df_annotations.loc[i] = annotations[i]
    df_annotations.set_index("name", inplace=True)
    
    # Write output to file
    df_annotations.to_csv(outfile)

    print "Done!!"
