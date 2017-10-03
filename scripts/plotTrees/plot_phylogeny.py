from __future__ import division
import sys
import json
import time
from ete3 import Tree, TreeStyle, AttrFace, NodeStyle, faces, CircleFace, TextFace
from Bio import SeqIO

#######################################################################
##### Functions for repairing multifurcations
def repair_multifurcations(T, L):
    """ Collapses branches shorter than 1/L (< 1 mutation) and creates multifurcations where needed """
    convert_short_branches_to_zero(T, L)
    delete_zero_length_internal_branches(T)
    delete_zero_length_leafs(T)
    
def convert_short_branches_to_zero(T, L):
    for node in T.traverse():
        if node.dist < 1. / float(L):
            node.dist = 0.0

def delete_zero_length_internal_branches(T):
    for node in T.traverse():
        if not node.is_leaf() and node.dist == 0.0 and node.name != "germline":
            node.delete()
            
def delete_zero_length_leafs(T):
    for node in T.traverse():
        if node.is_leaf() and node.dist == 0.0 and node.name != "germline":
            node.delete()

def delete_long_branches(T, cutoff):
    for node in T.traverse():
        if node.dist > cutoff:
#        if node.is_leaf() and node.dist > cutoff:
            node.delete()
            
#######################################################################

lineage_uid = sys.argv[1]

infile_root = "/Users/lime/Dropbox/quake/Bcell/selection/figures/treePlots/v2/trees/"
infile_tree = infile_root + str(lineage_uid) + "/fasttree.rep.nwk"
infile_aln = infile_root + str(lineage_uid) + "/alignment_ungapped_refined_germline.rep.fasta"
outfile = infile_root + str(lineage_uid) + "/fasttree.rep.nwk.clean.recolored.pdf"

start_time = time.time()

print lineage_uid
print infile_tree
print infile_aln
print outfile

# Get length of alignment
for seq_record in SeqIO.parse(infile_aln, "fasta"):
    L = len(seq_record) # length of alignment (nt)

T = Tree(infile_tree, format=1)
T.set_outgroup("germline")
repair_multifurcations(T, L)
# delete_long_branches(T, 0.2)
T.ladderize()

# T.dist = 0.0 # set germline distance to 0

# T.write(format=1, outfile=tree_file_name+".multifurc")

ts = TreeStyle()
# ts.mode = "c"
ts.scale = 500
ts.optimal_scale_level = "full"
# ts.arc_start = 180 # -180 
# ts.arc_span = 180 # 359
ts.show_leaf_name = False
ts.show_branch_length = False
ts.show_branch_support = False
# ts.root_opening_factor = 0.75
ts.draw_guiding_lines = False
ts.margin_left = 50
ts.margin_right = 50
ts.margin_top = 50
ts.margin_bottom = 50
ts.rotation = 0

path_to_sequence_string_uid_to_isotype_map="/Users/lime/Dropbox/quake/Bcell/selection/figures/treePlots/v2/Bcell_flu_high_res.sequences.isotypeDict.V6_Full.csv"
sequence_string_uid_to_isotype = {}
with open(path_to_sequence_string_uid_to_isotype_map) as f:
    for line in f:
        vals = line.rstrip().split()
        sequence_string_uid_to_isotype[vals[1]] = vals[2]

path_to_isotype_to_color_map = "/Users/lime/Dropbox/quake/Bcell/selection/figures/treePlots/v2/isotype_to_color_dict.json"
with open(path_to_isotype_to_color_map, 'rU') as f:
    isotype_to_color = json.load(f)

def my_layout(node):
    if node.is_leaf():
        node.img_style["size"] = 5
        if node.name != "germline":
            sequence_string_uid = node.name
            isotype = sequence_string_uid_to_isotype[sequence_string_uid]
            color = isotype_to_color[isotype]
            # isotype = node.name.split("~")[1]
            # color = "#666666"
            node.img_style["fgcolor"] = color
        else:
            color = "#000000"
            node.img_style["fgcolor"] = color
    else:
        node.img_style["size"] = 0
        
ts.layout_fn = my_layout

print "Rendering..."

T.render(outfile, w=1200, tree_style=ts)

print "Done!!"
print "Elapsed time (s)", time.time() - start_time
print
