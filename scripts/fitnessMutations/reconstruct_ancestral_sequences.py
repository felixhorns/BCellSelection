import sys
import shutil
import pickle
from Bio import SeqIO, Align, AlignIO, Phylo

sys.path.append("/local10G/rfhorns/resources/FitnessInference/prediction_src/")
import ancestral

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

if __name__ == "__main__":

    infile_tree = sys.argv[1]
    infile_aln = sys.argv[2]
    outfile_aln = sys.argv[3]
    # outfile_aln_pickle = sys.argv[4]

    print infile_tree
    print infile_aln

    T = load_tree(infile_tree)
    aln, _ = load_aln(infile_aln)

    # Do inference
    my_ancestral = ancestral.ancestral_sequences(T, aln)
    my_ancestral.calc_ancestral_sequences()

    # Write result to file
    shutil.copyfile(infile_aln, outfile_aln)
    with open(outfile_aln, 'a') as out:
        for node in my_ancestral.T.get_nonterminals():
            if node.name != None:
                out.write(">" + node.name + "\n")
                out.write(str(node.seq) + "\n")

    # Write alignment object to file (to keep likelihoods)
    # with open(outfile_aln_pickle, 'w') as out:
    #     pickle.dump(my_ancestral, out)

    print "Done!!"