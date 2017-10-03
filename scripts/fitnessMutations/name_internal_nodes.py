import sys
import ete3

def name_internal_nodes(T):
    i = 1
    for node in T.traverse():
        if not node.is_leaf():
            node.name = str(i) + "_"
            i += 1
    return None

if __name__ == "__main__":

    infile = sys.argv[1]
    outfile = sys.argv[2]

    print infile

    T = ete3.Tree(infile, format=1)
    T.set_outgroup("germline")
    name_internal_nodes(T)
    T.write(format=1, outfile=outfile)

    print "Done!!"