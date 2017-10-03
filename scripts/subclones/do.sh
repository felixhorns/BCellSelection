SEEDFILE=seedfile.txt

while read in ; do sbatch annotate_nodes_FayAndWusH.sh $in ; done < $SEEDFILE
