while read in; do echo "$in"; python plot_phylogeny.py "$in"; done < seedfile.txt
