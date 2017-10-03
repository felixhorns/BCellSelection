#!/bin/bash
#
#SBATCH --job-name annotate_FayAndWusH
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=6GB

source /local10G/rfhorns/resources/anaconda2/bin/activate /local10G/rfhorns/resources/anaconda2/
python /local10G/rfhorns/Bcell/flu_highres/scripts/subclones/annotate_nodes_FayAndWusH.py $1
