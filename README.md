# BCellSelection
Code associated with Horns et al. (2019). This code was developed to analyze the population genetic signatures of selection in human B cells.

## About
This repository contains workflows to perform phylogenetic and population genetic analysis of B cell lineages and code to reproduce the figures in Horns et al. (2019).

## Configuration

### Environment

Most code for analysis uses a Python 2.7 environment with numerous modules installed. This environment can be created and activated by running the following. This workflow assumes you have [Anaconda](https://conda.io/miniconda.html) installed.

1. Create a Python 2.7 environment from the environment.yml file:

```bash
conda env create -f environment.yml
```

2. Activate the new environment:

`source activate HornsSelection2019`.

### Environment for Snakemake workflows

Several workflows use the [snakemake](https://snakemake.readthedocs.io/en/stable/) tool to deploy computation in a cluster-computing environment. These workflows require a Python 3.7 environment. This environment can be created and activated by running the following.

1. Create a Python 3.7 environment from the environment_py37.yml file:

```bash
conda env create -f environment_py37.yml
```

2. Activate the new environment:

`source activate HornsSelection2019_py37`.

### External dependencies

External tools are used for several scripts. This section lists these tools and provides links to download and install them.

1. [MUSCLE 3.8.31](https://www.drive5.com/muscle/downloads.htm) for multiple sequence alignment.

2. [FastTree 2.1.7](http://www.microbesonline.org/fasttree/#Install) for phylogenetic reconstruction.

3. [FitnessInference](https://github.com/rneher/FitnessInference) package for inference of fitness based on the shape of phylogenies.

### Data

Preprocessed data is available for [download](http://bit.ly/2BL83JV). These data are:

1. Sequences (`sequences.tsv`)

2. Sequence features or annotations (`sequence_annotations.tsv`)

3. Clonal lineage features (`lineage_dynamics.csv`)

These data are the input files for the code to reproduce the figures.

## Contents

### Figures

[Jupyter notebooks](https://jupyter.org/) are provided for reproducing figures. The title of each notebook describes which figures the code is associated with.

### Workflows in `scripts`

#### `trees`

In this workflow, multiple seqence alignment and phylogenetic reconstruction are performed on clonal B cell lineages. The workflow is initialized using a seedfile (`seedfile.txt`), which provides a list of lineage identifiers (corresponding to `lineage_uid` in `sequence_annotations.tsv`). This workflow is designed to be deployed in a cluster-computing setting (e.g. on a SLURM cluster). It can be launched using snakemake. An example of how to launch the workflow using snakemake on a SLURM cluster is given in `do.sh`.

#### `plotTrees`

In this workflow, phylogenies are rendered as PDFs. The input is a Newick file representing the phylogeny and a fasta file representing the multiple sequence alignment supporting the phylogeny.

#### `subclones`

In this workflow, clades within trees are annotated with population genetic signatures of selection as quantified using the Fay and Wu's H metric. The input is a list of lineage identifiers (corresponding to `lineage_uid` in `sequence_annotations.tsv`).

#### `fitnessMutations`

In this workflow, fitness inference is performed based on the shape of phylogenies. This workflow is designed to be deployed in a cluster-computing setting (e.g. on a SLURM cluster). It can be launched using snakemake. An example of how to launch the workflow using snakemake on a SLURM cluster is given in `do.sh`.

## Disclaimer

This project is not maintained. Software is provided as is and requests for support may not be addressed. 

## Contact

If you have questions or comments, please contact Felix Horns at <rfhorns@gmail.com>.
