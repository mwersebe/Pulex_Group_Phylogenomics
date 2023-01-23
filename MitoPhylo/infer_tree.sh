#!/bin/bash

# Activate the conda environment with the needed tools:

source ~/miniconda3/etc/profile.d/conda.sh

conda activate BioInfoTools

# Clone a repo for concatentating sequences easily:

git clone https://github.com/marekborowiec/AMAS.git

cp AMAS/amas/AMAS.py .

# Concatentate the Genes and write a partition file:

./AMAS.py concat -i *.fa -f fasta -d dna -p Mito_Genes.txt -y raxml -u fasta -t Mito_Genes.fasta -c 8 -n 123

# Run IQtree to, infer the best model and partitioning scheme and then infer the ML tree:

iqtree -s Mito_Genes.fasta -p Mito_Genes.txt -m TESTMERGE -o Daphnia_retrocurvaClay -B 1000 -nt AUTO

# Convert the ML tree file to newick for reading into R.

mv Mito_Genes.txt.treefile Mito_Genes.nwk

## Done
