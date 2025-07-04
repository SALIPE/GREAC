#!/bin/bash

PROJECTHOME=~/Desktop/genomic-extractor/GREAC



TRAIN=~/Desktop/datasets/sars_cov2/train/kmers
GROUPNAME=sars

#TESTDIR=~/Desktop/datasets/bees/data/test
#TRAIN=~/Desktop/datasets/bees/data/train/kmers_9
#GROUPNAME=bees

# TESTDIR=~/Desktop/datasets/dengue/test
# TRAIN=~/Desktop/datasets/dengue/train/kmers
# GROUPNAME=dengue

# TESTDIR=~/Desktop/genomic-extractor/comparison_scripts/castor_hiv_data/variants/test
# TRAIN=~/Desktop/genomic-extractor/comparison_scripts/castor_hiv_data/variants/train/kmers
# GROUPNAME=hiv

# TESTDIR=~/Desktop/datasets/HBV/data/test
# TRAIN=~/Desktop/datasets/HBV/data/train/kmers
# GROUPNAME=hbv

cd $PROJECTHOME && julia --project src/GREAC.jl \
   --group-name $GROUPNAME \
   -w $1 fasta-regions -i $TRAIN 


