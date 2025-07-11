#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC

# INPUT=~/Desktop/datasets/sars_cov2/
# GROUPNAME=sars

INPUT=~/Desktop/datasets/dengue/
GROUPNAME=denv

# INPUT=~/Desktop/GREAC/comparison_scripts/castor_hiv_data/variants/train/kmers
# GROUPNAME=hiv

# INPUT=~/Desktop/datasets/HBV/data/train/kmers
# GROUPNAME=hbv

cd $PROJECTHOME && julia --project src/GREAC.jl \
   --group-name $GROUPNAME \
   -w $1 fasta-regions -i $INPUT 


