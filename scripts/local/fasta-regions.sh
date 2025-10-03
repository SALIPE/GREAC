#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC

# INPUT=~/Desktop/datasets/sars_cov2/
# GROUPNAME=sars

# INPUT=~/Desktop/datasets/dengue/
# GROUPNAME=denv

# INPUT=~/Desktop/datasets/mkpx/data
# GROUPNAME=monkeypox

# INPUT=~/Desktop/GREAC/study-cases/castor_hiv_data/variants
# GROUPNAME=hiv

# INPUT=~/Desktop/GREAC/study-cases/HBV/data
# GROUPNAME=hbv

GROUPNAME=$1
WINDOW=$2

cd $PROJECTHOME && julia --project src/GREAC.jl \
   --group-name $GROUPNAME \
   -w $WINDOW fasta-regions # -i $INPUT 


