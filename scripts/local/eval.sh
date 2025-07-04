#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor/GREAC

#TESTDIR=~/Desktop/datasets/test_voc/test
#TRAIN=~/Desktop/datasets/test_voc/train/kmers
#GROUPNAME=covid

TESTDIR=~/Desktop/datasets/sars_cov2/test
TRAIN=~/Desktop/datasets/sars_cov2/train/kmers
GROUPNAME=covid_2

# TESTDIR=~/Desktop/datasets/bees/test
# TRAIN=~/Desktop/datasets/bees/kmers
# GROUPNAME=bees

# TESTDIR=~/Desktop/datasets/denv_eq/test
# TRAIN=~/Desktop/datasets/denv_eq/train/kmers
# GROUPNAME=denv_eq

cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME -w $1 performance-evaluation --train-dir $TRAIN


