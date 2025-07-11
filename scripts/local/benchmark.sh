#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC


# TESTDIR=~/Desktop/datasets/sars_cov2/test
# TRAIN=~/Desktop/datasets/sars_cov2/train/kmers
# GROUPNAME=covid_2

#TESTDIR=~/Desktop/datasets/bees/data/test
#TRAIN=~/Desktop/datasets/bees/data/train/kmers_9
#GROUPNAME=bees

# TESTDIR=~/Desktop/datasets/dengue/test
# TRAIN=~/Desktop/datasets/dengue/train/kmers
# GROUPNAME=dengue

TRAIN=$1
TESTDIR=$2
GROUPNAME=$3
WINDOW=$4
METRIC=$5
KMER=$6
THRESHOLD=$7


# TESTDIR=~/Desktop/genomic-extractor/comparison_scripts/castor_hiv_data/variants/test
# TRAIN=~/Desktop/genomic-extractor/comparison_scripts/castor_hiv_data/variants/train/kmers
# GROUPNAME=hiv

# TESTDIR=~/Desktop/datasets/HBV/data/test
# TRAIN=~/Desktop/datasets/HBV/data/train/kmers
# GROUPNAME=hbv

cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME \
   -w $WINDOW benchmark --train-dir $TRAIN --test-dir $TESTDIR -m $METRIC --threshold $THRESHOLD \
   -o ./output-$KMER


# docker run --rm --cpus="4" -e JULIA_NUM_THREADS=4 \
#     -v $TRAIN:/train_dir \
#     -v $TESTDIR:/test_dir \
#     -v ./data/output:/output \
#     -v ~/.project_cache:/root/.project_cache \
#     greac:latest \
#     julia --project=/app /app/src/GREAC.jl \
#     --no-cache --group-name $GROUPNAME \
#     -w $WINDOW \
#     benchmark --train-dir /train_dir \
#     --test-dir /test_dir  \
#     -m $METRIC \
# # #     -o /output
# cd $PROJECTHOME &&  julia --project src/GREAC.jl --no-cache \
#    --group-name $GROUPNAME -w $WINDOW fit-parameters --train-dir $TRAIN --test-dir $TESTDIR 


