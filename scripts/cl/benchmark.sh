#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/genomic-extractor/GREAC
DATAHOME=/tmp2/felipe
DATASETS=/home/a61491/datasets

TESTDIR=$DATASETS/bees/data/test
TRAIN=$DATASETS/bees/data/train/kmers_9
GROUPNAME=bees

# TRAIN=$DATAHOME/sars_cov2/train/kmers
# TESTDIR=$DATAHOME/sars_cov2/test
# GROUPNAME=covid_2



cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME \
    -w $1 benchmark --train-dir $TRAIN --test-dir $TESTDIR -m $2 -o ./output

# cd $PROJECTHOME &&  julia --project src/GREAC.jl  --group-name $GROUPNAME -w 0.001 fit-parameters  --train-dir $TRAIN --test-dir $TESTDIR 

