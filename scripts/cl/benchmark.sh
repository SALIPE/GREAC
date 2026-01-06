#!/bin/sh

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/GREAC/GREAC
# DATAHOME=/tmp2/felipe
# DATASETS=/home/a61491/datasets

# TESTDIR=$DATASETS/bees/data/test
# TRAIN=$DATASETS/bees/data/train/kmers_9
# GROUPNAME=bees

TRAIN=$1
TESTDIR=$2
GROUPNAME=$3
WINDOW=$4
METRIC=$5
KMER=$6
THRESHOLD=$7


cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME \
   -w $WINDOW benchmark \
   --train-dir $TRAIN \
   --test-dir $TESTDIR \
   -m $METRIC \
   -k $KMER \
   --threshold $THRESHOLD \
   -o ./output-$KMER\
   --classifier


