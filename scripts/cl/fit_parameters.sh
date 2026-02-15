#!/bin/sh

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/GREAC/GREAC

TRAIN=$1
TESTDIR=$2
GROUPNAME=$3
WINDOW=$4
KMER=$5
REFERENCE=$6


cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME \
   -w $WINDOW fit-parameters  \
   --train-dir $TRAIN \
   --test-dir $TESTDIR \
   --k-len $KMER \
   --reference $REFERENCE \
   --usexgboost


