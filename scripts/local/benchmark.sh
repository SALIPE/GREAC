#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC

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
   -o ./output-$KMER \
   --classifier
