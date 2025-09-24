#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC


TRAIN=$1
TESTDIR=$2
GROUPNAME=$3
WINDOW=$4
METRIC=$5
KMER=$6
THRESHOLD=$7
REFERENCE=$8


cd $PROJECTHOME && julia --project src/GREAC.jl --group-name $GROUPNAME \
   -w $WINDOW file-classification --file $TESTDIR -m $METRIC --threshold $THRESHOLD \
   -o ./output-$KMER --reference $REFERENCE



