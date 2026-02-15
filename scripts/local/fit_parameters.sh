#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC

TRAIN=$1
TESTDIR=$2
GROUPNAME=$3
WINDOW=$4
KMER=$5


cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME \
   -w $WINDOW fit-parameters  \
   --train-dir $TRAIN \
   --test-dir $TESTDIR \
   --k-len $KMER \
   --classifier


