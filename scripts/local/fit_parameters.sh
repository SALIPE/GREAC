#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC

TRAIN=$1
TESTDIR=$2
GROUPNAME=$3
WINDOW=$4
REFERENCE=$5


cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME \
   -w $WINDOW fit-parameters  \
   --train-dir $TRAIN \
   --test-dir $TESTDIR \
   --reference $REFERENCE \
   --usexgboost


