#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC


GROUPNAME=$1
WINDOW=$2
TRAIN=$3
THRESHOLD=$4
REFERENCE=$5

cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME \
        -w $WINDOW \
        extract-features --train-dir $TRAIN 
