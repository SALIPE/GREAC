#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC


TRAIN=$1
GROUPNAME=$2
WINDOW=$3
METRIC=$4
KMER=$5
THRESHOLD=$6


cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache \
        --group-name $GROUPNAME \
        -w $WINDOW \
        extract-features \
        --train-dir $TRAIN \
        -k $KMER \
        --threshold $THRESHOLD
