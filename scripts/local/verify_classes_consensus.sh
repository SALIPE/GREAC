#!/bin/bash

PROJECTHOME=/home/salipe/Desktop/GitHub/rrm-genomic-extractor
# DIR_INPUT=/home/salipe/Desktop/GitHub/datasets/tutorial_data/VOCs
DIR_INPUT=/home/salipe/Desktop/GitHub/datasets/tutorial_data/consensus_classes
OUTPUT=$PROJECTHOME/output_dist/$1
POSDIR=$PROJECTHOME/positions/all


rm -r $OUTPUT
mkdir -p $OUTPUT

cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1 -p $POSDIR -o $OUTPUT



