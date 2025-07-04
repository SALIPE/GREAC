#!/bin/bash

PROJECTHOME=/home/salipe/Desktop/GitHub/rrm-genomic-extractor
DIR_INPUT=/home/salipe/Desktop/GitHub/datasets/tutorial_data/VOCs
OUTPUT=$PROJECTHOME/output_convergence


rm -r $OUTPUT
mkdir -p $OUTPUT
julia --project Main.jl -d $DIR_INPUT -w $1  -o $OUTPUT/euclidian_consensus\_$1% --convergence-analysis > $OUTPUT/convergence-analysis\_$1% 



