#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc


PROJECTHOME=/home/a61491/rrm-genomic-extractor
DIR_INPUT=/home/a61491/datasets/tutorial_data/VOCs
# DIR_INPUT=/home/a61491/datasets/consensus
OUTPUT=$PROJECTHOME/output_convergence


cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1  -o out --convergence-analysis > $OUTPUT/convergence-class-analysis\_$1.txt
# cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1  -o out --convergence-analysis > $OUTPUT/convergence-organims-analysis\_$1.txt



