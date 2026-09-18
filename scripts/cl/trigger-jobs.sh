#!/bin/bash

echo "Trigger $1"
qsub -q node17.q -v ORGANISMS=$1,REPEATS=50 ./kmer_sweep.sh