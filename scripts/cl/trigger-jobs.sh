#!/bin/bash

echo "Trigger $1"

qsub -q node19.q -v REPEATS=100 ./kmer_sweep.sh $1