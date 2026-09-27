#!/bin/bash

echo "Trigger $1"

qsub -q node19.q -v REPEATS=100 ./kmer_sweep.sh $1

qsub -q node19.q -v REPEAT_START=auto,REPEATS=100 ./batch_benchmark.sh sars 0.001 10 0.5 