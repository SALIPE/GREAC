#!/bin/bash

ORGANISM=$1

qsub -q node17.q -v ORGANISMS="$ORGANISM",REPEATS=50 ./kmer_sweep.sh