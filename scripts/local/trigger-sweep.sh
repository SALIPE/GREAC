#!/bin/bash

# REPEAT_START=auto REPS=61 ./kmer_sweep.sh hbv 0.0015 
# REPEAT_START=auto REPS=100 ./kmer_sweep.sh mkpx 0.001 
# REPEAT_START=auto REPS=10 ./kmer_sweep.sh hiv 0.0015

REPEAT_START=auto REPS=100 ./batch_benchmark.sh denv 0.002 8 0.6 
#REPEAT_START=auto REPS=100 ./batch_benchmark.sh sars 0.001 10 0.5 