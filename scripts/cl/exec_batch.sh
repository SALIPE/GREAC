#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECT=/home/a61491/rrm-genomic-extractor/scripts/cl/benchmark.sh


for metric in manhattan euclidian chisquared mahalanobis kld; do
    $PROJECT 0.002 $metric
    $PROJECT 0.004 $metric
    $PROJECT 0.006 $metric
    $PROJECT 0.008 $metric
    $PROJECT 0.01 $metric
done
