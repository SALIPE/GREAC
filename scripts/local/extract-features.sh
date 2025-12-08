#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC


# TRAIN=$1
# GROUPNAME=$2
# WINDOW=$3
# METRIC=$4
# THRESHOLD=$5
# REFERENCE=$6

TRAIN=/home/salipe/Desktop/datasets/mkpx/data/train/kmers
GROUPNAME=mpox
WINDOW=0.0001
METRIC=manhattan
KMER=9
THRESHOLD=0.6
REFERENCE=~/Desktop/datasets/refseq_mkpx/ncbi_dataset/data/GCF_000857045.1/GCF_000857045.1_ViralProj15142_genomic.fna


cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache \
        --group-name $GROUPNAME \
        -w $WINDOW \
        extract-features \
        --train-dir $TRAIN \
        -k $KMER \
        --threshold $THRESHOLD \
        --reference $REFERENCE 
