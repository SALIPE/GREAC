#!/bin/sh

PROJECTHOME=../GREAC

HBV=~/Desktop/datasets/HBV/data
REF_HBV=study-cases/HBV/refseq.fasta


GROUPNAME= hbv
WINDOW= 0.002
KMERSIZE= 7 
THRESHOLD= 0.55

TRAIN=$SOURCE/train/kmers
TESTDIR=$SOURCE/test
METRIC=manhattan

echo "📋 Parâmetros e caminhos configurados:"
echo "   - GROUPNAME: $GROUPNAME"
echo "   - WINDOW: $WINDOW"
echo "   - TRAIN: $TRAIN"
echo "   - TESTDIR: $TESTDIR"
echo "   - METRIC: $METRIC"


cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache \
    --group-name $GROUPNAME \
    -w $WINDOW benchmark \
    --train-dir $TRAIN \
    --test-dir $TESTDIR \
    --metric $METRIC \
    --threshold $THRESHOLD \
    -o ./output-$KMERSIZE \
    --reference $REF_TOTAL 


