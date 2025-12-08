#!/bin/sh

PROJECTHOME=../GREAC

TUTORIAL_DATA=../study-cases/HBV/HBV_tutorial.tar.gz

rm -r HBV_tutorial
mkdir HBV_tutorial
tar -xvf $TUTORIAL_DATA -C HBV_tutorial

REF_HBV=../study-cases/HBV/refseq.fasta


GROUPNAME=hbv
WINDOW=0.002
KMERSIZE=7 
THRESHOLD=0.7



TRAIN=../scripts/HBV_tutorial/train/kmers
TESTDIR=../scripts/HBV_tutorial/test
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
    --k-len 7 \
    --threshold $THRESHOLD \
    -o ./output-$KMERSIZE \
    --reference $REF_HBV 


