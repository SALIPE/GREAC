#!/bin/bash
# Sweep completo (k x window x threshold) via fit-parameters.
# Relatórios em GREAC/output-sweep-<GROUPNAME>/
set -u

PROJECTHOME=~/Desktop/GREAC/GREAC

TRAIN=$1
TESTDIR=$2
GROUPNAME=$3
WINDOW=$4
K_LIST=${5:-"5 6 7 8 9 10 11 12"}
WINDOW_MAX=${6:-0.003}
WINDOW_STEP=${7:-0.0005}
THRESHOLD_MIN=${8:-0.5}
THRESHOLD_MAX=${9:-0.8}
THRESHOLD_STEP=${10:-0.05}
REP=${11:-1}

cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME \
   -w $WINDOW fit-parameters \
   --train-dir $TRAIN \
   --test-dir $TESTDIR \
   -k $K_LIST \
   --window-max $WINDOW_MAX \
   --window-step $WINDOW_STEP \
   --threshold-min $THRESHOLD_MIN \
   --threshold-max $THRESHOLD_MAX \
   --threshold-step $THRESHOLD_STEP \
   --rep $REP \
   --classifier
