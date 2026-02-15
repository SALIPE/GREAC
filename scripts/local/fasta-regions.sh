#!/bin/bash

PROJECTHOME=~/Desktop/GREAC/GREAC


INPUT=$1
GROUPNAME=$2
WINDOW=$3

CMD="cd $PROJECTHOME && julia --project src/GREAC.jl \
   --group-name $GROUPNAME -w $WINDOW fasta-regions"

if [ -n "$INPUT" ]; then
    CMD="$CMD -i $INPUT"
fi

eval $CMD


