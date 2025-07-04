#!/bin/bash

PROJECTHOME=/home/salipe/Desktop/GitHub/rrm-genomic-extractor
VARIANT=$1
POSDIR=$PROJECTHOME/positions/$VARIANT
CONSENSUS_INPUT=$2
OUTPUT=$PROJECTHOME/out/$VARIANT

rm -r $PROJECTHOME/out/$VARIANT
mkdir -p $PROJECTHOME/out/$VARIANT
# Loop pelos arquivos no diretório
for arquivo in "$POSDIR"/*; do
    # Verifica se é um arquivo regular
    if [[ -f "$arquivo" ]]; then
        echo "Processando arquivo: $arquivo"
        nome_arquivo=$(basename "$arquivo")
        julia --project Main.jl -f $CONSENSUS_INPUT -w 0.01 -p $arquivo -o $OUTPUT/$nome_arquivo\_10% --variant-name $VARIANT
        # julia --project Main.jl -d $CONSENSUS_INPUT -w 0.01 -p $arquivo -o $OUTPUT/$nome_arquivo\_10% --variant-name $VARIANT

    else
        echo "Ignorando: $arquivo (não é um arquivo regular)"
    fi
done


