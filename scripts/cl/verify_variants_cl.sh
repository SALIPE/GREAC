#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/rrm-genomic-extractor
VARIANT=Gamma
POSDIR=$PROJECTHOME/positions/$VARIANT
CONSENSUS_INPUT=/home/a61491/datasets/tutorial_data/VOCs
# CONSENSUS_INPUT=/home/a61491/datasets/consensus
WND=0.01
OUTPUT=$PROJECTHOME/output_convergence


for arquivo in "$POSDIR"/*; do
    # Verifica se é um arquivo regular
    if [[ -f "$arquivo" ]]; then
        echo "Processando arquivo: $arquivo"
        nome_arquivo=$(basename "$arquivo")
        cd $PROJECTHOME && julia --project Main.jl -f $CONSENSUS_INPUT/$VARIANT.fasta -w $WND -p $arquivo -o $OUTPUT/distance-validation-$nome_arquivo\_$WND --variant-name $VARIANT
        # julia --project Main.jl -d $CONSENSUS_INPUT -w 0.01 -p $arquivo -o $OUTPUT/$nome_arquivo\_10% --variant-name $VARIANT

    else
        echo "Ignorando: $arquivo (não é um arquivo regular)"
    fi
done


