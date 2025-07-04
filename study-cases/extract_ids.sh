#!/bin/bash

# Uso: ./extrair_ids.sh arquivo.fasta

if [ $# -lt 1 ]; then
  echo "Uso: $0 <arquivo_fasta>"
  exit 1
fi

arquivo="$1"

## To extract sars
#grep '^>' "$arquivo" | sed 's/^>//' | cut -d'|' -f1

## To extract hiv
grep '^>' "$arquivo" | sed 's/^>//' 

