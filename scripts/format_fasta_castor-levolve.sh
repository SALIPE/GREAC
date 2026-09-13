#!/usr/bin/env sh

# Usage: ./format_fasta_castor-kevolve.sh <VARIANTE> <INPUT.fasta> [<OUTPUT.fasta>]

if [ "$#" -lt 2 ]; then
  echo "Uso: $0 <VARIANTE> <INPUT.fasta> [<OUTPUT.fasta>]"
  exit 1
fi

VARIANT="$1"
INPUT="$2"

if [ -n "$3" ]; then
  OUTPUT="$3"
else
  OUTPUT="formatted_${INPUT}"
fi

# Processa o FASTA: para cada linha de cabeçalho (que começa com '>'),
# substitui ">id" por ">id|VARIANT". Linhas de sequência permanecem inalteradas.
#sed "s/^>\(.*\)/>\1|${VARIANT}/" "$INPUT" > "$OUTPUT"
sed "s/^>\([^[:space:]]*\).*/>\1|${VARIANT}/" "$INPUT" > "$OUTPUT"

echo "Arquivo formatado salvo em: $OUTPUT"