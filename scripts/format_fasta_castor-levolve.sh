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

# Processa o FASTA: cada linha de cabeçalho (que começa com '>') fica no formato
#     >id_CLASSE|CLASSE
# com um único '|' separando id e classe, porque é assim que o CASTOR-KRFE e o
# KEVOLVE leem o rótulo. No id não há espaços nem pipes: ambos viram '_',
# inclusive os pipes que a própria classe tenha ("SARS-CoV-2|Alpha" ->
# "SARS-CoV-2_Alpha"). A classe não é repetida quando o cabeçalho já termina
# nela, então reprocessar um arquivo já formatado não muda nada.
#   >accn PQ481563|SARS-CoV-2|Alpha  ->  >accn_PQ481563_SARS-CoV-2_Alpha|SARS-CoV-2|Alpha
#   >accn KC692509 type1|type1       ->  >accn_KC692509_type1|type1
# Linhas de sequência permanecem inalteradas.
awk -v variante="$VARIANT" '
  function normalizar(t) {
    gsub(/[|[:space:]]+/, "_", t)
    sub(/^_+/, "", t)
    sub(/_+$/, "", t)
    return t
  }
  BEGIN { vid = normalizar(variante) }     # a classe como aparece dentro do id
  /^>/ {
    id = normalizar(substr($0, 2))
    # tira TODAS as repetições da classe no fim do id, para acrescentar uma só
    while (id != "") {
      if (id == vid) { id = ""; break }
      corte = length(id) - length(vid) - 1
      if (corte > 0 && substr(id, corte + 1) == "_" vid) id = substr(id, 1, corte)
      else break
    }
    print ">" (id == "" ? "" : id "_") vid "|" variante
    next
  }
  { print }
' "$INPUT" > "$OUTPUT"

echo "Arquivo formatado salvo em: $OUTPUT"