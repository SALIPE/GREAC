#!/bin/bash
# Sweep 3-D (k x window x threshold) por organismo.
# Resolve os caminhos do dataset a partir do GROUPNAME e delega para
# fit-parameters, que varre as três dimensões num único processo Julia.
# Saída: GREAC/output-sweep-<GROUPNAME>/parameter_sweep_<GROUPNAME>.csv
#        GREAC/output-sweep-<GROUPNAME>/parameter_summary_<GROUPNAME>.csv
#
# Repetições pelo ambiente:
#   REPS=10 $0 hbv 0.002                      # sweep novo, arquiva o CSV anterior
#   REPEAT_START=42 REPS=10 $0 hiv 0.002      # continua na repetição 42
#   REPEAT_START=auto REPS=10 $0 hiv 0.002    # continua de onde o CSV parou
# Qualquer REPEAT_START diferente de 1 acrescenta ao CSV, e é recusado se a
# repetição já estiver gravada: repeti-la com outro split substituiria os
# resultados antigos. Uma repetição cortada no meio fica com o que terminou; o
# auto segue para a próxima em vez de misturar dois splits no mesmo número.
set -u

DATASETS=~/Desktop/datasets/new_sars_denv/denv_sars_cov_2

FIT=$(dirname "$0")/fit_parameters.sh
BALANCEDATASET=~/Desktop/Fasta-splitter/FastaSplitter

if [ $# -lt 2 ]; then
    echo "❌ Erro: Argumentos insuficientes"
    echo "Uso: $0 <GROUPNAME> <WINDOW> [K_LIST] [WMAX] [WSTEP] [TMIN] [TMAX] [TSTEP]"
    echo "Ex:  $0 hbv 0.002 '5 6 7 8 9 10' "
    exit 1
fi

GROUPNAME=$1
WINDOW=$2
K_LIST=${3:-"5 6 7 8 9 10"}
WINDOW_MAX=${4:-0.002}
WINDOW_STEP=${5:-0.0005}
THRESHOLD_MIN=${6:-0.5}
THRESHOLD_MAX=${7:-0.8}
THRESHOLD_STEP=${8:-0.05}
REPS=${REPS:-49}
REPEAT_START=${REPEAT_START:-1}   # 1 = sweep novo | N = continua em N | auto

if [ "$REPEAT_START" != "auto" ] && ! [[ "$REPEAT_START" =~ ^[1-9][0-9]*$ ]]; then
    echo "❌ Erro: REPEAT_START deve ser um inteiro >= 1 ou 'auto', recebido: $REPEAT_START"
    exit 1
fi

case $GROUPNAME in
    denv)      
    SOURCE=$DATASETS/denv/data_clean
    SOURCE_TESTE=$DATASETS/denv_ncbi
    ;;
    sars)      
    SOURCE=$DATASETS/sars_cov2/data_clean 
    SOURCE_TESTE=$DATASETS/sars_ncbi
    ;;
    *)
        echo "❌ Erro: GROUPNAME desconhecido: $GROUPNAME"
        exit 1
        ;;
esac

if [ ! -d "$SOURCE" ]; then
    echo "❌ Erro: Diretório SOURCE não existe: $SOURCE"
    exit 1
fi

TRAIN=$SOURCE/train
TESTDIR=$SOURCE_TESTE/test

SWEEP_CSV=~/Desktop/GREAC/GREAC/output-sweep-$GROUPNAME/parameter_sweep_$GROUPNAME.csv

# Maior repetição já gravada (0 se não houver CSV)
LAST=0
if [ -s "$SWEEP_CSV" ]; then
    LAST=$(awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "rep") c = i; next }
                    c && $c + 0 > m { m = $c + 0 } END { print m + 0 }' "$SWEEP_CSV")
fi

case "$REPEAT_START" in
    1)    START=1 ;;
    auto) START=$((LAST + 1)) ;;
    *)
        if (( REPEAT_START <= LAST )); then
            echo "❌ Erro: REPEAT_START=$REPEAT_START, mas o CSV já tem até a repetição $LAST."
            echo "   Use REPEAT_START=$((LAST + 1)) ou REPEAT_START=auto."
            exit 1
        fi
        START=$REPEAT_START
        ;;
esac
END=$((START + REPS - 1))

echo "📋 Sweep 3-D configurado:"
echo "   - GROUPNAME: $GROUPNAME"
echo "   - SOURCE:    $SOURCE"
echo "   - K_LIST:    $K_LIST"
echo "   - WINDOW:    $WINDOW -> $WINDOW_MAX (passo $WINDOW_STEP)"
echo "   - THRESHOLD: $THRESHOLD_MIN -> $THRESHOLD_MAX (passo $THRESHOLD_STEP)"
echo "   - REPS:      $START..$END"

# Só um sweep novo arquiva o CSV anterior; uma continuação acrescenta a ele.
if [ "$START" = "1" ] && [ -f "$SWEEP_CSV" ]; then
    mv "$SWEEP_CSV" "${SWEEP_CSV%.csv}_$(date +%Y%m%d%H%M%S).csv"
    echo "📦 Sweep anterior (reps 1..$LAST) arquivado"
fi

for REP in $(seq "$START" "$END"); do
    echo "=== Repetição $REP ($START..$END)"

    if [ -x "$BALANCEDATASET/run_maxtrain.sh" ]; then
        $BALANCEDATASET/run_maxtrain.sh $SOURCE
    else
        echo "⚠️  $BALANCEDATASET/run_maxtrain.sh não encontrado, usando o split atual"
    fi

    $FIT "$TRAIN" "$TESTDIR" "$GROUPNAME" "$WINDOW" "$K_LIST" \
        "$WINDOW_MAX" "$WINDOW_STEP" "$THRESHOLD_MIN" "$THRESHOLD_MAX" "$THRESHOLD_STEP" \
        "$REP"
done

echo "Resultados: GREAC/output-sweep-$GROUPNAME/"
