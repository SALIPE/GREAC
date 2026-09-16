#!/bin/bash
# Sweep 3-D (k x window x threshold) por organismo.
# Resolve os caminhos do dataset a partir do GROUPNAME e delega para
# fit-parameters, que varre as três dimensões num único processo Julia.
# Saída: GREAC/output-sweep-<GROUPNAME>/parameter_optimization_<data>.csv
#        GREAC/output-sweep-<GROUPNAME>/best_parameters_<data>.csv
set -u

DATASETS=~/Desktop/datasets/original

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
WINDOW_MAX=${4:-0.003}
WINDOW_STEP=${5:-0.0005}
THRESHOLD_MIN=${6:-0.5}
THRESHOLD_MAX=${7:-0.8}
THRESHOLD_STEP=${8:-0.05}
REPS=10
case $GROUPNAME in
    denv)      SOURCE=$DATASETS/denv/data_clean ;;
    hbv)       SOURCE=$DATASETS/HBV/data_clean ;;
    hiv)       SOURCE=$DATASETS/hiv/data_clean ;;
    sars)      SOURCE=$DATASETS/sars_cov2/data_clean ;;
    monkeypox) SOURCE=$DATASETS/mkpx/data_clean ;;
    bees[0-9]*)
        chr=${GROUPNAME#bees}
        if (( chr >= 1 && chr <= 16 )); then
            SOURCE=$DATASETS/bees/data_$chr
        else
            echo "❌ Erro: Número fora do intervalo permitido (1–16): $chr"
            exit 1
        fi
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
TESTDIR=$SOURCE/test

echo "📋 Sweep 3-D configurado:"
echo "   - GROUPNAME: $GROUPNAME"
echo "   - SOURCE:    $SOURCE"
echo "   - K_LIST:    $K_LIST"
echo "   - WINDOW:    $WINDOW -> $WINDOW_MAX (passo $WINDOW_STEP)"
echo "   - THRESHOLD: $THRESHOLD_MIN -> $THRESHOLD_MAX (passo $THRESHOLD_STEP)"
echo "   - REPS:      $REPS"

# Arquiva o CSV do sweep anterior: as repetições acumulam no mesmo arquivo,
# então sem isso o summary misturaria experimentos diferentes.
SWEEP_CSV=~/Desktop/GREAC/GREAC/output-sweep-$GROUPNAME/parameter_sweep_$GROUPNAME.csv
if [ -f "$SWEEP_CSV" ]; then
    mv "$SWEEP_CSV" "${SWEEP_CSV%.csv}_$(date +%Y%m%d%H%M%S).csv"
    echo "📦 Sweep anterior arquivado"
fi

for REP in $(seq 1 $REPS); do
    echo "=== Repetição $REP de $REPS"

    # novo split a cada repetição: é daqui que vem a variância entre as reps.
    # Dentro de uma repetição o split é fixo, então as combinações continuam
    # comparáveis entre si.
    if [ -x "$BALANCEDATASET/test.sh" ]; then
        $BALANCEDATASET/test.sh $SOURCE
    else
        echo "⚠️  $BALANCEDATASET/test.sh não encontrado, usando o split atual"
    fi

    $FIT "$TRAIN" "$TESTDIR" "$GROUPNAME" "$WINDOW" "$K_LIST" \
        "$WINDOW_MAX" "$WINDOW_STEP" "$THRESHOLD_MIN" "$THRESHOLD_MAX" "$THRESHOLD_STEP" \
        "$REP"
done

echo "Resultados: GREAC/output-sweep-$GROUPNAME/"
