#!/bin/bash
source ~/.py-venv/bin/activate

export RUST_BACKTRACE=full

DENGUE=~/Desktop/datasets/dengue
HBV=~/Desktop/datasets/HBV/data
BEES=~/Desktop/datasets/bees/data
SARS=~/Desktop/datasets/sars_cov2
MONKEYPOX=~/Desktop/datasets/monkeypox
HIV=~/Desktop/GREAC/study-cases/castor_hiv_data/variants

GREAC=~/Desktop/GREAC/scripts/local/benchmark.sh
BALANCEDATASET=~/Desktop/Fasta-splitter/FastaSplitter

REF_HIV=~/Desktop/GREAC/study-cases/castor_hiv_data/hiv1_refseq.fasta
REF_HBV=~/Desktop/datasets/HBV/refseq.fasta
REF_DENV=~/Desktop/datasets/denv/refseq.fasta
REF_SARS=~/Desktop/datasets/tutorial_data/reference/SARS-CoV2_wuhan_refseq.fasta
REF_MONKEYPOX=~/Desktop/datasets/monkeypox-raw/refseq.fasta
REF_BEES=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group1.fasta

if [ $# -lt 4 ]; then
    echo "❌ Erro: Argumentos insuficientes"
    echo "Uso: $0 <GROUPNAME> <WINDOW> <KMERSIZE> <THRESHOLD>"
    exit 1
fi

GROUPNAME=$1
WINDOW=$2
KMERSIZE=$3 
THRESHOLD=$4

echo "📋 Parâmetros recebidos:"
echo "   - GROUPNAME: $GROUPNAME"
echo "   - WINDOW: $WINDOW"

case $GROUPNAME in
    denv)
        SOURCE=$DENGUE
        echo "✅ Dataset DENGUE selecionado: $SOURCE"
        ;;
    hbv)
        SOURCE=$HBV
        echo "✅ Dataset HBV selecionado: $SOURCE"
        ;;
    bees)
        SOURCE=$BEES
        echo "✅ Dataset BEES selecionado: $SOURCE"
        ;;
    hiv)
        SOURCE=$HIV
        echo "✅ Dataset HIV selecionado: $SOURCE"
        ;;
    sars)
        SOURCE=$SARS
        echo "✅ Dataset SARS selecionado: $SOURCE"
        ;;
    monkeypox)
        SOURCE=$MONKEYPOX
        echo "✅ Dataset MONKEYPOX selecionado: $SOURCE"
        ;;
    *)
        echo "❌ Erro: GROUPNAME deve ser 'denv' ou 'hbv'"
        exit 1
        ;;
esac


function get_kmers_denv() {
    
    for variant in type1 type2 type3 type4; do
        
        gramep get-only-kmers \
            --rpath $REF_DENV \
            --spath $SOURCE/train/$variant.fasta \
            --save-path $SOURCE/train/kmers/ \
            --word $KMERSIZE \
            --step 1
        
        mv $SOURCE/train/$variant.fasta kmers/$variant/$variant.fasta
    done
}

function get_kmers_hbv() {
    
    for variant in A B C D E F; do
        gramep get-only-kmers \
            --rpath $REF_HBV \
            --spath $SOURCE/train/$variant.fasta \
            --save-path $SOURCE/train/kmers/ \
            --word $KMERSIZE \
            --step 1 
        
        mv $SOURCE/train/$variant.fasta kmers/$variant/$variant.fasta
    done
}

function get_kmers_monkeypox() {
    
    for variant in A B C D E F; do
        gramep get-only-kmers \
            --rpath $REF_MONKEYPOX \
            --spath $SOURCE/train/$variant.fasta \
            --save-path $SOURCE/train/kmers/ \
            --word $KMERSIZE \
            --step 1 -d ALL
        
        mv $SOURCE/train/$variant.fasta kmers/$variant/$variant.fasta
    done
}

function get_kmers_hiv() {
    
    for variant in HIV1_A HIV1_B HIV1_C HIV1_D HIV1_F HIV1_G; do
        gramep get-only-kmers \
            --rpath $REF_HIV \
            --spath $SOURCE/train/$variant.fasta \
            --save-path $SOURCE/train/kmers/ \
            --word $KMERSIZE \
            --step 1
        
        mv $SOURCE/train/$variant.fasta kmers/$variant/$variant.fasta
    done
}

function get_kmers_sars() {
    
    for variant in Alpha Beta Delta Epsilon Eta Gamma Iota Kappa Lambda Omicron; do
        gramep get-only-kmers \
            --rpath $REF_SARS \
            --spath $SOURCE/train/$variant.fasta \
            --save-path $SOURCE/train/kmers/ \
            --word $KMERSIZE \
            --step 1
        
        mv $SOURCE/train/$variant.fasta kmers/$variant/$variant.fasta
    done
}


function get_kmers_bees() {
    
    for variant in M_Group1 C_Group1; do
        gramep get-only-kmers \
            --rpath $REF_BEES \
            --spath $SOURCE/train/$variant.fasta \
            --save-path $SOURCE/train/kmers/ \
            --word $KMERSIZE  \
            --step 1 -d ALL
        
        mv $SOURCE/train/$variant.fasta kmers/$variant/$variant.fasta
    done
}

TRAIN=$SOURCE/train/kmers
TESTDIR=$SOURCE/test
METRIC=manhattan

echo "📁 Caminhos configurados:"
echo "   - TRAIN: $TRAIN"
echo "   - TESTDIR: $TESTDIR"
echo "   - METRIC: $METRIC"

echo "🔍 Verificando existência dos diretórios..."
if [ ! -d "$SOURCE" ]; then
    echo "❌ Erro: Diretório SOURCE não existe: $SOURCE"
    exit 1
fi


for i in {1..90}; do
    echo "Iteração $i de 100"
    
    $BALANCEDATASET/test.sh $SOURCE

    cd $SOURCE/train
    mkdir -p kmers

     case $GROUPNAME in
        denv)
            get_kmers_denv
            ;;
        hbv)
            get_kmers_hbv
            ;;
        bees)
            get_kmers_bees
            ;;
        hiv)
            get_kmers_hiv
            ;;
        sars)
            get_kmers_sars
            ;;
        monkeypox)
            #get_kmers_monkeypox
            ;;
        *)
            echo "❌ Erro: GROUPNAME inválido: $GROUPNAME"
            exit 1
            ;;
    esac
    
    $GREAC $TRAIN $TESTDIR $GROUPNAME $WINDOW $METRIC $KMERSIZE $THRESHOLD
    
done

