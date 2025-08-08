#!/bin/bash
source ~/.py-venv/bin/activate

export RUST_BACKTRACE=full

DENGUE=~/Desktop/datasets/dengue
HBV=~/Desktop/datasets/HBV/data
SARS=~/Desktop/datasets/sars_cov2
MONKEYPOX=~/Desktop/datasets/monkeypox
HIV=~/Desktop/GREAC/study-cases/castor_hiv_data/variants

BEES_1=~/Desktop/datasets/bees/data_1
BEES_2=~/Desktop/datasets/bees/data_2
BEES_3=~/Desktop/datasets/bees/data_3
BEES_4=~/Desktop/datasets/bees/data_4
BEES_5=~/Desktop/datasets/bees/data_5
BEES_6=~/Desktop/datasets/bees/data_6
BEES_7=~/Desktop/datasets/bees/data_7
BEES_8=~/Desktop/datasets/bees/data_8
BEES_9=~/Desktop/datasets/bees/data_9
BEES_10=~/Desktop/datasets/bees/data_10
BEES_11=~/Desktop/datasets/bees/data_11
BEES_12=~/Desktop/datasets/bees/data_12
BEES_13=~/Desktop/datasets/bees/data_13
BEES_14=~/Desktop/datasets/bees/data_14
BEES_15=~/Desktop/datasets/bees/data_15
BEES_16=~/Desktop/datasets/bees/data_16

GREAC=~/Desktop/GREAC/scripts/local/benchmark.sh
BALANCEDATASET=~/Desktop/Fasta-splitter/FastaSplitter

REF_HIV=~/Desktop/GREAC/study-cases/castor_hiv_data/hiv1_refseq.fasta
REF_HBV=~/Desktop/datasets/HBV/refseq.fasta
REF_DENV=~/Desktop/datasets/denv/refseq.fasta
REF_SARS=~/Desktop/datasets/tutorial_data/reference/SARS-CoV2_wuhan_refseq.fasta
REF_MONKEYPOX=~/Desktop/datasets/monkeypox-raw/refseq.fasta

REF_BEES_1=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group1.fasta
REF_BEES_2=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group2.fasta
REF_BEES_3=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group3.fasta
REF_BEES_4=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group4.fasta
REF_BEES_5=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group5.fasta
REF_BEES_6=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group6.fasta
REF_BEES_7=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group7.fasta
REF_BEES_8=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group8.fasta
REF_BEES_9=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group9.fasta
REF_BEES_10=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group10.fasta
REF_BEES_11=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group11.fasta
REF_BEES_12=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group12.fasta
REF_BEES_13=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group13.fasta
REF_BEES_14=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group14.fasta
REF_BEES_15=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group15.fasta
REF_BEES_16=~/Desktop/datasets/bees/GCA_000002195.1_Amel_4.5_genomic_Group16.fasta

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
    bees[0-9]*)
        if [[ $GROUPNAME =~ ^bees([0-9]+)$ ]]; then
            chr="${BASH_REMATCH[1]}"
            if (( chr >= 1 && chr <= 16 )); then

                REF_BEES="REF_BEES_${chr}"
                REF_TOTAL=$(eval echo \$$REF_BEES)

                SOURCE_VAR="BEES_${chr}"
                SOURCE=$(eval echo \$$SOURCE_VAR)
                echo "✅ Dataset BEES selecionado: $SOURCE"
            else
                echo "❌ Erro: Número fora do intervalo permitido (1–16): $chr"
                exit 1
            fi
        else
            echo "❌ Erro: Formato inválido para bees: $GROUPNAME"
            exit 1
        fi
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

function get_kmers_bees_chr() {
    local chr=$1  # valor numérico recebido como argumento

    # Construção dinâmica das variáveis M_Group e C_Group
    local m_group_var="M_Group${chr}"
    local c_group_var="C_Group${chr}"
    local ref_bees_var="REF_BEES_${chr}"

    # Acessando os valores com eval
    local m_group_val=$(eval echo \$$m_group_var)
    local c_group_val=$(eval echo \$$c_group_var)
    local ref_bees_val=$(eval echo \$$ref_bees_var)

    for variant in $m_group_val $c_group_val; do
        gramep get-only-kmers \
            --rpath "$ref_bees_val" \
            --spath "$SOURCE/train/${variant}.fasta" \
            --save-path "$SOURCE/train/kmers/" \
            --word "$KMERSIZE" \
            --step 1 -d ALL

        mv "$SOURCE/train/${variant}.fasta" "kmers/${variant}/${variant}.fasta"
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


for i in {1..1}; do
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
        hiv)
            get_kmers_hiv
            ;;
        sars)
            get_kmers_sars
            ;;
        monkeypox)
            # get_kmers_monkeypox
            ;;
        bees[0-9]*)
            if [[ $GROUPNAME =~ ^bees([0-9]+)$ ]]; then
                chr="${BASH_REMATCH[1]}"
                get_kmers_bees_chr $chr
            else
                echo "❌ Erro: Formato inválido para bees: $GROUPNAME"
                exit 1
            fi
            ;;
        *)
            echo "❌ Erro: GROUPNAME inválido: $GROUPNAME"
            exit 1
            ;;
    esac

    
    $GREAC $TRAIN $TESTDIR $GROUPNAME $WINDOW $METRIC $KMERSIZE $THRESHOLD $REF_TOTAL 
    
done

