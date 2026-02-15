module ClassificationModel

include("RegionExtraction.jl")

using FLoops, .RegionExtraction, LinearAlgebra, Statistics, StatsBase, Distributions, DecisionTree, Serialization
export ClassificationModel

struct MultiClassModel
    classes::Vector{String}
    class_string_probs::Dict{String,Vector{Float64}}
    variant_stats::Dict{String,Dict{Symbol,Float64}}
    kmerset::Set{String}
    regions::Vector{Tuple{Int,Int}}
end

function calculate_mash_distance(jaccard_similarity::Float64, kmer_size::Int)::Float64
    if jaccard_similarity <= 0.0
        return 0 # A distância é infinita se não há similaridade
    end

    # A fórmula do Mash para a distância
    # Baseada no modelo de mutação de Poisson
    mash_dist = -log(2 * jaccard_similarity / (1 + jaccard_similarity)) / kmer_size

    return mash_dist
end

function calculate_mash_pvalue(jaccard_similarity::Float64, kmer_size::Int;
    sketch_size::Int=256,
    alphabet_size::Int=4)
    if !(0.0 <= jaccard_similarity <= 1.0)
        throw(ArgumentError("A similaridade de Jaccard deve estar entre 0 e 1."))
    end

    # 1. Estimar o número de hashes partilhados (x)
    # Deve ser um inteiro para os cálculos da distribuição
    x = round(Int, jaccard_similarity * sketch_size)

    # Se não há hashes em comum, a similaridade não é significativa (p-value = 1)
    if x == 0
        return 1.0
    end

    # 2. Calcular a probabilidade 'm' de um único hash aleatório ser partilhado
    # Usamos BigFloat para precisão numérica com k-mers grandes, pois r^k pode ficar enorme
    # A probabilidade de um k-mer aleatório específico ser selecionado num único sorteio de hash
    # é 1 / (alphabet_size^kmer_size)
    p_single_kmer = 1.0 / (BigFloat(alphabet_size)^kmer_size)

    # Probabilidade de um hash específico NÃO aparecer no outro esboço de tamanho 'sketch_size'
    prob_not_in_sketch = (1.0 - p_single_kmer)^sketch_size

    # Probabilidade 'm' de um hash de um esboço estar presente no outro esboço por acaso
    m = 1.0 - prob_not_in_sketch

    # 3. Criar a distribuição binomial para a hipótese nula
    # B(n, p) onde n=sketch_size e p=m
    dist = Binomial(sketch_size, Float64(m)) # Converte m de volta para Float64

    # 4. Calcular o p-value.
    # É a probabilidade de obter 'x' ou MAIS sucessos.
    # Isto corresponde à função de distribuição cumulativa complementar (ccdf) em x-1.
    # P(X >= x) = 1 - P(X <= x-1) = 1 - cdf(dist, x-1) = ccdf(dist, x-1)
    p_value = ccdf(dist, x - 1)

    return p_value
end

function fitMulticlass(
    kmerset::Set{String},
    meta_data::Dict{String,Int},
    byte_seqs::Dict{String,Vector{Base.CodeUnits}},
    regions::Vector{Tuple{Int,Int}},
    xg_model_name::String,
    use_xg::Bool=false
)::MultiClassModel

    class_string_probs = Dict{String,Vector{Float64}}()
    variant_stats = Dict{String,Dict{Symbol,Float64}}()
    regions_len = length(regions)

    X = Vector{Vector{Float64}}()
    y_str = String[]


    for class in keys(meta_data)

        class_seqs::Vector{Base.CodeUnits} = byte_seqs[class]
        println("Calculating $class probabilities")
        get_class_appearences = Base.Fix1(def_kmer_classes_probs, (regions, class_seqs))

        # Analyze the n behaviour based on all the data from the set
        @floop for kmer in collect(kmerset)
            kmer_seq_histogram = get_class_appearences(kmer)

            @reduce(
                kmer_distribution = zeros(UInt64, regions_len) .+ kmer_seq_histogram
            )
        end

        # Process F_wr
        class_freq = kmer_distribution ./ (length(kmerset) * length(class_seqs))
        class_string_probs[class] = class_freq

        in_group::Vector{Float64} = zeros(Float64, length(class_seqs))

        intern_X = Vector{Vector{Float64}}(undef, length(class_seqs))
        intern_y = fill(class, length(class_seqs))

        # Analyze each sequence in the variant set
        @floop for s in eachindex(class_seqs)
            seq::Base.CodeUnits = class_seqs[s]


            seq_distribution, kmer_counts = sequence_kmer_distribution_optimized(regions, seq, collect(kmerset))
            # seq_distribution = kmer_distribution ./ length(kmerset)

            diverg::Vector{Float64} = zeros(length(seq_distribution) - 1)
            @inbounds for i in 1:(length(seq_distribution)-1)
                diverg[i] = abs((seq_distribution[i+1] - seq_distribution[i]))
            end

            # manhattan distance for interval trust
            d = sum(abs.(seq_distribution - class_freq))
            in_group[s] = d

            if (use_xg)
                intern_X[s] = vcat(
                    [d, 0],
                    seq_distribution,
                    diverg,
                    kmer_counts
                )
            end
        end

        stats = Dict(
            :mu => mean(in_group),
            :sigma => std(in_group),
        )

        if (use_xg)
            for x_seq in intern_X
                x_seq[2] = gaussian_membership(stats[:mu], stats[:sigma], x_seq[1])
            end

            append!(X, intern_X)
            append!(y_str, intern_y)
        end

        variant_stats[class] = stats
    end

    if (use_xg)
        labels = unique(y_str)
        label2int = Dict(label => i for (i, label) in enumerate(labels))

        y = [label2int[cls] for cls in y_str]
        X_mat = convert(Matrix{Float32}, hcat(X...)')
        # dtrain = DMatrix(X_mat, label=y)

        model = DecisionTree.build_forest(y, X_mat)

        # Salvar modelo e mapeamento de labels
        model_data = Dict(
            "model" => model,
            "labels" => labels,
            "label2int" => label2int,
            "int2label" => Dict(v => k for (k, v) in label2int)
        )

        mkpath(dirname(xg_model_name))
        serialize(xg_model_name, model_data)

        # model = xgboost(
        #     dtrain,
        #     num_class=length(labels),
        #     objective="multi:softprob"
        # )
        # # objective="multi:softmax")

        # XGBoost.save(model, xg_model_name)
    end

    return MultiClassModel(
        [class for (class, _) in meta_data],
        class_string_probs,
        variant_stats,
        kmerset,
        regions)
end

function gaussian_membership(
    mean::Float64,
    std::Float64,
    d::Float64
)
    return exp(-((d - mean)^2) / (2 * std^2))
end


function predict_membership(
    parameters::Tuple{MultiClassModel,
        Union{Nothing,String},
        Bool,
        String},
    input::Tuple{Vector{Float64},Base.CodeUnits,Vector{Int32}})::Tuple{String,Dict{String,Float64}}

    model, metric, use_xg, xg_model_name = parameters
    X, seq, kmer_counts = input
    classification = Dict{String,Float64}()

    if (use_xg)
        # modelo_carregado = Booster(DMatrix[])
        # XGBoost.load!(modelo_carregado, xg_model_name)
        if !isfile(xg_model_name)
            error("Arquivo de modelo não encontrado: $xg_model_name")
        end

        modelo_carregado = deserialize(xg_model_name)
    end

    diverg::Vector{Float64} = zeros(length(X) - 1)
    @inbounds for i in 1:(length(X)-1)
        diverg[i] = abs((X[i+1] - X[i]))
    end

    for i in eachindex(model.classes)
        c = model.classes[i]
        class_freq = model.class_string_probs[c]
        stats = model.variant_stats[c]


        d = metrics_options(model, metric, class_freq, X)
        memb::Float64 = gaussian_membership(stats[:mu], stats[:sigma], d)

        if (use_xg)
            labels_int = collect(1:length(modelo_carregado["labels"]))
            i_novo = [vcat(
                [d, memb],
                X,
                diverg,
                kmer_counts
            )]
            i_novo_mat = convert(Matrix{Float32}, hcat(i_novo...)')
            # y_pred_int = XGBoost.predict(modelo_carregado, DMatrix(i_novo_mat))
            y_pred_proba = DecisionTree.apply_forest_proba(
                modelo_carregado["model"],
                i_novo_mat,
                labels_int)

            label_idx = findfirst(==(c), modelo_carregado["labels"])

            classification[c] = y_pred_proba[1, label_idx]
            # classification[c] = y_pred_int[i]
        else
            classification[c] = memb
        end

    end

    return argmax(classification), classification
end

function entropy(frequencias::Vector{Float64})::Float64
    entropia = -sum(p * log2(p) for p in frequencias if p > 0)
    return entropia
end


function def_kmer_classes_probs(
    seq_data::Tuple{Vector{Tuple{Int,Int}},Vector{Base.CodeUnits}},
    kmer::String)::Vector{UInt64}

    regions, sequences = seq_data

    regions_len = length(regions)
    fn_occursin = Base.Fix1(RegionExtraction.occursinKmerBit, codeunits(kmer))

    @floop for seq in sequences
        local_seq_histogram = zeros(UInt64, regions_len)
        seq_len = length(seq)

        for i in eachindex(regions)
            init_pos, end_pos = regions[i]

            if (end_pos > seq_len)
                end_pos = seq_len
            end

            wndw_buffer = @view seq[init_pos:end_pos]

            if fn_occursin(wndw_buffer)
                local_seq_histogram[i] += 1
            end
        end

        @reduce(
            seq_histogram = zeros(UInt64, regions_len) .+ local_seq_histogram
        )

    end

    return seq_histogram
end

function sequence_kmer_distribution_optimized(
    regions::Vector{Tuple{Int,Int}},
    seq::Base.CodeUnits,
    kmerset::Vector{String}
)::Tuple{Vector{Float64},Vector{Int32}}

    # Pre-process kmers for efficient search
    kmer_set = Set(codeunits(kmer) for kmer in kmerset)
    kmer_histogram = Dict{Vector{UInt8},Int32}((kmer, 0) for kmer in kmer_set)

    kmer_length = length(kmerset[1])

    kmer_distribution::Vector{UInt64} = zeros(UInt64, length(regions))
    seq_len = length(seq)

    @inbounds for (region_idx, (init_pos, end_pos)) in enumerate(regions)
        # Ajusta end_pos se necessário
        actual_end = min(end_pos, seq_len)

        # Extrai todos os kmers possíveis da região em uma única passada
        region_kmers = Set{Vector{UInt8}}()
        for i in init_pos:(actual_end-kmer_length+1)
            kmer_candidate = @view seq[i:(i+kmer_length-1)]
            kmer_vector = Vector{UInt8}(kmer_candidate)

            if kmer_vector in kmer_set
                push!(region_kmers, kmer_vector)
                kmer_histogram[kmer_vector] += 1

            end
        end

        kmer_distribution[region_idx] = length(region_kmers)
    end

    return kmer_distribution ./ length(kmerset), collect(values(kmer_histogram))
end


function predict_raw(
    parameters::Tuple{MultiClassModel,Union{Nothing,String}},
    X::Vector{Float64})::Tuple{String,Dict{String,Float64}}

    model, metric = parameters

    dists = Dict{String,Float64}([(class, zero(Float64)) for class in model.classes])
    for c in model.classes

        # Get the class's precomputed conditional frequencies
        class_freqs = model.class_string_probs[c]
        dists[c] = metrics_options(model, metric, class_freqs, X)
    end

    return argmin(dists), dists
end

function metrics_options(
    model,
    metric::Union{Nothing,String},
    class_freqs,
    X::Vector{Float64}
)

    if (isnothing(metric))
        metric = "manhattan"
    end
    epsilon = 1e-6
    if metric == "manhattan"

        # Manhattan distance
        return sum(abs.(X - class_freqs))

    elseif metric == "euclidian"
        # Euclidian distance
        return sqrt(sum((X - class_freqs) .^ 2))

    elseif metric == "mahalanobis"
        # Need repair, is not receivnig model here
        train_data = hcat([model.class_string_probs[c] for c in model.classes]...)

        covariance = cov(train_data; dims=2)
        inv_covariance = inv(covariance + epsilon * I(size(covariance, 1)))

        # Mahalanobis distance requires inverse covariance matrix
        if inv_covariance === nothing
            error("Mahalanobis metric requires inverse covariance matrix")
        end

        delta = X - class_freqs
        return sqrt(delta' * inv_covariance * delta)

    elseif metric == "chisquared"
        # Chi-squared distance
        return sum((X - class_freqs) .^ 2 ./ (class_freqs .+ 1e-9))

    else
        error("Unsupported metric: $metric")
    end


end


end