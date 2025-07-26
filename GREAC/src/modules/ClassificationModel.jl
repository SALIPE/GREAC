module ClassificationModel

include("RegionExtraction.jl")

using FLoops, .RegionExtraction, LinearAlgebra, Statistics, StatsBase, XGBoost
export ClassificationModel

struct MultiClassModel
    classes::Vector{String}
    class_string_probs::Dict{String,Vector{Float64}}
    variant_stats::Dict{String,Dict{Symbol,Float64}}
    kmerset::Set{String}
    regions::Vector{Tuple{Int,Int}}
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
    for (class, _) in meta_data

        class_seqs::Vector{Base.CodeUnits} = byte_seqs[class]
        println("Calculating $class probabilities")
        get_class_appearences = Base.Fix1(def_kmer_classes_probs, (regions, class_seqs))

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
        @floop for s in eachindex(class_seqs)
            seq::Base.CodeUnits = class_seqs[s]

            seq_distribution = sequence_kmer_distribution_optimized(regions, seq, collect(kmerset)) ./ length(kmerset)

            #manhttan distance for interval trust
            d = sum(abs.(seq_distribution - class_freq))
            in_group[s] = d

            if (use_xg)
                diverg::Vector{Float64} = zeros(length(seq_distribution) - 1)
                @inbounds for i in 1:(length(seq_distribution)-1)
                    diverg[i] = abs((seq_distribution[i+1] - seq_distribution[i]))
                end
                intern_X[s] = vcat(
                    [d, 0],
                    seq_distribution,
                    diverg
                )
            end
        end

        stats = Dict(
            :mu => mean(in_group),
            :sigma => std(in_group)
        )

        if (use_xg)
            for x_seq in intern_X
                x_seq[2] = gaussian_membership(stats, x_seq[1])
            end

            append!(X, intern_X)
            append!(y_str, intern_y)
        end

        variant_stats[class] = stats
    end


    if (use_xg)
        labels = unique(y_str)
        label2int = Dict(label => i - 1 for (i, label) in enumerate(labels))

        y = [label2int[cls] for cls in y_str]
        X_mat = convert(Matrix{Float32}, hcat(X...)')
        dtrain = DMatrix(X_mat, label=y)

        model = xgboost(dtrain,
            # num_round=20,
            # max_depth=30,
            # η=0.5,
            num_class=length(labels),
            objective="multi:softprob")
        # objective="multi:softmax")

        XGBoost.save(model, xg_model_name)
    end

    return MultiClassModel(
        [class for (class, _) in meta_data],
        class_string_probs,
        variant_stats,
        kmerset,
        regions)
end

function gaussian_membership(
    stats::Dict{Symbol,Float64},
    d::Float64
)
    mean = stats[:mu]
    std = stats[:sigma]
    return exp(-((d - mean)^2) / (2 * std^2))
end

function predict_membership(
    parameters::Tuple{MultiClassModel,Union{Nothing,String},String},
    X::Vector{Float64})::Tuple{String,Dict{String,Float64}}


    use_xg::Bool = false
    model, metric, xg_model_name = parameters

    classification = Dict{String,Float64}()

    if (use_xg)
        modelo_carregado = Booster(DMatrix[])
        XGBoost.load!(modelo_carregado, xg_model_name)

        diverg::Vector{Float64} = zeros(length(X) - 1)
        @inbounds for i in 1:(length(X)-1)
            diverg[i] = abs((X[i+1] - X[i]))
        end

    end

    for i in eachindex(model.classes)
        c = model.classes[i]
        class_freq = model.class_string_probs[c]
        stats = model.variant_stats[c]
        d = metrics_options(model, metric, class_freq, X)
        memb::Float64 = gaussian_membership(stats, d)

        if (use_xg)
            i_novo = [vcat(
                [d, memb],
                X,
                diverg
            )]
            i_novo_mat = convert(Matrix{Float32}, hcat(i_novo...)')
            y_pred_int = XGBoost.predict(modelo_carregado, DMatrix(i_novo_mat))

            classification[c] = y_pred_int[i]
        else
            classification[c] = (memb * d) + (1 - d)
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
)::Vector{UInt64}

    # Pré-processa os kmers para busca mais eficiente
    kmer_set = Set(codeunits(kmer) for kmer in kmerset)
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
            end
        end

        kmer_distribution[region_idx] = length(region_kmers)
    end

    return kmer_distribution
end

function sequence_kmer_distribution(
    regions::Vector{Tuple{Int,Int}},
    seq::Base.CodeUnits,
    kmerset::Vector{String}
)::Vector{UInt64}

    @floop for kmer in kmerset
        local_seq_histogram = zeros(UInt64, length(regions))
        seq_len = length(seq)

        for i in eachindex(regions)
            init_pos, end_pos = regions[i]

            if (end_pos > seq_len)
                end_pos = seq_len
            end

            wndw_buffer = @view seq[init_pos:end_pos]

            if RegionExtraction.occursinKmerBit(codeunits(kmer), wndw_buffer)
                local_seq_histogram[i] += 1
            end
        end

        @reduce(
            kmer_distribution = zeros(UInt64, length(regions)) .+ local_seq_histogram
        )
    end
    return kmer_distribution
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

function kld(
    class_freqs,
    X::Vector{Float64}
)
    #  Kullback-Leibler (KL) divergence
    Q_norm = X ./ sum(X)
    P_norm = class_freqs ./ sum(class_freqs)

    # Smooth to avoid zeros
    P_smoothed = P_norm .+ 1e-6
    P_smoothed = P_smoothed ./ sum(P_smoothed)

    # Compute KL(Q || P_smoothed)
    return sum(q * (log(q) - log(p)) for (q, p) in zip(Q_norm, P_smoothed) if q > 0)
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

    elseif metric == "kld"
        return kld(class_freqs, X)
    else
        error("Unsupported metric: $metric")
    end


end


end