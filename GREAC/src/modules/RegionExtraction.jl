module RegionExtraction

include("DataIO.jl")

using .DataIO,
    Serialization,
    AbstractFFTs,
    DSP,
    Statistics,
    Normalization,
    FASTX,
    FLoops

export RegionExtraction

function filterRegions(
    variantDirPath::String,
    wnwPercent::Float32,
    groupName::String,
    wndwSize::Int,
    maxSeqLen::Int
)::Vector{Tuple{Int,Int}}

    variantDirs::Vector{String} = readdir(variantDirPath)
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"

    hist_collection::Vector{Vector{UInt64}} = []

    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        cache::Union{Nothing,Tuple{String,Tuple{Vector{UInt16},BitArray}}} = DataIO.load_cache("$cachdir/$(variant)_outmask.dat")
        push!(hist_collection, cache[2][1])
    end

    min_size = minimum(length, hist_collection)


    desvio = zeros(Float64, min_size)

    @inbounds for i in 1:min_size
        values_at_point_i = [s[i] for s in hist_collection]
        desvio[i] = std(values_at_point_i)
    end

    marked = falses(min_size)

    @inbounds for i in eachindex(desvio)
        if (desvio[i] > mean(desvio))
            init = max(1, i - 2)
            endi = min(min_size, i + 2)
            while (init <= endi)
                end_pos = min(maxSeqLen, init + wndwSize - 1)
                marked[init:end_pos] .= true
                init += 1
            end
        end
    end

    extracted_regions = Vector{Tuple{Int,Int}}()

    init_pos = 1
    current::Bool = false

    for i in eachindex(marked)
        if (marked[i] && !current)
            current = true
            init_pos = i
        elseif (!marked[i] && current)
            push!(extracted_regions, (init_pos, i - 1))
            init_pos = i
            current = false
        end
    end

    if (current)
        push!(extracted_regions, (init_pos, length(marked)))
    end
    return extracted_regions
end


function regionsConjuction(
    variantDirPath::String,
    wnwPercent::Float32,
    groupName::String
)::Vector{Tuple{Int,Int}}

    variantDirs::Vector{String} = readdir(variantDirPath)
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"

    hit_region::Union{Nothing,BitArray} = nothing

    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        cache::Union{Nothing,Tuple{String,Tuple{Vector{UInt16},BitArray}}} = DataIO.load_cache("$cachdir/$(variant)_outmask.dat")

        marked = cache[2][2]

        if (isnothing(hit_region))
            hit_region = marked
        else
            for i in eachindex(marked)
                if (ismissing(hit_region[i]))
                    hit_region[i] = marked[i]
                elseif (marked[i] && !hit_region[i])
                    hit_region[i] = true
                end

            end
        end
    end

    extracted_regions = Vector{Tuple{Int,Int}}()

    init_pos = 1
    current::Bool = false

    for i in eachindex(hit_region)
        if (hit_region[i] && !current)
            current = true
            init_pos = i
        elseif (!hit_region[i] && current)
            push!(extracted_regions, (init_pos, i - 1))
            init_pos = i
            current = false
        end
    end

    if (current)
        push!(extracted_regions, (init_pos, length(hit_region)))
    end

    return extracted_regions

end

function generate_kmers(k::Int)
    nucleotides = ['A', 'C', 'T', 'G']
    sketch = Dict{String,Int}()
    function gerar_combinacoes(prefixo::String, tamanho_restante::Int)
        if tamanho_restante == 0
            sketch[prefixo] = 0
            return
        end

        for nucleotide in nucleotides
            gerar_combinacoes(prefixo * nucleotide, tamanho_restante - 1)
        end
    end

    gerar_combinacoes("", k)

    return sketch
end

function rolling_hash_kmers(
    sequence::String,
    kmer_hash_map::Dict{UInt64,Tuple{String,Int32}},
    k_len::Int)::Dict{UInt64,Tuple{String,Int32}}

    seq_len = length(sequence)

    update_hashmap = Dict{UInt64,Tuple{String,Int32}}()

    if seq_len < k_len
        error("K-mer size must be greater than the inputed sequence")
    end
    if length(keys(kmer_hash_map)) == 0
        error("0 k-mers found!")
    end

    base = UInt64(5)

    power = UInt64(1)
    for i in 1:(k_len-1)
        power *= base
    end


    # Have to init the hash
    current_hash = UInt64(0)
    @inbounds for i in 1:k_len
        current_hash = current_hash * base + UInt64(sequence[i])
    end

    if haskey(kmer_hash_map, current_hash)
        update_hashmap[current_hash] = kmer_hash_map[current_hash]
        kmer, value = update_hashmap[current_hash]
        update_hashmap[current_hash] = (kmer, value + 1)
    end

    @inbounds for i in (k_len+1):seq_len
        current_hash = current_hash - UInt64(sequence[i-k_len]) * power
        current_hash = current_hash * base + UInt64(sequence[i])


        if haskey(kmer_hash_map, current_hash)
            if haskey(update_hashmap, current_hash)
                kmer, value = update_hashmap[current_hash]
                update_hashmap[current_hash] = (kmer, value + 1)
            else
                update_hashmap[current_hash] = kmer_hash_map[current_hash]
                kmer, value = update_hashmap[current_hash]
                update_hashmap[current_hash] = (kmer, value + 1)
            end
        end
    end

    return update_hashmap

end

function max_entropy(kmers::Dict{String,Int32})
    # Extrair valores e ordenar
    data = collect(values(kmers))
    total = sum(data)
    sort!(data, rev=true)

    # Normalizar dados
    normalized_data = data ./ total


    function entropy(probs::Vector{Float64})
        -sum(p * log(p) for p in probs if p > 0.0)
    end

    # Calcular curva de entropia 
    entropy_curve = map(1:length(normalized_data)-1) do s
        # Região A: probs[1:s], probabilidade total P_A
        p_a = sum(normalized_data[1:s])
        h_a = if p_a > 0.0
            p_a_data = normalized_data[1:s] ./ p_a
            entropy(p_a_data)
        else
            0.0
        end

        # Região B: probs[s+1:end], probabilidade total P_B
        p_b = sum(normalized_data[s+1:end])
        h_b = if p_b > 0.0
            p_b_data = normalized_data[s+1:end] ./ p_b
            entropy(p_b_data)
        else
            0.0
        end

        h_a + h_b
    end

    # Encontrar índice de máxima entropia
    max_entropy_idx = argmax(entropy_curve)
    threshold = max_entropy_idx
    frequency = data[threshold]
    @show frequency

    return frequency
end

function get_exclusive_kmers(
    k_len::Int,
    variantDirPath::String)::Set{String}

    sketch::Dict{String,Int} = generate_kmers(k_len)

    kmer_hash_map = Dict{UInt64,Tuple{String,Int32}}() # 4^klen

    for kmer in keys(sketch)
        h = compute_hash(kmer)
        if !haskey(kmer_hash_map, h)
            kmer_hash_map[h] = (kmer, Int32(0))
        end
    end

    all_exclusive_kmers = Set{String}()

    variantDirs::Vector{String} = readdir(variantDirPath)
    variant_kmers = Dict{String,Set{String}}()


    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        println("Getting $variant k-mers")
        var_hash = kmer_hash_map

        sequences::Vector{String} = DataIO.loadStringSequences("$variantDirPath/$variant")
        for seq in sequences
            var_hash = rolling_hash_kmers(seq, var_hash, k_len)
        end

        # A busca das mais informativas tem que ser aqui em comparação com a referencia
        kmer_dict = Dict{String,Int32}()
        for kmer_freq in values(var_hash)
            kmer_dict[kmer_freq[1]] = kmer_freq[2]
        end

        if length(keys(kmer_dict)) <= 1
            error("Insufficient k-mers found!")
        else
            @info "Found $(length(keys(kmer_dict))) kmers for $variant"
        end

        freq = max_entropy(kmer_dict)
        filter!(e -> e[2] >= freq, kmer_dict)
        variant_kmers[variant] = Set(keys(kmer_dict))

        for kmer_values in values(var_hash)
            union!(all_exclusive_kmers, Set([kmer_values[1]]))
        end

    end

    # Measure all te intersections
    kmers_in_multiple_variants = Set{String}()
    variant_list = collect(keys(variant_kmers))

    @inbounds for i in eachindex(variant_list)
        for j in (i+1):length(variant_list)
            v1 = variant_list[i]
            v2 = variant_list[j]
            shared = intersect(variant_kmers[v1], variant_kmers[v2])
            union!(kmers_in_multiple_variants, shared)
        end
    end

    union_exclusive = setdiff(all_exclusive_kmers, kmers_in_multiple_variants)
    return all_exclusive_kmers
end


# Function to extract discriminatives features from each class
function extractFeaturesTemplate(
    wnwPercent::Float32,
    groupName::String,
    variantDirPath::String,
    k_len::Int,
    histogramThreshold::Float16=Float16(0.8))

    @info "Threads:" Threads.nthreads()
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"

    try
        mkpath(cachdir)
    catch e
        @error "create cache direcotry failed" exception = (e, catch_backtrace())
    end

    variantDirs::Vector{String} = readdir(variantDirPath)

    outputs = Vector{Tuple{String,Tuple{Vector{UInt16},BitArray}}}(undef, length(variantDirs))

    kmerset::Set{String} = get_exclusive_kmers(k_len, variantDirPath)
    DataIO.save_cache("$cachdir/kmerset.dat", kmerset)

    @info "Found: " kmerset
    # kmerset = Set{String}()

    # for variant in variantDirs
    #     variantKmers = DataIO.read_pickle_data("$variantDirPath/$variant/$(variant)_ExclusiveKmers.sav")
    #     union!(kmerset, Set(variantKmers))
    # end


    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        println("Processing $variant")
        cache_path = "$cachdir/$(variant)_outmask.dat"

        cache::Union{Nothing,Tuple{String,Tuple{Vector{UInt16},BitArray}}} = DataIO.load_cache(cache_path)

        if !isnothing(cache)
            @info "Using cached data from $cache_path"
            outputs[v] = cache
        else

            sequences::Vector{String} = DataIO.loadStringSequences("$variantDirPath/$variant")
            minSeqLength::UInt64 = minimum(map(length, sequences))
            wnwSize::UInt64 = ceil(UInt64, minSeqLength * wnwPercent)
            kmer_vector::Vector{String} = collect(kmerset)
            @info "Window size = $(Int(wnwSize))"

            data::Tuple{String,Tuple{Vector{UInt16},BitArray}} = (variant, _wndwExlcusiveKmersHistogram_bytes(kmer_vector, wnwSize, sequences, histogramThreshold))
            outputs[v] = data
            DataIO.save_cache(cache_path, data)
        end
    end

end


function compute_hash(s::String)::UInt64
    h = UInt64(0)
    # base = UInt64(4^length(s))
    base = UInt64(5)

    for char in s
        h = h * base + UInt64(char)
    end

    return h
end

function getOccursin_rolling_hash(
    sequence::String,
    kmer_hash_map::Dict{UInt64,Vector{String}},
    k_len::Int
)::Vector{Int}
    seq_len = length(sequence)
    positions = Int[]

    if seq_len < k_len
        return positions
    end

    base = UInt64(5)

    # Calculate base^(k_len-1) for rolling hash
    power = UInt64(1)
    for i in 1:(k_len-1)
        power *= base
    end

    # Calculate initial hash
    current_hash = UInt64(0)
    @inbounds for i in 1:k_len
        current_hash = current_hash * base + UInt64(sequence[i])
    end

    # Check first k-mer
    if haskey(kmer_hash_map, current_hash)
        kmer_candidate = SubString(sequence, 1, k_len)
        if String(kmer_candidate) in kmer_hash_map[current_hash]
            push!(positions, 1)
        end
    end

    # Roll through sequence
    @inbounds for i in (k_len+1):seq_len
        # Rolling hash: remove leftmost char, add rightmost char
        current_hash = current_hash - UInt64(sequence[i-k_len]) * power
        current_hash = current_hash * base + UInt64(sequence[i])

        # Check if hash matches any k-mer
        if haskey(kmer_hash_map, current_hash)
            kmer_candidate = SubString(sequence, i - k_len + 1, i)
            if String(kmer_candidate) in kmer_hash_map[current_hash]
                push!(positions, i - k_len + 1)
            end
        end
    end

    return positions
end



function _wndwExlcusiveKmersHistogram_bytes(
    exclusiveKmers::Vector{String},
    wndwSize::UInt64,
    sequences::Vector{String},
    histogramThreshold::Float16
)::Tuple{Vector{UInt16},BitArray}

    @assert all(≤(wndwSize), length.(exclusiveKmers)) "All k-mers must be ≤ window size"

    k_len::Int = length(exclusiveKmers[1])
    maxSeqLen = maximum(length, sequences)
    total_windows = maxSeqLen - wndwSize + 1

    kmer_hash_map = Dict{UInt64,Vector{String}}()

    for kmer in exclusiveKmers
        h = compute_hash(kmer)
        if !haskey(kmer_hash_map, h)
            kmer_hash_map[h] = String[]
        end
        push!(kmer_hash_map[h], kmer)
    end

    for (key, value) in kmer_hash_map
        if length(value) > 1
            @show key value
        end
    end


    @floop for seq in sequences
        positions = getOccursin_rolling_hash(seq, kmer_hash_map, k_len)

        window_coverage = falses(total_windows)

        @inbounds for pos in positions
            start_window = max(1, pos - Int(wndwSize) + k_len)
            end_window = min(total_windows, pos)

            if start_window <= end_window
                window_coverage[start_window:end_window] .= true
            end
        end

        @reduce(histogram = zeros(UInt32, total_windows) .+ UInt32.(window_coverage))
    end

    histogram_u16 = UInt16.(min.(histogram, typemax(UInt16)))
    threshold_count = UInt32(ceil(length(sequences) * histogramThreshold))
    threshold_mincount = UInt32(ceil(length(sequences) * 0.15))

    marked = falses(maxSeqLen)
    h_len = length(histogram)
    # hist = lowpass_filter(histogram_u16)
    @inbounds for i in eachindex(histogram)
        if (histogram[i] >= threshold_count)
            # end_pos = min(maxSeqLen, i + Int(wndwSize) - 1)
            # marked[i:end_pos] .= true
            init = max(1, i - 2)
            endi = min(h_len, i + 2)
            while (init <= endi)
                end_pos = min(maxSeqLen, init + Int(wndwSize) - 1)
                marked[init:end_pos] .= true
                init += 1
            end
        end
    end
    return histogram_u16, marked
end

function lowpass_filter(hist)
    freq = rfftfreq(length(hist))

    var_fft = rfft(hist)
    # Generate frequency bins for real FFT
    # For rfft, we only need frequencies up to Nyquist frequency
    freq = rfftfreq(length(hist))
    @inbounds for i in eachindex(freq)
        if freq[i] >= 0.01
            var_fft[i] = 0
        end
    end
    return irfft(var_fft, length(hist))
end

function getOccursin(
    sequence::String,
    kmerset::Vector{String})::Vector{Int}

    occurrences_pos = Set{Int}()
    @inbounds for i in eachindex(kmerset)
        regex = Regex("\\Q$(kmerset[i])\\E")
        match_positions = [m.offset for m in eachmatch(regex, sequence, overlap=true)]
        union!(occurrences_pos, Set(match_positions))
    end
    return collect(occurrences_pos)
end

function occursinKmer(
    kmer::String,
    windowBuffer::Union{SubArray,Vector{UInt8}}
)::Bool
    return occursin(Regex("\\Q$kmer\\E"), String(windowBuffer))
end

#=
Codeunits "regex", beacuse is a byte comparison 
=#
function occursinKmerBit(
    kmer::Base.CodeUnits,
    windowBuffer::Union{SubArray,Vector{UInt8}}
)::Bool
    wlen = length(windowBuffer)
    klen = length(kmer)

    @inbounds for i in 1:(wlen-klen+1)
        match = true
        for j in 1:klen
            (windowBuffer[i+j-1] ≠ kmer[j]) && (match = false; break)
        end
        match && return true
    end
    return false

end

end


