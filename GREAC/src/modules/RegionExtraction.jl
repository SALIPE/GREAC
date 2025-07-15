module RegionExtraction

include("DataIO.jl")

using .DataIO,
    Serialization,
    AbstractFFTs,
    DSP,
    FASTX,
    FLoops

export RegionExtraction


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


# Function to extract discriminatives features from each class
function extractFeaturesTemplate(
    wnwPercent::Float32,
    groupName::String,
    variantDirPath::String,
    histogramThreshold::Float16=Float16(0.8)
)

    @info "Threads:" Threads.nthreads()
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"

    try
        mkpath(cachdir)
    catch e
        @error "create cache direcotry failed" exception = (e, catch_backtrace())
    end

    variantDirs::Vector{String} = readdir(variantDirPath)


    outputs = Vector{Tuple{String,Tuple{Vector{UInt16},BitArray}}}(undef, length(variantDirs))
    # varKmer = Dict{String,Vector{String}}()
    kmerset = Set{String}()

    for variant in variantDirs
        variantKmers = DataIO.read_pickle_data("$variantDirPath/$variant/$(variant)_ExclusiveKmers.sav")
        # varKmer[variant] = variantKmers
        union!(kmerset, Set(variantKmers))
    end


    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        println("Processing $variant")
        cache_path = "$cachdir/$(variant)_outmask.dat"

        cache::Union{Nothing,Tuple{String,Tuple{Vector{UInt16},BitArray}}} = DataIO.load_cache(cache_path)

        if !isnothing(cache)
            @info "Using cached data from $cache_path"
            outputs[v] = cache
        else

            sequences::Vector{String} = DataIO.loadStringSequences("$variantDirPath/$variant/$variant.fasta")
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
    base = UInt64(257)

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

    base = UInt64(257)

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


