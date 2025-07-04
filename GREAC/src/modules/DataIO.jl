module DataIO
using FASTX, Pickle, Serialization, BioSequences

export DataIO

EIIP_NUCLEOTIDE = Dict{Char,Float32}([
    ('A', 0.1260),
    ('G', 0.0806),
    ('T', 0.1335),
    ('C', 0.1340)])

EIIP_AMINOACID = Dict{Char,Float32}([
    ('*', 0.0000),
    ('A', 0.0373),
    ('B', 0.0000),#non-exist
    ('C', 0.0829),
    ('D', 0.1263),
    ('E', 0.0058),
    ('F', 0.0946),
    ('G', 0.0050),
    ('H', 0.0242),
    ('I', 0.0000),
    ('J', 0.0000),#non-exist
    ('K', 0.0371),
    ('L', 0.0000),
    ('M', 0.0823),
    ('N', 0.0036),
    ('O', 0.0000),#non-exist 
    ('P', 0.0198),
    ('Q', 0.0761),
    ('R', 0.0959),
    ('S', 0.0829),
    ('T', 0.0941),
    ('U', 0.0000),#non-exist 
    ('V', 0.0057),
    ('W', 0.0548),
    ('X', 0.0000),#non-exist 
    ('Y', 0.0516),
    ('Z', 0.0000),#non-exist 
])

function padRNA(rna::LongSequence{RNAAlphabet{4}})
    pad_length = (3 - (length(rna) % 3)) % 3
    return rna * LongSequence{RNAAlphabet{4}}(repeat('N', pad_length))
end


function sequence2AminNumSerie(
    sequence::AbstractString
)::Vector{Float32}

    dna = LongSequence{DNAAlphabet{4}}(sequence)
    rna = convert(LongSequence{RNAAlphabet{4}}, dna)
    amn = BioSequences.translate(padRNA(rna))
    eiip = Vector{Float32}(undef, length(amn))

    @inbounds for (i, c) in enumerate(amn)
        eiip[i] = EIIP_AMINOACID[Char(c)]
    end
    return eiip
end

function sequence2NumericalSerie(
    sequence::AbstractString
)::Vector{Float32}

    eiip = Vector{Float32}(undef, length(sequence))
    @inbounds for (i, c) in enumerate(sequence)
        eiip[i] = EIIP_NUCLEOTIDE[Char(c)]
    end
    return eiip
end


function readRegionsFromBed(
    bedFilePath::String
)::Vector{Tuple{String,Int,Int}}

    regioes = Vector{Tuple{String,Int,Int}}()
    open(bedFilePath) do bed
        for linha in readlines(bed)
            if strip(linha) != ""  # Ignora linhas vazias
                campos = split(linha)
                cromossomo::String = campos[1]
                inicio::String = campos[2] # Começo da região (0-indexed)
                fim::String = campos[3]     # Fim da região (1-indexed)
                push!(regioes, (cromossomo, parse(Int, inicio), parse(Int, fim)))
            end
        end

    end
    return regioes
end

function readVectorFromFile(file::String, T::Type)::Vector{T}

    arr = Vector{T}()
    open(file, "r") do fi
        for linha in eachline(fi)
            push!(arr, parse(T, linha))
        end
    end
    return arr
end


function read_pickle_data(file_name::AbstractString)
    load_pickle = Pickle.load(file_name)
    return load_pickle
end

function save_cache(cache_path::String, data)
    try
        open(cache_path, "w") do io
            serialize(io, data)
        end
        @info "Cache saved successfully: $cache_path"
    catch e
        @error "Failed to save cache" exception = (e, catch_backtrace())
    end
end

function load_cache(cache_path::String)
    try
        isfile(cache_path) || return nothing
        open(io -> deserialize(io), cache_path, "r")
    catch e
        @error "Cache loading failed" exception = (e, catch_backtrace())
        nothing
    end
end

function loadStringSequences(
    file::String
)::Vector{String}

    sequences = String[]
    for record in open(FASTAReader, file)
        push!(sequences, sequence(String, record))
    end
    return sequences
end

function loadCodeUnitsSequences(
    file::String
)::Vector{Base.CodeUnits}

    sequences = Vector{Base.CodeUnits}()
    for record in open(FASTAReader, file)
        push!(sequences, codeunits(sequence(String, record)))
    end
    return sequences
end
function loadCodeUnitsSequences(
    file::String,
    chunk_init::Int,
    chunk_end::Int
)::Vector{Tuple{String,Base.CodeUnits}}

    sequences = Vector{Tuple{String,Base.CodeUnits}}()
    for (i, record) in enumerate(open(FASTAReader, file))
        if i >= chunk_init && i <= chunk_end
            id = identifier(record)
            seq = sequence(String, record)
            push!(sequences, (String(id), codeunits(seq)))
        end
    end
    return sequences
end

function countSequences(
    file::String
)::Int

    count = 0
    for record in open(FASTAReader, file)
        count += 1
    end
    return count
end

function createFastaRegionFile(
    filePath::String,
    outFile::String,
    regions::Vector{Tuple{Int,Int}}
)
    open(outFile, "w") do io
        writer = FASTX.FASTA.Writer(io)
        open(FASTX.FASTA.Reader, filePath) do reader
            for record in reader
                sequence_str = FASTX.FASTA.sequence(record)
                header = FASTX.FASTA.identifier(record)
                new_header = "$(header)_reducted"
                concatenated_seq = ""

                for (start_pos, end_pos) in regions
                    if start_pos >= 1 && end_pos <= length(sequence_str) && start_pos <= end_pos
                        subseq = sequence_str[start_pos:end_pos]
                        concatenated_seq *= subseq
                    else
                        @warn "Invalid region ($start_pos, $end_pos) for sequence $(header) of length $(length(sequence_str))"
                    end
                end
                new_record = FASTX.FASTA.Record(new_header, concatenated_seq)
                write(writer, new_record)
            end
        end
        close(writer)
    end
end


end
