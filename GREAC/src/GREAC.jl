module GREAC

include("modules/DataIO.jl")
include("modules/Report.jl")
include("modules/RegionExtraction.jl")
include("modules/ClassificationModel.jl")

using FLoops,
    FASTX,
    LinearAlgebra,
    Normalization,
    Statistics,
    ArgParse,
    .DataIO,
    .RegionExtraction,
    .ClassificationModel,
    .Report,
    CSV,
    DataFrames,
    Dates

export GREAC


function greacClassificationFile(
    file_path::String,
    outputdir::Union{Nothing,String},
    wnwPercent::Float32,
    groupName::String,
    metric::Union{Nothing,String},
    use_xg::Bool
)

    model_name::String = "$(homedir())/.project_cache/$groupName/$wnwPercent/$groupName-multiclass"
    modelCachedFile = "$(homedir())/.project_cache/$groupName/$wnwPercent/kmers_distribution.dat"
    model::Union{Nothing,ClassificationModel.MultiClassModel} = DataIO.load_cache(modelCachedFile)

    kmerset::Vector{String} = collect(model.kmerset)
    regions::Vector{Tuple{Int,Int}} = model.regions

    classification_probs = Dict{String,Vector{Tuple{String,String,Dict{String,Float64}}}}()
    # predict_raw predict_membership (model, metric)
    classify = Base.Fix1(ClassificationModel.predict_membership, (model, metric, use_xg, model_name))

    total = DataIO.countSequences(file_path)

    chunk_size = 10000
    chunk_init::Int = 1

    @info "Classyfing $total sequences:"
    classifications = Vector{Tuple{String,String,Dict{String,Float64}}}(undef, total)

    while chunk_init <= total

        chunk_end = min(chunk_init + chunk_size - 1, total)
        current_chunk_size = chunk_end - chunk_init + 1

        classeqs::Vector{Tuple{String,Base.CodeUnits}} = DataIO.loadCodeUnitsSequences(file_path, chunk_init, chunk_end)

        inner_classifications = Vector{Tuple{String,String,Dict{String,Float64}}}(undef, current_chunk_size)

        @floop for local_idx in 1:current_chunk_size
            id, seq::Base.CodeUnits = classeqs[local_idx]

            kmer_distribution, kmer_counts = ClassificationModel.sequence_kmer_distribution_optimized(
                regions, seq, kmerset
            )
            seq_distribution = kmer_distribution ./ length(kmerset)

            if !iszero(sum(seq_distribution))
                cl, memberships = classify((seq_distribution, seq, kmer_counts))
                inner_classifications[local_idx] = (id, cl, memberships)
            end
        end

        classifications[chunk_init:chunk_end] = inner_classifications

        @info "Chunk processed $chunk_init - $chunk_end"
        chunk_init = chunk_end + 1
    end

    classification_probs[class] = classifications

    if !isnothing(outputdir)
        MEMBERSHIPS = "$outputdir/classifications_$groupName.csv"
        mkpath(outputdir)

        open(MEMBERSHIPS, "a") do io

            if filesize(MEMBERSHIPS) == 0
                types = join(model.classes, ",")
                write(io, "id," * types * ",predicted_label")
            end

            for (key, value) in classification_probs
                for i in eachindex(value)
                    try
                        id, classification = value[i]
                        line = "\n$id,"
                        for (_, cl, probability) in classification
                            line = line * "$(round(probability, digits=4)),"
                        end
                        line = line * "$cl"
                        write(io, line)
                    catch e
                        # @warn "Erro encontrado: $e"
                        continue
                    end
                end
            end
        end
    end
end


function greacClassification(
    folderPath::String,
    outputdir::Union{Nothing,String},
    wnwPercent::Float32,
    groupName::String,
    metric::Union{Nothing,String},
    use_xg::Bool
)

    model_name::String = "$(homedir())/.project_cache/$groupName/$wnwPercent/$groupName-multiclass"
    modelCachedFile = "$(homedir())/.project_cache/$groupName/$wnwPercent/kmers_distribution.dat"
    model::Union{Nothing,ClassificationModel.MultiClassModel} = DataIO.load_cache(modelCachedFile)

    y_true = String[]
    y_pred = String[]
    kmerset::Vector{String} = collect(model.kmerset)
    regions::Vector{Tuple{Int,Int}} = model.regions

    classification_probs = Dict{String,Vector{Tuple{String,String,Dict{String,Float64}}}}()
    # predict_raw predict_membership (model, metric)
    classify = Base.Fix1(ClassificationModel.predict_membership, (model, metric, use_xg, model_name))

    for class in model.classes
        file_path::String = "$folderPath/$class"
        total = 0

        try
            total = DataIO.countSequences(file_path)
        catch
            file_path = "$folderPath/$class.fasta"
            total = DataIO.countSequences(file_path)
        end


        chunk_size = 10000
        chunk_init::Int = 1

        @info "Classyfing $class $total sequences:"
        classifications = Vector{Tuple{String,String,Dict{String,Float64}}}(undef, total)
        local_y_pred = String[]
        local_y_true = String[]

        while chunk_init <= total

            chunk_end = min(chunk_init + chunk_size - 1, total)
            current_chunk_size = chunk_end - chunk_init + 1

            classeqs::Vector{Tuple{String,Base.CodeUnits}} = DataIO.loadCodeUnitsSequences(file_path, chunk_init, chunk_end)

            inner_classifications = Vector{Tuple{String,String,Dict{String,Float64}}}(undef, current_chunk_size)
            inner_y_pred = Vector{String}(undef, current_chunk_size)

            @floop for local_idx in 1:current_chunk_size
                id, seq::Base.CodeUnits = classeqs[local_idx]

                kmer_distribution, kmer_counts = ClassificationModel.sequence_kmer_distribution_optimized(
                    regions, seq, kmerset
                )
                seq_distribution = kmer_distribution ./ length(kmerset)

                if !iszero(sum(seq_distribution))
                    cl, memberships = classify((seq_distribution, seq, kmer_counts))
                    inner_classifications[local_idx] = (id, cl, memberships)
                    inner_y_pred[local_idx] = cl
                else
                    inner_y_pred[local_idx] = ""
                end
            end

            classifications[chunk_init:chunk_end] = inner_classifications

            for pred in inner_y_pred
                if !(pred == "")
                    push!(local_y_pred, pred)
                    push!(local_y_true, class)
                end
            end

            @info "Chunk processed $chunk_init - $chunk_end"
            chunk_init = chunk_end + 1
        end

        append!(y_pred, local_y_pred)
        append!(y_true, local_y_true)
        classification_probs[class] = classifications

    end

    results = compute_variant_metrics(model.classes, y_true, y_pred)


    @info "f1 = " results[:macro][:f1]

    if !isnothing(outputdir)
        RESULTS_CSV = "$outputdir/benchmark_results_$groupName.csv"
        MEMBERSHIPS = "$outputdir/classifications_$groupName.csv"
        mkpath(outputdir)

        open(MEMBERSHIPS, "a") do io

            if filesize(MEMBERSHIPS) == 0
                types = join(model.classes, ",")
                write(io, "id," * types * ",predicted_label,true_label")
            end
            for (key, value) in classification_probs
                for i in eachindex(value)
                    try
                        id, cl, classification = value[i]
                        line = "\n$id,"
                        for (_, probability) in classification
                            line = line * "$(round(probability, digits=4)),"
                        end
                        line = line * "$cl,$key"
                        write(io, line)
                    catch e
                        # @warn "Erro encontrado: $e"
                        continue
                    end
                end
            end
        end

        open(RESULTS_CSV, "a") do io
            # Write header if file is empty/new
            if filesize(RESULTS_CSV) == 0
                types = join(model.classes, ",")
                write(io, "wndwPercent,metric,windows,kmerset,final_len," * types * ",macro_f1,macro_precision,macro_recall,cm\n")
            end

            # Format data components
            cm = replace(string(results[:confusion_matrix]), "\n" => " | ")
            perclass = join([v[:f1] for (k, v) in results[:per_class]], ",")
            # Create CSV line
            line = join([
                    escape_string(string(wnwPercent)),
                    escape_string(string(metric)),
                    length(model.regions),
                    length(model.kmerset),
                    escape_string(string(count_region_length(model.regions))),
                    perclass,
                    results[:macro][:f1],
                    results[:macro][:precision],
                    results[:macro][:recall],
                    cm
                ], ",")

            write(io, line * "\n")
        end

        # Report.generate_report_pdf(
        #     wnwPercent,
        #     groupName,
        #     model,
        #     outputdir,
        #     results)
    end
    return results[:macro][:f1]
end

function compute_variant_metrics(
    variants::Vector{String},
    y_true::Vector{String},
    y_pred::Vector{String}
)
    n_classes = length(variants)

    length(y_true) == length(y_pred) || error("Input vectors must have the same length")
    all(v in variants for v in y_true) || error("Invalid variant in y_true")
    all(v in variants for v in y_pred) || error("Invalid variant in y_pred")

    # Create confusion matrix
    class_idx = Dict(v => i for (i, v) in enumerate(variants))
    cm = zeros(Int, n_classes, n_classes)

    for (t, p) in zip(y_true, y_pred)
        i = class_idx[t]
        j = class_idx[p]
        cm[i, j] += 1
    end

    metrics = Dict{String,Dict}()

    total_tp = 0
    total_fp = 0
    total_fn = 0

    for (i, variant) in enumerate(variants)
        tp = cm[i, i]
        fp = sum(cm[:, i]) - tp
        fn = sum(cm[i, :]) - tp

        # For micro metrics
        total_tp += tp
        total_fp += fp
        total_fn += fn

        precision = (tp + fp) == 0 ? 0.0 : tp / (tp + fp)
        recall = (tp + fn) == 0 ? 0.0 : tp / (tp + fn)
        f1 = (precision + recall) ≈ 0.0 ? 0.0 : 2 * ((precision * recall) / (precision + recall))
        support = sum(cm[i, :])

        metrics[variant] = Dict(
            :precision => round(precision, digits=4),
            :recall => round(recall, digits=4),
            :f1 => round(f1, digits=4),
            :support => support
        )
    end

    macro_precision = mean([m[:precision] for m in values(metrics)])
    macro_recall = mean([m[:recall] for m in values(metrics)])
    macro_f1 = mean([m[:f1] for m in values(metrics)])

    # micro_precision = (total_tp + total_fp) == 0 ? 0.0 : total_tp / (total_tp + total_fp)
    # micro_recall = (total_tp + total_fn) == 0 ? 0.0 : total_tp / (total_tp + total_fn)
    # micro_f1 = (micro_precision + micro_recall) ≈ 0.0 ? 0.0 :
    #            2 * ((micro_precision * micro_recall) / (micro_precision + micro_recall))

    return Dict(
        :confusion_matrix => cm,
        :classes => variants,
        :per_class => metrics,
        :macro => Dict(
            :precision => round(macro_precision, digits=4),
            :recall => round(macro_recall, digits=4),
            :f1 => round(macro_f1, digits=4)
        ),
        # :micro => Dict(
        #     :precision => round(micro_precision, digits=4),
        #     :recall => round(micro_recall, digits=4),
        #     :f1 => round(micro_f1, digits=4)
        # )
    )
end


function getKmersDistributionPerClass(
    wnwPercent::Float32,
    groupName::String,
    variantDirPath::String,
    use_xg::Bool,
    k_len::Int
)
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"
    model_name::String = "$(homedir())/.project_cache/$groupName/$wnwPercent/$groupName-multiclass"

    try
        mkpath(cachdir)
        @info "Cache dir created"
    catch e
        @error "create cache directory failed" exception = (e, catch_backtrace())
    end

    model::Union{Nothing,ClassificationModel.MultiClassModel} = DataIO.load_cache("$cachdir/kmers_distribution.dat")

    if !isnothing(model)
        @info "Model already processed from cached data from $cachdir"
    else

        variantDirs::Vector{String} = readdir(variantDirPath)

        # kmerset::Set{String} = RegionExtraction.get_exclusive_kmers(k_len, variantDirPath)
        kmerset::Set{String} = DataIO.load_cache("$cachdir/kmerset.dat")

        # kmerset = Set{String}()

        # @simd for variant in variantDirs
        #     variantKmers = DataIO.read_pickle_data("$variantDirPath/$variant/$(variant)_ExclusiveKmers.sav")
        #     union!(kmerset, Set(variantKmers))
        # end

        meta_data = Dict{String,Int}()
        byte_seqs = Dict{String,Vector{Base.CodeUnits}}()

        for variant in variantDirs
            try
                byte_seqs[variant] = DataIO.loadCodeUnitsSequences("$variantDirPath/$variant/$variant.fasta")
            catch
                byte_seqs[variant] = DataIO.loadCodeUnitsSequences("$variantDirPath/$variant")
            end

            meta_data[variant] = length(byte_seqs[variant])
        end

        @info meta_data
        @info "Total k-mer set" length(kmerset)

        distribution::ClassificationModel.MultiClassModel = ClassificationModel.fitMulticlass(
            kmerset,
            meta_data,
            byte_seqs,
            RegionExtraction.regionsConjuction(variantDirPath, wnwPercent, groupName),
            # RegionExtraction.filterRegions(variantDirPath, wnwPercent, groupName, Int(win_size), Int(maxSeqLen)),
            model_name,
            use_xg)

        DataIO.save_cache("$cachdir/kmers_distribution.dat", distribution)
        return distribution
    end
end


function count_region_length(regions)::Int
    total_length = 0
    for (i, e) in regions
        total_length += e - i
    end
    return total_length
end

function havein(pos, regions)

    for (i, e) in regions
        if pos >= i && pos <= e
            return pos, i, e
        end
    end
    return pos, 0, 0
end

function fitParameters(
    args,
    groupName::String,
    window::Float32,
)
    kmer::Int = args["k-len"]
    current_f1 = 0.0
    current_w = window
    current_metric = "manhattan"
    current_threshold = 0.5

    results = DataFrame(
        timestamp=DateTime[],
        window=Float32[],
        threshold=Float16[],
        kmer=Int[],
        metric=String[],
        f1_score=Float64[]
    )

    while window <= 0.003
        @info ">> Window " window
        threshold::Float16 = 0.5

        while threshold <= 0.8
            rm("$(homedir())/.project_cache/$(groupName)/$window";
                recursive=true, force=true)

            try
                RegionExtraction.extractFeaturesTemplate(
                    window,
                    groupName,
                    args["train-dir"],
                    kmer,
                    threshold
                )

                getKmersDistributionPerClass(
                    window,
                    groupName,
                    args["train-dir"],
                    args["classifier"],
                    kmer,
                )

                f1 = greacClassification(
                    args["test-dir"],
                    nothing,
                    window,
                    groupName,
                    current_metric,
                    args["classifier"]
                )

                push!(results, (
                    now(),
                    window,
                    threshold,
                    kmer,
                    current_metric,
                    f1
                ))

                if f1 > current_f1
                    current_f1 = f1
                    current_w = window
                    current_threshold = threshold
                    @info "New Best:" current_f1 current_w current_threshold kmer
                end

            catch e
                @error "Error during iteraction" exception = (e, catch_backtrace())
            end

            threshold += Float16(0.05)
        end

        window += Float32(0.0005)
    end

    @info "Best Parameters:" current_f1 current_w current_threshold kmer


    timestamp_str = Dates.format(now(), "yyyymmdd")
    output_dir = "./output-$kmer/reports-$groupName"
    mkpath(output_dir)

    csv_filename = "$(output_dir)/parameter_optimization_$(timestamp_str).csv"
    CSV.write(csv_filename, results)

    best_result = DataFrame(
        parameter=["window", "threshold", "kmer", "metric", "f1_score"],
        value=[string(current_w), string(current_threshold),
            string(kmer), current_metric, string(current_f1)]
    )

    best_filename = "$(output_dir)/best_parameters_$(timestamp_str).csv"
    CSV.write(best_filename, best_result)

    @info "Reports:" csv_filename best_filename

    return (
        f1=current_f1,
        window=current_w,
        threshold=current_threshold,
        kmer=kmer,
        metric=current_metric,
        all_results=results
    )
end


function add_benchmark_args!(settings)
    s = settings["benchmark"]
    @add_arg_table! s begin
        "--train-dir"
        help = "Training dataset path"
        required = true
        "--test-dir"
        help = "Test dataset path"
        required = true
        "--reference"
        help = "Reference sequence path"
        required = false
        "-k", "--k-len"
        help = "K-mer K value"
        required = true
        arg_type = Int
        "-m", "--metric"
        help = "Metric used for classification"
        required = false
        range_tester = (x -> x in ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"])
        "--threshold"
        help = "Window theshold consideration"
        required = false
        arg_type = Float16
        "-o", "--output-directory"
        help = "Where the files go"
        required = false
        "--use-gramep"
        help = "Use K-mer set from GRAMEP"
        action = :store_true
        "--classifier"
        help = "Classify sequences using XGBoost"
        action = :store_true
    end
end

function add_classification_args!(settings)
    s = settings["file-classification"]
    @add_arg_table! s begin
        "--file"
        help = "Test dataset path"
        required = true
        "-m", "--metric"
        help = "Metric used for classification"
        required = false
        range_tester = (x -> x in ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"])
        "--threshold"
        help = "Window theshold consideration"
        required = false
        arg_type = Float16
        "-o", "--output-directory"
        help = "Where the files go"
        required = false
        "--classifier"
        help = "Classify sequences using Random Forest"
        action = :store_true
    end
end

function add_extract_features_args!(settings)
    s = settings["extract-features"]
    @add_arg_table! s begin
        "--train-dir"
        help = "Training dataset path"
        required = true
        "--threshold"
        help = "Window theshold consideration"
        arg_type = Float16
        required = false
        "-k", "--k-len"
        help = "K-mer K value"
        required = true
        arg_type = Int
        "--classifier"
        help = "Classify sequences using XGBoost"
        action = :store_true
    end
end

function add_fit_parameters_args!(settings)
    s = settings["fit-parameters"]
    @add_arg_table! s begin
        "--train-dir"
        help = "Training dataset path"
        required = true
        "--test-dir"
        help = "Test dataset path"
        required = true
        "-k", "--k-len"
        help = "K-mer K value"
        required = false
        arg_type = Int
        "--classifier"
        help = "Classify sequences using XGBoost"
        action = :store_true
    end
end

function add_fasta_regions_args!(settings)
    s = settings["fasta-regions"]
    @add_arg_table! s begin
        "-i", "--input"
        help = "Training dataset path"
        required = false
    end
end

function handle_file_classification(args,
    groupName::String,
    window::Float32)
    @info "Starting File Classification " args window groupName

    greacClassificationFile(
        args["file"],
        args["output-directory"],
        window,
        groupName,
        args["metric"],
        args["classifier"]
    )
end


function handle_benchmark(args,
    groupName::String,
    window::Float32)
    @info "Starting benchmark" args window groupName
    @info "Starting model extraction"
    RegionExtraction.extractFeaturesTemplate(
        window,
        groupName,
        args["train-dir"],
        args["k-len"],
        args["reference"],
        args["use-gramep"],
        args["threshold"]
    )
    distribution = getKmersDistributionPerClass(
        window,
        groupName,
        args["train-dir"],
        args["classifier"],
        args["k-len"],
    )

    @info "Starting classification evaluation"
    greacClassification(
        args["test-dir"],
        args["output-directory"],
        window,
        groupName,
        args["metric"],
        args["classifier"]
    )
end

function extract_features(args,
    groupName::String,
    window::Float32)
    @info "Starting model extraction" args window groupName
    RegionExtraction.extractFeaturesTemplate(
        window,
        groupName,
        args["train-dir"],
        args["k-len"],
        args["reference"],
        args["use-gramep"],
        args["threshold"]
    )
    distribution = getKmersDistributionPerClass(
        window,
        groupName,
        args["train-dir"],
        args["classifier"],
        args["k-len"],
    )
end

function handle_extract_file_reads(args,
    groupName::String,
    wnwPercent::Float32)

    variantDirPath::Union{Nothing,String} = args["input"]
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"
    model::Union{Nothing,ClassificationModel.MultiClassModel} = DataIO.load_cache("$cachdir/kmers_distribution.dat")

    if isnothing(model)
        error("Model not found!")
    end

    if !isnothing(variantDirPath)
        variantDirs::Vector{String} = readdir(variantDirPath)
        for variant in variantDirs
            DataIO.createFastaRegionFile(
                "$variantDirPath/$variant/$variant.fasta",
                "$variantDirPath/$variant/$variant-extracted.fasta",
                model.regions)
        end
    end

    REGIONS = "regions_$groupName.bed"

    open(REGIONS, "w") do io
        for (init_pos, end_pos) in model.regions
            write(io, "\n$groupName\t$init_pos\t$end_pos")
        end
    end

    @info "Region BED file saved at: $REGIONS"

end

function julia_main()::Cint

    settings = ArgParseSettings(
        description="Genome Regions Extractor and Classifier",
        version="0.1",
        add_version=true
    )

    # Global options (apply to all commands)
    @add_arg_table! settings begin
        "--no-cache"
        help = "Remove cached files"
        action = :store_true
        "--group-name"
        help = "Process Group Name"
        required = true
        arg_type = String
        "-w", "--window"
        help = "Sliding window percent size"
        arg_type = Float32
        required = true
    end

    # Create subcommand structure
    @add_arg_table! settings begin
        ("benchmark", action=:command,
            help="Benchmark extract features model and classify creating and print confusion matrix")
        ("extract-features", action=:command,
            help="Exract features from k-mer set")
        ("fit-parameters", action=:command,
            help="Fit better params")
        ("fasta-regions", action=:command,
            help="Create file fasta region reads from extract")
        ("file-classification", action=:command,
            help="Classify file of sequences outputting the txt file")
    end

    # Add arguments for each subcommand
    add_benchmark_args!(settings)
    add_extract_features_args!(settings)
    add_fit_parameters_args!(settings)
    add_fasta_regions_args!(settings)
    add_classification_args!(settings)
    parsed_args = parse_args(settings)


    try

        if parsed_args["no-cache"]
            rm("$(homedir())/.project_cache/$(parsed_args["group-name"])/$(parsed_args["window"])"; recursive=true, force=true)
        end

        if parsed_args["%COMMAND%"] == "fit-parameters"
            fitParameters(parsed_args["fit-parameters"], parsed_args["group-name"], parsed_args["window"])
        elseif parsed_args["%COMMAND%"] == "extract-features"
            extract_features(parsed_args["extract-features"], parsed_args["group-name"], parsed_args["window"])
        elseif parsed_args["%COMMAND%"] == "benchmark"
            handle_benchmark(parsed_args["benchmark"], parsed_args["group-name"], parsed_args["window"])
        elseif parsed_args["%COMMAND%"] == "fasta-regions"
            handle_extract_file_reads(parsed_args["fasta-regions"], parsed_args["group-name"], parsed_args["window"])
        elseif parsed_args["%COMMAND%"] == "file-classification"
            handle_file_classification(parsed_args["file-classification"], parsed_args["group-name"], parsed_args["window"])
        end
    catch e
        @error "Error processing command" exception = (e, catch_backtrace())
    end
    return 0
end
end

GREAC.julia_main()
