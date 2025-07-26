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
    BenchmarkTools,
    .DataIO,
    .RegionExtraction,
    .ClassificationModel,
    .Report

export GREAC


function greacClassification(
    folderPath::String,
    outputdir::Union{Nothing,String},
    wnwPercent::Float32,
    groupName::String,
    metric::Union{Nothing,String}
)

    model_name::String = "$(homedir())/.project_cache/$groupName/$wnwPercent/$groupName-multiclass.xgb"
    modelCachedFile = "$(homedir())/.project_cache/$groupName/$wnwPercent/kmers_distribution.dat"
    model::Union{Nothing,ClassificationModel.MultiClassModel} = DataIO.load_cache(modelCachedFile)

    classification_probs = Dict{String,Vector{Tuple{String,Dict{String,Float64}}}}()
    # predict_raw predict_membership (model, metric)
    classify = Base.Fix1(ClassificationModel.predict_membership, (model, metric, model_name))

    y_true = String[]
    y_pred = String[]
    kmerset::Vector{String} = collect(model.kmerset)
    regions::Vector{Tuple{Int,Int}} = model.regions

    for class in model.classes
        file_path::String = "$folderPath/$class.fasta"
        total = DataIO.countSequences(file_path)

        chunk_size = 10000
        chunk_init::Int = 1

        @info "Classyfing $class $total sequences:"
        classifications = Vector{Tuple{String,Dict{String,Float64}}}(undef, total)
        local_y_pred = String[]
        local_y_true = String[]

        while chunk_init <= total

            chunk_end = min(chunk_init + chunk_size - 1, total)
            current_chunk_size = chunk_end - chunk_init + 1

            classeqs::Vector{Tuple{String,Base.CodeUnits}} = DataIO.loadCodeUnitsSequences(file_path, chunk_init, chunk_end)

            inner_classifications = Vector{Tuple{String,Dict{String,Float64}}}(undef, current_chunk_size)
            inner_y_pred = Vector{String}(undef, current_chunk_size)

            @floop for local_idx in 1:current_chunk_size
                id, seq::Base.CodeUnits = classeqs[local_idx]

                kmer_distribution = ClassificationModel.sequence_kmer_distribution_optimized(
                    regions, seq, kmerset
                )
                seq_distribution = kmer_distribution ./ length(kmerset)

                if !iszero(sum(seq_distribution))
                    cl, memberships = classify(seq_distribution)
                    inner_classifications[local_idx] = (id, memberships)
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

    Report.generate_report_pdf(
        wnwPercent,
        groupName,
        model,
        classification_probs,
        results)

    @info "f1 = " results[:macro][:f1]

    if !isnothing(outputdir)
        RESULTS_CSV = "$outputdir/benchmark_results_$groupName.csv"
        MEMBERSHIPS = "$outputdir/memberships_$groupName.txt"
        mkpath(outputdir)

        open(MEMBERSHIPS, "w") do io

            for (key, value) in classification_probs
                write(io, "\n\n########### $(uppercase(key)) ############")

                for i in eachindex(value)

                    try
                        id, classification = value[i]
                        write(io, "\n--- Classificação $id ---\n")

                        for (class_name, probability) in classification
                            write(io, "$class_name: $(round(probability, digits=4)) \n")
                        end
                    catch e
                        @warn "Erro encontrado: $e"
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
    variantDirPath::String
)
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"
    model_name::String = "$(homedir())/.project_cache/$groupName/$wnwPercent/$groupName-multiclass.xgb"

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

        kmerset = Set{String}()

        @simd for variant in variantDirs
            variantKmers = DataIO.read_pickle_data("$variantDirPath/$variant/$(variant)_ExclusiveKmers.sav")
            union!(kmerset, Set(variantKmers))
        end

        meta_data = Dict{String,Int}()
        byte_seqs = Dict{String,Vector{Base.CodeUnits}}()

        # win_size = zero(UInt64)
        # maxSeqLen = zero(UInt64)

        for variant in variantDirs
            byte_seqs[variant] = DataIO.loadCodeUnitsSequences("$variantDirPath/$variant/$variant.fasta")
            meta_data[variant] = length(byte_seqs[variant])
            # minSeqLength::UInt64 = minimum(map(length, byte_seqs[variant]))
            # maxSeqLength::UInt64 = maximum(map(length, byte_seqs[variant]))
            # wnwSize::UInt64 = ceil(UInt64, minSeqLength * wnwPercent)

            # if (win_size == zero(UInt64) || wnwSize < win_size)
            #     win_size = wnwSize
            # end

            # if (maxSeqLen < maxSeqLength)
            #     maxSeqLen = maxSeqLength
            # end
        end

        @info meta_data
        # @show Int(win_size)
        # @show Int(maxSeqLen)

        distribution::ClassificationModel.MultiClassModel = ClassificationModel.fitMulticlass(
            kmerset,
            meta_data,
            byte_seqs,
            RegionExtraction.regionsConjuction(variantDirPath, wnwPercent, groupName),
            # RegionExtraction.filterRegions(variantDirPath, wnwPercent, groupName, Int(win_size), Int(maxSeqLen)),
            model_name)

        DataIO.save_cache("$cachdir/kmers_distribution.dat", distribution)
        return distribution
    end
end

function export_sars_pos(
    groupName::String,
    wnwPercent::Float32,
    distribution::ClassificationModel.MultiClassModel)

    RESULTS_CSV = "regions_val_$groupName.csv"

    open(RESULTS_CSV, "a") do io
        if filesize(RESULTS_CSV) == 0
            write(io, "wndwPercent,found,finalLength\n")
        end
        line = join([
                escape_string(string(wnwPercent)),
                escape_string(string(sars_pos_val(distribution))),
                escape_string(string(count_region_length(distribution.regions)))
            ], ",")

        write(io, line * "\n")
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

function sars_pos_val(model::ClassificationModel.MultiClassModel)::Int
    pos = [670,
        913,
        1059,
        2790,
        3037,
        3267,
        4184,
        5386,
        5648,
        5730,
        8393,
        8986,
        9053,
        9344,
        9534,
        9866,
        9867,
        10029,
        10135,
        10447,
        10449,
        11201,
        11537,
        11674,
        11750,
        12160,
        12880,
        13195,
        14257,
        14408,
        14676,
        15240,
        15279,
        15451,
        15714,
        16176,
        16466,
        16935,
        17039,
        17236,
        17410,
        18163,
        18171,
        23040,
        23055,
        23063,
        23075,
        23271,
        23403,
        23525,
        23593,
        23599,
        23604,
        23673,
        23709,
        23731,
        23948,
        24138,
        24424,
        24469,
        24503,
        24506,
        24748,
        24914,
        25000,
        25088,
        25416,
        25469,
        25563,
        25584,
        26060,
        26149,
        26270,
        26709,
        26767,
        28272,
        28512,
        28724,
        28881,
        28882,
        28883,
        28913,
        28916]

    havepos = Base.Fix2(havein, model.regions)

    found = map(havepos, pos)
    count = 0

    for (pos, i, _) in found
        if i != 0
            count += 1
        end
    end
    return count
end


function fitParameters(
    args,
    groupName::String,
    window::Float32
)
    current_f1 = 0
    current_w = window
    current_metric = "manhattan"
    current_threhold = 0.5

    while window <= 0.0025

        threhold::Float16 = 0.5
        # while threhold <= 0.9
        rm("$(homedir())/.project_cache/$(groupName)/$window"; recursive=true, force=true)
        try
            RegionExtraction.extractFeaturesTemplate(
                window,
                groupName,
                args["train-dir"],
                threhold)

            getKmersDistributionPerClass(
                window,
                groupName,
                args["train-dir"]
            )

            f1 = greacClassification(
                args["test-dir"],
                nothing,
                window,
                groupName,
                current_metric,
            )
            if f1 > current_f1
                current_f1 = f1
                current_w = window
                @info "New Best:" current_f1, current_w
            end
        catch e
            @error e
        end
        window += Float32(0.0005)
    end
    @info current_f1, current_w
end


function add_benchmark_args!(settings)
    s = settings["benchmark"]
    @add_arg_table! s begin
        "--train-dir"
        help = "Training dataset path"
        required = true
        "-m", "--metric"
        help = "Metric used for classification"
        required = false
        range_tester = (x -> x in ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"])
        "--test-dir"
        help = "Test dataset path"
        required = true
        "--threshold"
        help = "Window theshold consideration"
        required = false
        arg_type = Float16
        "-o", "--output-directory"
        help = "Where the files go"
        required = false
    end
end

function add_extract_features_args!(settings)
    s = settings["extract-features"]
    @add_arg_table! s begin
        "--train-dir"
        help = "Training dataset path"
        required = true
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
    end
end

function add_fasta_regions_args!(settings)
    s = settings["fasta-regions"]
    @add_arg_table! s begin
        "-i", "--input"
        help = "Training dataset path"
        required = true
    end
end

function add_performance_args!(settings)
    s = settings["performance-evaluation"]
    @add_arg_table! s begin
        "--train-dir"
        help = "Training dataset path"
        required = true
    end
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
        args["threshold"]
    )
    distribution = getKmersDistributionPerClass(
        window,
        groupName,
        args["train-dir"]
    )

    # export_sars_pos(groupName, window, distribution)
    @info "Starting classification evaluation"
    greacClassification(
        args["test-dir"],
        args["output-directory"],
        window,
        groupName,
        args["metric"]
    )
end

function extract_features(args,
    groupName::String,
    window::Float32)
    @info "Starting model extraction" args window groupName
    RegionExtraction.extractFeaturesTemplate(
        window,
        groupName,
        args["train-dir"]
    )
    distribution = getKmersDistributionPerClass(
        window,
        groupName,
        args["train-dir"]
    )
end

function handle_extract_file_reads(args,
    groupName::String,
    wnwPercent::Float32)

    variantDirPath::String = args["input"]
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"
    model::Union{Nothing,ClassificationModel.MultiClassModel} = DataIO.load_cache("$cachdir/kmers_distribution.dat")

    if isnothing(model)
        error("Model not found!")
    end

    variantDirs::Vector{String} = readdir(variantDirPath)

    for variant in variantDirs
        DataIO.createFastaRegionFile(
            "$variantDirPath/$variant/$variant.fasta",
            "$variantDirPath/$variant/$variant-extracted.fasta",
            model.regions)
    end

end

function handle_performance_evaluation(args, groupName::String)
    @info "Starting performance evaluation"

    # Configure benchmark parameters
    params = BenchmarkTools.DEFAULT_PARAMETERS
    params.seconds = 60      # Longer time budget for stable results
    params.evals = 1         # Better for multithreaded functions
    params.gcsample = true   # Collect GC statistics

    train_dir = args["train-dir"]
    windows = Float32[0.002, 0.004, 0.006, 0.008]

    # Storage for results
    results = Dict{String,Any}()

    # Benchmark each window configuration
    for (i, window) in enumerate(windows)

        # Feature extraction benchmark
        bench_feat = @benchmarkable(
                         RegionExtraction.extractFeaturesTemplate(w, $groupName, $train_dir),
                         setup = (GC.gc(); w = $window),  # Ensure clean state and fixed window
                         teardown = (GC.gc()),
                         evals = 1,
                         samples = 20
                     ) |> tune! |> run

        # Model fitting benchmark
        bench_model = @benchmarkable(
                          getKmersDistributionPerClass(w, $groupName, $train_dir),
                          setup = (GC.gc(); w = $window),
                          teardown = (GC.gc()),
                          evals = 1,
                          samples = 20
                      ) |> tune! |> run

        results["window_$i"] = Dict(
            :window => window,
            :feature => bench_feat,
            :model => bench_model
        )
    end

    open("benchmark_summary_$(groupName).txt", "w") do io
        println(io, "Performance Evaluation Report")
        println(io, "Group: ", groupName)
        println(io, "Threads Available: ", Threads.nthreads(), "\n")

        for (k, v) in results
            println(io, "\n", "-"^50)
            println(io, "Window: ", v[:window])

            println(io, "\nFeature Extraction:")
            show(io, MIME"text/plain"(), v[:feature])

            println(io, "\n\nModel Fitting:")
            show(io, MIME"text/plain"(), v[:model])

            println(io, "\nMemory Summary:")
            println(io, "  Feature - Allocs: ", v[:feature].allocs)
            println(io, "  Model - Allocs: ", v[:model].allocs)
        end
    end

    @info "Benchmark results saved to benchmark_summary_$(groupName).txt"
    return results
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
        ("performance-evaluation", action=:command,
            help="Evaluate function performance using bechmark tools")
        ("fit-parameters", action=:command,
            help="Fit better params")
        ("fasta-regions", action=:command,
            help="Create file fasta region reads from extract")
    end

    # Add arguments for each subcommand
    add_performance_args!(settings)
    add_benchmark_args!(settings)
    add_extract_features_args!(settings)
    add_fit_parameters_args!(settings)
    add_fasta_regions_args!(settings)
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
        elseif parsed_args["%COMMAND%"] == "performance-evaluation"
            handle_performance_evaluation(parsed_args["performance-evaluation"], parsed_args["group-name"])
        end
    catch e
        @error "Error processing command" exception = (e, catch_backtrace())
    end
    return 0
end
end

GREAC.julia_main()
