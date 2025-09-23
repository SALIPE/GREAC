module Report

include("DataIO.jl")

using Plots, .DataIO
using Plots: pdf

export Report

function generate_report_pdf(
    wnwPercent::Float32,
    groupName::String,
    model,
    classification_probs,
    results,
    output_file="analysis_report.pdf")
    # Set the backend to PDF-compatible backend
    gr()

    # Create the first plot: Variants Regions Behavior
    plt1 = plot(title="Variants Regions Behavior - $(length(model.regions)) windows",
        titlefont=10,
        guidefontsize=8,
        legend=:outertop,
        legendcolumns=5,
        dpi=300,
        size=(800, 600))

    for (class, distribution) in model.class_string_probs
        plot!(plt1, distribution, label=class)
    end

    ylabel!(plt1, "Frequency (Fwr)")
    xlabel!(plt1, "Extracted Regions (Features)")

    # Create the second plot: K-mer presence - Windows Histogram
    plt2 = plot(title="K-mer presence - Windows Histogram",
        titlefont=10,
        guidefontsize=8,
        legend=:outertop,
        legendcolumns=5,
        dpi=300,
        size=(800, 600))

    for class in model.classes
        try
            variant, (hist, mask) = DataIO.load_cache(
                "$(homedir())/.project_cache/$groupName/$wnwPercent/$(class)_outmask.dat")
            plot!(plt2, hist, label=class)
        catch e
            @warn "Failed to load data for class $class: $e"
        end
    end

    ylabel!(plt2, "Sequence presence amount")
    xlabel!(plt2, "Windows")

    summary_plot = create_summary_text_plot(
        groupName, wnwPercent, results, model)


    # Create membership summary plot
    # membership_plot = create_membership_detailed_plot(classification_probs)

    # Create confusion matrix heatmap
    cm_plot = create_confusion_matrix_plot(results[:confusion_matrix], model.classes)

    # Create a combined layout with plots and text summaries
    combined_plot = plot(plt1, plt2, cm_plot, summary_plot,
        layout=(5, 1),
        size=(800, 2400),
        margin=8Plots.mm)

    try
        savefig(combined_plot, output_file)
        println("PDF report successfully generated: $output_file")
        return output_file
    catch e
        error("Failed to save PDF: $e")
    end
end


function create_summary_text_plot(
    groupName,
    wnwPercent,
    results,
    model)
    summary_text = """
    BENCHMARK RESULTS SUMMARY
    ========================
    Group: $groupName
    Window Percentage: $wnwPercent
    Windows: $(length(model.regions))
    K-mer Set Size: $(length(model.kmerset))
    Final Length: $(count_region_length(model.regions))

    MACRO METRICS:
    - F1 Score: $(round(results[:macro][:f1], digits=4))
    - Precision: $(round(results[:macro][:precision], digits=4))
    - Recall: $(round(results[:macro][:recall], digits=4))

    PER-CLASS F1 SCORES:
    """

    for (class_name, metrics) in results[:per_class]
        summary_text *= "$class_name: $(round(metrics[:f1], digits=4))\n"
    end

    text_plot = plot(xlims=(0, 1), ylims=(0, 1),
        title="Benchmark Results Summary",
        titlefont=12,
        showaxis=false,
        grid=false,
        legend=false,
        size=(800, 400))

    annotate!(text_plot, 0.05, 0.95, text(summary_text, :left, :top, 8, :black))

    return text_plot
end
function create_membership_detailed_plot(classification_probs)
    # Build comprehensive membership text
    membership_text = "CLASSIFICATION MEMBERSHIPS - DETAILED RESULTS\n"
    membership_text *= "="^50 * "\n\n"

    for (key, value) in classification_probs
        membership_text *= "########### $(uppercase(key)) ############\n"

        for i in eachindex(value)
            try
                id, classification = value[i]
                membership_text *= "--- Classification $id ---\n"

                # Sort probabilities for better readability
                sorted_probs = sort(collect(classification), by=x -> x[2], rev=true)

                for (class_name, probability) in sorted_probs
                    prob_str = string(round(probability, digits=4))
                    membership_text *= "$class_name: $prob_str\n"
                end
                membership_text *= "\n"

            catch e
                @warn "Error found: $e"
                continue
            end
        end
        membership_text *= "\n"
    end

    # Create text plot with scrollable content
    # Split text into manageable chunks for display
    lines = split(membership_text, '\n')
    max_lines_per_plot = 40  # Adjust based on readability

    if length(lines) <= max_lines_per_plot
        # Single plot if content fits
        display_text = join(lines, '\n')
        text_plot = plot(xlims=(0, 1), ylims=(0, 1),
            title="Classification Memberships - Detailed",
            titlefont=12,
            showaxis=false,
            grid=false,
            legend=false,
            size=(800, 600))

        annotate!(text_plot, 0.02, 0.98, text(display_text, :left, :top, 7, :black))
    else
        # Multiple sections if content is too long
        display_text = join(lines[1:min(max_lines_per_plot, length(lines))], '\n')
        if length(lines) > max_lines_per_plot
            display_text *= "\n\n... (Content truncated for display) ..."
            display_text *= "\nTotal categories: $(length(classification_probs))"
            display_text *= "\nTotal samples: $(sum(length(v) for v in values(classification_probs)))"
        end

        text_plot = plot(xlims=(0, 1), ylims=(0, 1),
            title="Classification Memberships - Detailed",
            titlefont=12,
            showaxis=false,
            grid=false,
            legend=false,
            size=(800, 600))

        annotate!(text_plot, 0.02, 0.98, text(display_text, :left, :top, 7, :black))
    end

    return text_plot
end


function create_confusion_matrix_plot(confusion_matrix, class_names)
    # Ensure we're working with the matrix data directly
    cm_data = Matrix(confusion_matrix)  # Convert to Matrix if not already

    # Create heatmap
    cm_plot = heatmap(cm_data,
        title="Confusion Matrix",
        titlefont=12,
        xlabel="Predicted Class",
        ylabel="True Class",
        color=:Blues,
        size=(600, 500),
        dpi=300)

    # Add class labels if available
    if length(class_names) == size(cm_data, 1) == size(cm_data, 2)
        # Create tick positions and labels
        n_classes = length(class_names)
        tick_positions = 1:n_classes

        plot!(cm_plot,
            xticks=(tick_positions, class_names),
            yticks=(tick_positions, class_names),
            xrotation=45)
    else
        # Use generic class labels if names don't match matrix dimensions
        n_classes = size(cm_data, 1)
        generic_names = ["Class_$i" for i in 1:n_classes]
        tick_positions = 1:n_classes

        plot!(cm_plot,
            xticks=(tick_positions, generic_names),
            yticks=(tick_positions, generic_names),
            xrotation=45)
    end

    # Add text annotations with values
    max_val = maximum(cm_data)
    for i in 1:size(cm_data, 1)
        for j in 1:size(cm_data, 2)
            # Choose text color based on background intensity
            text_color = cm_data[i, j] > max_val * 0.5 ? :white : :black
            annotate!(cm_plot, j, i, text(string(cm_data[i, j]), 10, text_color, :center))
        end
    end

    return cm_plot
end


function count_region_length(regions)::Int
    total_length = 0
    for (i, e) in regions
        total_length += e - i
    end
    return total_length
end

end