#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function row_key(row)
    return Tuple(row.terms)
end

function block_ordinals(problem)
    ownership = Dict{Int,Set{Int}}()
    for (block_index, block) in enumerate(problem.psd_blocks)
        for row in block.coordinates, (ordinal, _) in row.terms
            push!(get!(ownership, ordinal, Set{Int}()), block_index)
        end
    end
    return ownership
end

function print_top_rows(problem, values, stats; count=15)
    ownership = block_ordinals(problem)
    order = sortperm(stats; by=stat -> -abs(stat.residual_high_precision))
    println("top_equality_rows")
    println("row\tterms\tresidual_high_precision\tbackward_error\trow_scaled_residual\tsummation_difference\tdominant_term")
    for row_index in order[1:min(count, length(order))]
        row = problem.equalities[row_index]
        stat = stats[row_index]
        dominant = isempty(row.terms) ? nothing : first(sort(
            row.terms;
            by=term -> -abs(Float64(last(term)) * values[first(term)]),
        ))
        dominant_text = if isnothing(dominant)
            "none"
        else
            ordinal, coefficient = something(dominant)
            blocks = sort!(collect(get(ownership, ordinal, Set{Int}())))
            contribution = Float64(coefficient) * values[ordinal]
            "v$(ordinal):coeff=$(coefficient):value=$(values[ordinal]):contribution=$(contribution):psd_blocks=$(join(blocks, ','))"
        end
        println(
            row_index, '\t', stat.term_count, '\t',
            stat.residual_high_precision, '\t', stat.backward_error, '\t',
            stat.row_scaled_residual, '\t', stat.summation_difference, '\t',
            dominant_text,
        )
    end
end

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: diagnose_ray_conditioning.jl MODEL.mof.json[.gz] RAY.variables.tsv",
    )
    problem = extract_exact_problem(args[1])
    moi_indices, values = read_ray_values(args[2])
    moi_indices == collect(1:length(values)) || error("ray indices are not canonical")
    stats = equality_conditioning(problem, values)
    block_scales = psd_block_scales(problem, values)
    pivots = private_pivot_rows(problem)
    homogeneous_groups = Dict{Any,Vector{Int}}()
    affine_groups = Dict{Any,Vector{Int}}()
    for (row_index, (row, offset)) in enumerate(zip(
        problem.equalities,
        problem.equality_offsets,
    ))
        push!(get!(homogeneous_groups, row_key(row), Int[]), row_index)
        push!(get!(affine_groups, (row_key(row), offset), Int[]), row_index)
    end
    homogeneous_duplicates = filter(
        pair -> length(last(pair)) > 1,
        collect(homogeneous_groups),
    )
    affine_duplicates = filter(
        pair -> length(last(pair)) > 1,
        collect(affine_groups),
    )

    println("variable_count\t", problem.variable_count)
    println("equality_count\t", length(problem.equalities))
    println("psd_dimensions\t", join((block.dimension for block in problem.psd_blocks), ','))
    println("ray_max_abs\t", maximum(abs, values))
    println("max_equality_residual_high_precision\t", maximum(
        stat -> abs(stat.residual_high_precision), stats,
    ))
    println("max_equality_residual_normalized_by_ray\t", maximum(
        stat -> abs(stat.residual_high_precision), stats,
    ) / maximum(abs, values))
    println("max_equality_backward_error\t", maximum(stat -> stat.backward_error, stats))
    println("max_row_scaled_residual\t", maximum(stat -> stat.row_scaled_residual, stats))
    println("max_float_summation_difference\t", maximum(stat -> stat.summation_difference, stats))
    println("rows_with_private_pivots\t", count(!isnothing, pivots))
    println("zero_homogeneous_rows\t", count(row -> isempty(row.terms), problem.equalities))
    println("homogeneous_duplicate_row_groups\t", length(homogeneous_duplicates))
    println("exact_affine_duplicate_row_groups\t", length(affine_duplicates))
    println("exact_affine_redundant_rows\t", sum(
        length(last(group)) - 1 for group in affine_duplicates;
        init=0,
    ))
    for (block_index, stat) in enumerate(block_scales)
        println(
            "psd_block_", block_index, '\t',
            "dimension=", stat.dimension,
            ";nonzero_coordinates=", stat.nonzero_coordinates,
            ";maximum=", stat.maximum,
            ";minimum_nonzero=", stat.minimum_nonzero,
        )
    end
    print_top_rows(problem, values, stats)
    return 0
end

exit(main())
