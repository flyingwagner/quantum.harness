#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: export_affine_matching_minor.jl MODEL.mof.json[.gz] OUTPUT.tsv",
    )
    model_path, output_path = args
    ispath(output_path) && error("refusing to overwrite output: $output_path")
    problem = extract_exact_problem(model_path)
    analysis = affine_peeling_analysis(problem)
    matching = coupled_structural_matching(problem, analysis)
    matching.full_row_structural_rank ||
        error("coupled affine core lacks a full structural row matching")

    column_position = Dict(
        variable_index => position
        for (position, (_, variable_index)) in
            enumerate(matching.row_to_column)
    )
    nonzeros = 0
    open(output_path, "w") do io
        println(
            io,
            "row_position\tcolumn_position\trow_index\tvariable_index\tnumerator\tdenominator",
        )
        for (row_position, (row_index, _)) in
            enumerate(matching.row_to_column)
            for (variable_index, coefficient) in
                problem.equalities[row_index].terms
                haskey(column_position, variable_index) || continue
                println(
                    io,
                    row_position,
                    '\t',
                    column_position[variable_index],
                    '\t',
                    row_index,
                    '\t',
                    variable_index,
                    '\t',
                    numerator(coefficient),
                    '\t',
                    denominator(coefficient),
                )
                nonzeros += 1
            end
            if row_position % 500 == 0
                println(
                    "progress\tmatching_minor_rows\t",
                    row_position,
                    '/',
                    length(matching.row_to_column),
                )
                flush(stdout)
            end
        end
    end
    println("dimension\t", length(matching.row_to_column))
    println("nonzero_coefficients\t", nonzeros)
    println("optimizer_invoked\tfalse")
    println("output\t", abspath(output_path))
    return 0
end

exit(main())
