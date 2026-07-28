#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: export_affine_coupled_core.jl MODEL.mof.json[.gz] OUTPUT.tsv",
    )
    model_path, output_path = args
    ispath(output_path) && error("refusing to overwrite output: $output_path")
    problem = extract_exact_problem(model_path)
    analysis = affine_peeling_analysis(problem)
    isempty(analysis.coupled_unique_row_indices) &&
        error("affine system has no coupled core")
    column_position = Dict(
        variable_index => position
        for (position, variable_index) in
            enumerate(analysis.coupled_column_indices)
    )
    nonzeros = 0
    open(output_path, "w") do io
        println(
            io,
            "row_position\tcolumn_position\trow_index\tvariable_index\tnumerator\tdenominator",
        )
        for (row_position, row_index) in
            enumerate(analysis.coupled_unique_row_indices)
            for (variable_index, coefficient) in
                problem.equalities[row_index].terms
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
                    "progress\tcoupled_core_rows\t",
                    row_position,
                    '/',
                    length(analysis.coupled_unique_row_indices),
                )
                flush(stdout)
            end
        end
    end
    println("rows\t", length(analysis.coupled_unique_row_indices))
    println("columns\t", length(analysis.coupled_column_indices))
    println("nonzero_coefficients\t", nonzeros)
    println("optimizer_invoked\tfalse")
    println("output\t", abspath(output_path))
    return 0
end

exit(main())
