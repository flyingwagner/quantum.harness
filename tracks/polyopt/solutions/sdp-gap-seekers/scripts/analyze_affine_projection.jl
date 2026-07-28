#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) == 1 || error(
        "usage: analyze_affine_projection.jl MODEL.mof.json[.gz]",
    )
    problem = extract_exact_problem(args[1])
    analysis = affine_peeling_analysis(problem)
    println("variable_count\t", problem.variable_count)
    println("original_equality_rows\t", length(problem.equalities))
    println("unique_homogeneous_rows\t", length(analysis.unique_row_indices))
    println("duplicate_rows_removed\t", analysis.duplicate_rows_removed)
    println("peeled_unique_rows\t", analysis.peeled_unique_row_count)
    println("peeled_original_rows\t", analysis.peeled_original_row_count)
    println(
        "coupled_unique_rows\t",
        length(analysis.coupled_unique_row_indices),
    )
    println("coupled_original_rows\t", length(analysis.coupled_row_indices))
    println("coupled_columns\t", length(analysis.coupled_column_indices))
    println(
        "columns_outside_coupled_core\t",
        analysis.columns_outside_coupled_core,
    )
    println("optimizer_invoked\tfalse")
    return isempty(analysis.coupled_unique_row_indices) ? 0 : 2
end

exit(main())
