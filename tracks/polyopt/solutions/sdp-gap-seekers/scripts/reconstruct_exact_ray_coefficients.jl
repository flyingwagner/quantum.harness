#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) == 4 || error(
        "usage: reconstruct_exact_ray_coefficients.jl " *
        "MODEL.mof.json[.gz] INPUT_EXACT_RAY.tsv " *
        "COEFFICIENT_TOLERANCE OUTPUT_EXACT_RAY.tsv",
    )
    model_path, input_path, tolerance_text, output_path = args
    coefficient_tolerance = parse(Float64, tolerance_text)
    coefficient_tolerance >= 0 ||
        error("coefficient tolerance must be nonnegative")

    problem = extract_exact_problem(
        model_path;
        coefficient_tolerance=coefficient_tolerance,
    )
    input_ray = read_exact_ray_values(input_path)
    length(input_ray) == problem.variable_count ||
        error("model/ray variable count mismatch")
    initial_residuals = exact_residuals(problem, input_ray)
    corrected, unresolved, corrections = correct_with_private_pivots(
        problem,
        input_ray,
    )
    isempty(unresolved) ||
        error(
            "coefficient reconstruction left " *
            "$(length(unresolved)) unresolved affine rows",
        )
    final_residuals = exact_residuals(problem, corrected)
    all(iszero, final_residuals) ||
        error("coefficient reconstruction left nonzero affine residuals")
    improvement = objective_improvement(problem, corrected)
    improvement > 0 ||
        error("coefficient reconstruction lost objective improvement")
    write_exact_ray(output_path, corrected)

    println("variable_count\t", problem.variable_count)
    println("equality_count\t", length(problem.equalities))
    println("coefficient_tolerance\t", coefficient_tolerance)
    println(
        "initial_nonzero_exact_residuals\t",
        count(!iszero, initial_residuals),
    )
    println("private_pivot_corrections\t", length(corrections))
    println("unresolved_equalities\t", length(unresolved))
    println(
        "corrected_nonzero_exact_residuals\t",
        count(!iszero, final_residuals),
    )
    println("objective_improvement_exact\t", improvement)
    println("optimizer_invoked\tfalse")
    println("exact_ray_path\t", abspath(output_path))
    return 0
end

exit(main())
