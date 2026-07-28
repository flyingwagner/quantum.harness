#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) in (5, 6) || error(
        "usage: project_affine_ray_exact.jl MODEL.mof.json[.gz] " *
        "RAY.variables.tsv RATIONAL_TOLERANCE ZERO_THRESHOLD EXACT_RAY.tsv " *
        "[COEFFICIENT_TOLERANCE]",
    )
    model_path, ray_path, tolerance_text, zero_text, output_path = args[1:5]
    rational_tolerance = parse(Float64, tolerance_text)
    zero_threshold = parse(Float64, zero_text)
    coefficient_tolerance =
        length(args) == 6 ? parse(Float64, args[6]) : 0.0

    problem = extract_exact_problem(
        model_path;
        coefficient_tolerance=coefficient_tolerance,
    )
    moi_indices, ray = read_ray_values(ray_path)
    moi_indices == collect(1:length(ray)) ||
        error("ray indices are not canonical")
    length(ray) == problem.variable_count ||
        error("model/ray variable count mismatch")
    candidate, scale, pivot = normalize_rational_ray(
        ray;
        rational_tolerance=rational_tolerance,
        zero_threshold=zero_threshold,
    )
    analysis = affine_peeling_analysis(problem)
    initial = exact_residuals(problem, candidate)
    coupled_nonzero = count(
        row_index -> !iszero(initial[row_index]),
        analysis.coupled_unique_row_indices,
    )
    coupled_nonzero == 0 ||
        error("coupled affine core is not exactly satisfied")
    corrected, unresolved, corrections = correct_with_affine_peeling(
        problem,
        candidate,
        analysis,
    )
    final = exact_residuals(problem, corrected)
    final_nonzero = count(!iszero, final)
    final_nonzero == 0 ||
        error("exact affine projection left nonzero residuals")
    improvement = objective_improvement(problem, corrected)
    improvement > 0 ||
        error("exact affine projection lost objective improvement")
    write_exact_ray(output_path, corrected)

    println("variable_count\t", problem.variable_count)
    println("equality_count\t", length(problem.equalities))
    println("normalization\tmax_abs_equals_one")
    println("source_scale\t", scale)
    println("normalization_pivot\t", pivot)
    println("rational_tolerance\t", rational_tolerance)
    println("zero_threshold\t", zero_threshold)
    println("coefficient_tolerance\t", coefficient_tolerance)
    println("initial_nonzero_exact_residuals\t", count(!iszero, initial))
    println("initial_coupled_nonzero_exact_residuals\t", coupled_nonzero)
    println("peeling_corrections\t", length(corrections))
    println(
        "nonzero_peeling_corrections\t",
        count(correction -> !iszero(last(correction)), corrections),
    )
    println("unresolved_original_rows\t", length(unresolved))
    println("corrected_nonzero_exact_residuals\t", final_nonzero)
    println("objective_improvement_exact\t", improvement)
    println("objective_improvement_float\t", Float64(improvement))
    println("optimizer_invoked\tfalse")
    println("exact_ray_path\t", abspath(output_path))
    return 0
end

exit(main())
