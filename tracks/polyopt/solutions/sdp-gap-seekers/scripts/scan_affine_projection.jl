#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) >= 4 || error(
        "usage: scan_affine_projection.jl MODEL.mof.json[.gz] " *
        "RAY.variables.tsv ZERO_THRESHOLD RATIONAL_TOLERANCE...",
    )
    model_path, ray_path, zero_text = args[1:3]
    zero_threshold = parse(Float64, zero_text)
    rational_tolerances = parse.(Float64, args[4:end])

    problem = extract_exact_problem(model_path)
    moi_indices, ray = read_ray_values(ray_path)
    moi_indices == collect(1:length(ray)) ||
        error("ray indices are not canonical")
    length(ray) == problem.variable_count ||
        error("model/ray variable count mismatch")
    analysis = affine_peeling_analysis(problem)

    print(
        "rational_tolerance\tzero_threshold\tinitial_nonzero\t",
        "initial_coupled_nonzero\tcorrection_nonzero\t",
        "correction_max_abs\tfinal_nonzero\tobjective_improvement\t",
        "psd_minimum\tpsd_minima\toptimizer_invoked\n",
    )
    flush(stdout)
    for rational_tolerance in rational_tolerances
        candidate, _, _ = normalize_rational_ray(
            ray;
            rational_tolerance=rational_tolerance,
            zero_threshold=zero_threshold,
        )
        initial = exact_residuals(problem, candidate)
        corrected, _, corrections = correct_with_affine_peeling(
            problem,
            candidate,
            analysis,
        )
        final = exact_residuals(problem, corrected)
        minima = float_psd_minima(problem, corrected)
        correction_maximum = maximum(
            correction -> abs(Float64(last(correction))),
            corrections;
            init=0.0,
        )
        println(
            rational_tolerance,
            '\t',
            zero_threshold,
            '\t',
            count(!iszero, initial),
            '\t',
            count(
                row_index -> !iszero(initial[row_index]),
                analysis.coupled_unique_row_indices,
            ),
            '\t',
            count(correction -> !iszero(last(correction)), corrections),
            '\t',
            correction_maximum,
            '\t',
            count(!iszero, final),
            '\t',
            Float64(objective_improvement(problem, corrected)),
            '\t',
            minimum(minima),
            '\t',
            join(minima, ','),
            "\tfalse",
        )
        flush(stdout)
    end
    return 0
end

exit(main())
