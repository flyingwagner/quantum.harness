#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) in (2, 3) || error(
        "usage: postprocess_gap_ray.jl MODEL.mof.json[.gz] RAY.variables.tsv [EXACT_RAY.tsv]",
    )
    # At this TFIM calibration point the only non-integer affine coefficients
    # are ±1/2, ±γ, and ±2γ with γ=201/800. A 1e-14 reconstruction tolerance
    # maps their Float64 spellings to those intended rationals. The emitted
    # coefficient table below makes that model conversion auditable.
    problem = extract_exact_problem(args[1]; coefficient_tolerance=1e-14)
    moi_indices, ray = read_ray_values(args[2])
    moi_indices == collect(1:length(ray)) || error("ray indices are not canonical")
    length(ray) == problem.variable_count || error("model/ray variable count mismatch")

    candidate, scale, pivot = normalize_rational_ray(
        ray;
        rational_tolerance=1e-12,
        zero_threshold=1e-20,
    )
    initial_residuals = exact_residuals(problem, candidate)
    corrected, unresolved, corrections = correct_with_private_pivots(problem, candidate)
    corrected_residuals = exact_residuals(problem, corrected)
    minima = float_psd_minima(problem, corrected)
    improvement = objective_improvement(problem, corrected)
    psd_proofs = [rigorous_psd_proof(block, corrected) for block in problem.psd_blocks]
    coefficients = sort(unique(
        coefficient for row in problem.equalities for (_, coefficient) in row.terms
    ))
    expected_coefficients = BigRational[
        -4, -2, -1, -201//400, -1//2, -201//800,
        201//800, 1//2, 201//400, 1, 2, 4,
    ]
    problem.variable_count == 23_949 || error("unexpected TFIM variable count")
    length(problem.equalities) == 2_705 || error("unexpected TFIM equality count")
    [block.dimension for block in problem.psd_blocks] == [211, 50, 11, 14] ||
        error("unexpected TFIM PSD block inventory")
    coefficients == expected_coefficients ||
        error("unexpected TFIM exact coefficient inventory")

    println("variable_count\t", problem.variable_count)
    println("equality_count\t", length(problem.equalities))
    println("psd_dimensions\t", join((block.dimension for block in problem.psd_blocks), ','))
    println("normalization\tmax_abs_equals_one")
    println("source_scale\t", scale)
    println("normalization_pivot\t", pivot)
    println("rational_tolerance\t1.0e-12")
    println("zero_threshold\t1.0e-20")
    println("coefficient_reconstruction_tolerance\t1.0e-14")
    println("exact_coefficients\t", join(coefficients, ','))
    println("initial_nonzero_exact_residuals\t", count(!iszero, initial_residuals))
    println("private_pivot_corrections\t", length(corrections))
    println("unresolved_equalities\t", length(unresolved))
    println("corrected_nonzero_exact_residuals\t", count(!iszero, corrected_residuals))
    println("objective_improvement_positive\t", improvement > 0)
    println("objective_improvement_exact\t", improvement)
    println("objective_improvement_float\t", Float64(improvement))
    println("float_psd_minima\t", join(minima, ','))
    for (block_index, proof) in enumerate(psd_proofs)
        println(
            "psd_proof_$block_index\t",
            proof.status,
            ";active_dimension=", proof.active_dimension,
            ";exact_zero_rows=", join(proof.exact_zero_rows, ','),
            ";minimum_pivot_lower=", proof.minimum_pivot_lower,
        )
    end
    strict = isempty(unresolved) && all(iszero, corrected_residuals) &&
             improvement > 0 && all(proof -> proof.proved, psd_proofs)
    println("strict_certificate\t", strict)
    if length(args) == 3
        write_exact_ray(args[3], corrected)
        println("exact_ray_path\t", abspath(args[3]))
    end
    return strict ? 0 : 2
end

exit(main())
