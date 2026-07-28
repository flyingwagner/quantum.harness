#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) >= 5 || error(
        "usage: scan_affine_projection_psd.jl MODEL.mof.json[.gz] " *
        "RAY.variables.tsv ZERO_THRESHOLD BLOCK[,BLOCK...] " *
        "RATIONAL_TOLERANCE...",
    )
    model_path, ray_path, zero_text, block_text = args[1:4]
    zero_threshold = parse(Float64, zero_text)
    block_indices = parse.(Int, split(block_text, ','))
    length(unique(block_indices)) == length(block_indices) ||
        error("PSD block indices must be unique")
    rational_tolerances = parse.(Float64, args[5:end])

    problem = extract_exact_problem(model_path)
    all(index -> index in eachindex(problem.psd_blocks), block_indices) ||
        error("PSD block index is out of range")
    moi_indices, ray = read_ray_values(ray_path)
    moi_indices == collect(1:length(ray)) ||
        error("ray indices are not canonical")
    length(ray) == problem.variable_count ||
        error("model/ray variable count mismatch")
    analysis = affine_peeling_analysis(problem)

    print(
        "rational_tolerance\tzero_threshold\tfinal_nonzero\t",
        "objective_improvement",
    )
    for block_index in block_indices
        print(
            "\tblock_",
            block_index,
            "_status\tblock_",
            block_index,
            "_proved\tblock_",
            block_index,
            "_proved_indefinite\tblock_",
            block_index,
            "_failed_pivot\tblock_",
            block_index,
            "_failed_pivot_upper",
        )
    end
    println("\toptimizer_invoked")
    flush(stdout)

    for rational_tolerance in rational_tolerances
        candidate, _, _ = normalize_rational_ray(
            ray;
            rational_tolerance=rational_tolerance,
            zero_threshold=zero_threshold,
        )
        initial = exact_residuals(problem, candidate)
        any(
            row_index -> !iszero(initial[row_index]),
            analysis.coupled_unique_row_indices,
        ) && error("coupled affine core is not exactly satisfied")
        corrected, _, _ = correct_with_affine_peeling(
            problem,
            candidate,
            analysis,
        )
        final_nonzero = count(
            !iszero,
            exact_residuals(problem, corrected),
        )
        print(
            rational_tolerance,
            '\t',
            zero_threshold,
            '\t',
            final_nonzero,
            '\t',
            Float64(objective_improvement(problem, corrected)),
        )
        for block_index in block_indices
            proof = rigorous_psd_proof(
                problem.psd_blocks[block_index],
                corrected;
                precision=256,
            )
            print(
                '\t',
                proof.status,
                '\t',
                proof.proved,
                '\t',
                proof.proved_indefinite,
                '\t',
                proof.failed_pivot,
                '\t',
                proof.failed_pivot_upper,
            )
        end
        println("\tfalse")
        flush(stdout)
    end
    return 0
end

exit(main())
