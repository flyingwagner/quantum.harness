#!/usr/bin/env julia

using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) in (2, 3, 4, 5) || error(
        "usage: audit_exact_ray_psd.jl MODEL.mof.json[.gz] EXACT_RAY.tsv " *
        "[KERNEL_BLOCK KERNEL_SUPPORTS] [COEFFICIENT_TOLERANCE]",
    )
    coefficient_tolerance =
        length(args) in (3, 5) ? parse(Float64, args[end]) : 0.0
    problem = extract_exact_problem(
        args[1];
        coefficient_tolerance=coefficient_tolerance,
    )
    ray = read_exact_ray_values(args[2])
    has_kernel = length(args) in (4, 5)
    kernel_block = has_kernel ? parse(Int, args[3]) : nothing
    kernel_supports = has_kernel ? [
        parse.(Int, split(group, ','))
        for group in split(args[4], ';')
    ] : Vector{Vector{Int}}()
    length(ray) == problem.variable_count ||
        error("model/ray variable count mismatch")
    residuals = exact_residuals(problem, ray)
    improvement = objective_improvement(problem, ray)
    println("variable_count\t", problem.variable_count)
    println("equality_count\t", length(problem.equalities))
    println("nonzero_exact_residuals\t", count(!iszero, residuals))
    println("objective_improvement_exact\t", improvement)
    println("objective_improvement_positive\t", improvement > 0)
    println("coefficient_tolerance\t", coefficient_tolerance)
    for (block_index, block) in enumerate(problem.psd_blocks)
        println(
            "progress\tblock\t",
            block_index,
            "\tdimension\t",
            block.dimension,
        )
        flush(stdout)
        if block_index == kernel_block
            kernel = zeros(
                BigRational,
                block.dimension,
                length(kernel_supports),
            )
            for (kernel_index, support) in enumerate(kernel_supports)
                kernel[support, kernel_index] .= 1
            end
            proof = rigorous_psd_proof_with_exact_kernel(
                GapRayPostprocess.exact_block_matrix(block, ray),
                kernel,
                last.(kernel_supports);
                precision=256,
            )
        else
            proof = rigorous_psd_proof(block, ray; precision=256)
        end
        println(
            "block\t",
            block_index,
            "\tproved\t",
            proof.proved,
            "\tproved_indefinite\t",
            proof.proved_indefinite,
            "\tstatus\t",
            proof.status,
            "\tactive_dimension\t",
            proof.active_dimension,
            "\tfailed_pivot\t",
            proof.failed_pivot,
            "\tfailed_pivot_lower\t",
            proof.failed_pivot_lower,
            "\tfailed_pivot_upper\t",
            proof.failed_pivot_upper,
            "\tminimum_pivot_lower\t",
            proof.minimum_pivot_lower,
        )
        if proof.proved_indefinite && isnothing(kernel_block)
            matrix = GapRayPostprocess.exact_block_matrix(block, ray)
            active = [
                row
                for row in axes(matrix, 1)
                if any(!iszero, @view matrix[row, :])
            ]
            order = something(proof.failed_pivot)
            principal = matrix[active[1:order], active[1:order]]
            determinant = det(principal)
            println(
                "negative_principal_minor\tblock\t",
                block_index,
                "\torder\t",
                order,
                "\tnumerator\t",
                numerator(determinant),
                "\tdenominator\t",
                denominator(determinant),
            )
        end
        flush(stdout)
    end
    println("optimizer_invoked\tfalse")
    return 0
end

exit(main())
