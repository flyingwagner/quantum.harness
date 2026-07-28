#!/usr/bin/env julia

using LinearAlgebra
import Clarabel
import JuMP

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

const MOI = JuMP.MOI

function main(args=ARGS)
    length(args) == 4 || error(
        "usage: find_affine_psd_repair.jl MODEL.mof.json[.gz] " *
        "PSD_BLOCK RATIONAL_TOLERANCE EXACT_DIRECTION.tsv",
    )
    model_path, block_text, tolerance_text, output_path = args
    block_index = parse(Int, block_text)
    rational_tolerances = parse.(Float64, split(tolerance_text, ','))
    all(>=(0), rational_tolerances) ||
        error("rational tolerances must be nonnegative")

    problem = extract_exact_problem(model_path)
    block_index in eachindex(problem.psd_blocks) ||
        error("PSD block index is out of range")
    analysis = affine_peeling_analysis(problem)
    coupled = Set(analysis.coupled_column_indices)
    peel_pivots = Set(last.(analysis.peel_order))
    free_columns = Set(setdiff(
        1:problem.variable_count,
        union(coupled, peel_pivots),
    ))
    block = problem.psd_blocks[block_index]
    block_columns = sort!(unique(
        column
        for row in block.coordinates
        for (column, _) in row.terms
        if column in free_columns
    ))
    isempty(block_columns) &&
        error("target PSD block has no affine-kernel free columns")

    basis_directions = Vector{Vector{BigRational}}()
    basis_matrices = Vector{Matrix{BigRational}}()
    for (basis_index, column) in enumerate(block_columns)
        seed = zeros(BigRational, problem.variable_count)
        seed[column] = 1
        direction, _, _ = correct_with_affine_peeling(
            problem,
            seed,
            analysis,
        )
        # `column` is outside the coupled core and every reverse-peeling pivot
        # is exact by construction. Auditing all original rows for every basis
        # vector is redundant and prohibitively expensive; the rationalized
        # combined direction is audited once below.
        push!(basis_directions, direction)
        push!(
            basis_matrices,
            GapRayPostprocess.exact_block_matrix(block, direction),
        )
        if basis_index % 25 == 0
            println(
                "progress\tbasis\t",
                basis_index,
                '/',
                length(block_columns),
            )
            flush(stdout)
        end
    end

    model = JuMP.Model(Clarabel.Optimizer)
    JuMP.set_silent(model)
    JuMP.set_optimizer_attribute(model, "tol_gap_abs", 1e-10)
    JuMP.set_optimizer_attribute(model, "tol_gap_rel", 1e-10)
    JuMP.set_optimizer_attribute(model, "tol_feas", 1e-10)
    JuMP.set_optimizer_attribute(model, "max_iter", 300)
    JuMP.@variable(model, weights[eachindex(block_columns)])
    JuMP.@variable(model, margin)
    matrix = [
        JuMP.AffExpr(0.0)
        for _ in 1:block.dimension, _ in 1:block.dimension
    ]
    for basis_index in eachindex(block_columns)
        basis_matrix = basis_matrices[basis_index]
        for column in axes(basis_matrix, 2), row in 1:column
            coefficient = basis_matrix[row, column]
            iszero(coefficient) && continue
            JuMP.add_to_expression!(
                matrix[row, column],
                Float64(coefficient),
                weights[basis_index],
            )
            if row != column
                JuMP.add_to_expression!(
                    matrix[column, row],
                    Float64(coefficient),
                    weights[basis_index],
                )
            end
        end
    end
    JuMP.@constraint(
        model,
        sum(matrix[index, index] for index in 1:block.dimension) == 1,
    )
    shifted = copy(matrix)
    for index in 1:block.dimension
        JuMP.add_to_expression!(shifted[index, index], -1.0, margin)
    end
    JuMP.@constraint(model, Symmetric(shifted) in JuMP.PSDCone())
    JuMP.@objective(model, Max, margin)
    JuMP.optimize!(model)

    termination = JuMP.termination_status(model)
    primal = JuMP.primal_status(model)
    println("termination_status\t", termination)
    println("primal_status\t", primal)
    println("optimizer_invoked\ttrue")
    termination in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL) ||
        error("repair search did not return an optimization candidate")
    primal in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT) ||
        error("repair search did not return a primal candidate")
    numerical_candidate_only =
        termination != MOI.OPTIMAL || primal != MOI.FEASIBLE_POINT
    println("numerical_candidate_only\t", numerical_candidate_only)
    numerical_margin = JuMP.objective_value(model)
    println("numerical_margin\t", numerical_margin)

    for rational_tolerance in rational_tolerances
        rational_weights = [
            rationalize(
                BigInt,
                JuMP.value(weights[index]);
                tol=rational_tolerance,
            )
            for index in eachindex(block_columns)
        ]
        direction = zeros(BigRational, problem.variable_count)
        for (weight, basis_direction) in
            zip(rational_weights, basis_directions)
            iszero(weight) && continue
            direction .+= weight .* basis_direction
        end
        proof = rigorous_psd_proof(block, direction; precision=256)
        println(
            "reconstruction\ttolerance\t",
            rational_tolerance,
            "\tnonzero_weights\t",
            count(!iszero, rational_weights),
            "\tstatus\t",
            proof.status,
            "\tproved\t",
            proof.proved,
            "\tproved_indefinite\t",
            proof.proved_indefinite,
        )
        flush(stdout)
        proof.proved || continue
        residuals = exact_residuals(problem, direction)
        all(iszero, residuals) ||
            error("rationalized repair direction left an affine residual")
        println(
            "direction_objective_exact\t",
            objective_improvement(problem, direction),
        )
        println("selected_rational_tolerance\t", rational_tolerance)
        write_exact_ray(output_path, direction)
        println("exact_direction_path\t", abspath(output_path))
        return 0
    end
    error("no rationalized repair direction is rigorously PSD")
end

exit(main())
