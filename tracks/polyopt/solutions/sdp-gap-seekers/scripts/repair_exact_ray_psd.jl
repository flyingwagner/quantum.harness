#!/usr/bin/env julia

using LinearAlgebra
import Clarabel
import JuMP

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

const MOI = JuMP.MOI

function main(args=ARGS)
    length(args) == 6 || error(
        "usage: repair_exact_ray_psd.jl MODEL.mof.json[.gz] BASE_EXACT_RAY.tsv " *
        "PSD_BLOCK WEIGHT_RADIUS RATIONAL_TOLERANCE[,TOLERANCE...] " *
        "REPAIRED_EXACT_RAY.tsv",
    )
    model_path, base_path, block_text, radius_text, tolerance_text, output_path =
        args
    block_index = parse(Int, block_text)
    weight_radius = parse(Float64, radius_text)
    weight_radius > 0 || error("weight radius must be positive")
    rational_tolerances = parse.(Float64, split(tolerance_text, ','))
    all(>=(0), rational_tolerances) ||
        error("rational tolerances must be nonnegative")

    problem = extract_exact_problem(model_path)
    block_index in eachindex(problem.psd_blocks) ||
        error("PSD block index is out of range")
    base = read_exact_ray_values(base_path)
    length(base) == problem.variable_count ||
        error("model/ray variable count mismatch")
    all(iszero, exact_residuals(problem, base)) ||
        error("base ray does not satisfy exact affine equalities")
    objective_improvement(problem, base) > 0 ||
        error("base ray does not improve the objective")

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

    base_matrix = GapRayPostprocess.exact_block_matrix(block, base)
    model = JuMP.Model(Clarabel.Optimizer)
    JuMP.set_silent(model)
    JuMP.set_optimizer_attribute(model, "tol_gap_abs", 1e-10)
    JuMP.set_optimizer_attribute(model, "tol_gap_rel", 1e-10)
    JuMP.set_optimizer_attribute(model, "tol_feas", 1e-10)
    JuMP.set_optimizer_attribute(model, "max_iter", 300)
    JuMP.@variable(
        model,
        -weight_radius <= weights[eachindex(block_columns)] <= weight_radius,
    )
    JuMP.@variable(model, margin)
    matrix = [
        JuMP.AffExpr(Float64(base_matrix[row, column]))
        for row in 1:block.dimension, column in 1:block.dimension
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
    println(
        "numerical_candidate_only\t",
        termination != MOI.OPTIMAL || primal != MOI.FEASIBLE_POINT,
    )
    println("numerical_margin\t", JuMP.objective_value(model))

    for rational_tolerance in rational_tolerances
        rational_weights = [
            rationalize(
                BigInt,
                JuMP.value(weights[index]);
                tol=rational_tolerance,
            )
            for index in eachindex(block_columns)
        ]
        repair = zeros(BigRational, problem.variable_count)
        for (weight, basis_direction) in
            zip(rational_weights, basis_directions)
            iszero(weight) && continue
            repair .+= weight .* basis_direction
        end
        trial = base .+ repair
        target_proof = rigorous_psd_proof(block, trial; precision=256)
        println(
            "reconstruction\ttolerance\t",
            rational_tolerance,
            "\tnonzero_weights\t",
            count(!iszero, rational_weights),
            "\ttarget_status\t",
            target_proof.status,
            "\ttarget_proved\t",
            target_proof.proved,
        )
        flush(stdout)
        target_proof.proved || continue
        all(iszero, exact_residuals(problem, trial)) ||
            error("repaired ray left an exact affine residual")
        objective_improvement(problem, trial) > 0 ||
            error("repaired ray lost objective improvement")
        proofs = [
            rigorous_psd_proof(candidate_block, trial; precision=256)
            for candidate_block in problem.psd_blocks
        ]
        for (candidate_index, proof) in enumerate(proofs)
            println(
                "block\t",
                candidate_index,
                "\tstatus\t",
                proof.status,
                "\tproved\t",
                proof.proved,
            )
        end
        all(proof -> proof.proved, proofs) || continue
        println(
            "objective_improvement_exact\t",
            objective_improvement(problem, trial),
        )
        println("selected_rational_tolerance\t", rational_tolerance)
        write_exact_ray(output_path, trial)
        println("repaired_exact_ray_path\t", abspath(output_path))
        return 0
    end
    error("no rationalized repaired ray passed every rigorous PSD block")
end

exit(main())
