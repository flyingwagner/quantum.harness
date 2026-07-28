#!/usr/bin/env julia

using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function exact_particular_solution(
    coefficients::Matrix{BigRational},
    right_hand_side::Vector{BigRational},
)
    row_count, column_count = size(coefficients)
    length(right_hand_side) == row_count ||
        error("exact system dimension mismatch")
    augmented = hcat(coefficients, right_hand_side)
    pivot_columns = Int[]
    pivot_row = 1
    for column in 1:column_count
        candidate = findfirst(
            row -> !iszero(augmented[row, column]),
            pivot_row:row_count,
        )
        isnothing(candidate) && continue
        selected_row = pivot_row - 1 + something(candidate)
        if selected_row != pivot_row
            augmented[pivot_row, :], augmented[selected_row, :] =
                copy(augmented[selected_row, :]),
                copy(augmented[pivot_row, :])
        end
        pivot = augmented[pivot_row, column]
        augmented[pivot_row, :] ./= pivot
        for row in 1:row_count
            row == pivot_row && continue
            factor = augmented[row, column]
            iszero(factor) && continue
            augmented[row, :] .-= factor .* augmented[pivot_row, :]
        end
        push!(pivot_columns, column)
        pivot_row += 1
        pivot_row > row_count && break
    end
    for row in pivot_row:row_count
        all(iszero, @view augmented[row, 1:column_count]) ||
            error("internal exact row-reduction failure")
        iszero(augmented[row, end]) ||
            error("exact kernel-repair system is inconsistent")
    end
    solution = zeros(BigRational, column_count)
    for (row, column) in enumerate(pivot_columns)
        solution[column] = augmented[row, end]
    end
    return solution, length(pivot_columns)
end

function exact_pivot_columns(coefficients::Matrix{BigRational})
    matrix = copy(coefficients)
    row_count, column_count = size(matrix)
    pivot_columns = Int[]
    pivot_row = 1
    for column in 1:column_count
        candidate = findfirst(
            row -> !iszero(matrix[row, column]),
            pivot_row:row_count,
        )
        isnothing(candidate) && continue
        selected_row = pivot_row - 1 + something(candidate)
        if selected_row != pivot_row
            matrix[pivot_row, :], matrix[selected_row, :] =
                copy(matrix[selected_row, :]),
                copy(matrix[pivot_row, :])
        end
        pivot = matrix[pivot_row, column]
        matrix[pivot_row, :] ./= pivot
        for row in 1:row_count
            row == pivot_row && continue
            factor = matrix[row, column]
            iszero(factor) && continue
            matrix[row, :] .-= factor .* matrix[pivot_row, :]
        end
        push!(pivot_columns, column)
        pivot_row += 1
        pivot_row > row_count && break
    end
    return pivot_columns
end

function main(args=ARGS)
    length(args) in (5, 6) || error(
        "usage: repair_exact_psd_kernel.jl MODEL.mof.json[.gz] " *
        "BASE_EXACT_RAY.tsv PSD_BLOCK KERNEL_SUPPORT REPAIRED_EXACT_RAY.tsv " *
        "[COEFFICIENT_TOLERANCE]",
    )
    model_path, base_path, block_text, support_text, output_path = args[1:5]
    coefficient_tolerance =
        length(args) == 6 ? parse(Float64, args[6]) : 0.0
    block_index = parse(Int, block_text)
    kernel_supports = [
        parse.(Int, split(group, ','))
        for group in split(support_text, ';')
    ]
    all(!isempty, kernel_supports) ||
        error("kernel supports must not be empty")
    flattened_support = reduce(vcat, kernel_supports)
    length(unique(flattened_support)) == length(flattened_support) ||
        error("kernel support indices must be unique across groups")

    problem = extract_exact_problem(
        model_path;
        coefficient_tolerance=coefficient_tolerance,
    )
    block_index in eachindex(problem.psd_blocks) ||
        error("PSD block index is out of range")
    block = problem.psd_blocks[block_index]
    all(index -> index in 1:block.dimension, flattened_support) ||
        error("kernel support index is out of range")
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
    block_columns = sort!(unique(
        column
        for row in block.coordinates
        for (column, _) in row.terms
        if column in free_columns
    ))
    kernel = zeros(
        BigRational,
        block.dimension,
        length(kernel_supports),
    )
    for (kernel_index, support) in enumerate(kernel_supports)
        kernel[support, kernel_index] .= 1
    end
    base_matrix = GapRayPostprocess.exact_block_matrix(block, base)
    right_hand_side = vec(-(base_matrix * kernel))

    basis_directions = Vector{Vector{BigRational}}()
    action_columns = Vector{Vector{BigRational}}()
    for (basis_index, column) in enumerate(block_columns)
        seed = zeros(BigRational, problem.variable_count)
        seed[column] = 1
        direction, _, _ = correct_with_affine_peeling(
            problem,
            seed,
            analysis,
        )
        direction_matrix =
            GapRayPostprocess.exact_block_matrix(block, direction)
        push!(basis_directions, direction)
        push!(action_columns, vec(direction_matrix * kernel))
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
    action = hcat(action_columns...)
    independent_rows = exact_pivot_columns(Matrix(transpose(action)))
    independent_action = action[independent_rows, :]
    independent_right_hand_side = right_hand_side[independent_rows]
    normal = independent_action * transpose(independent_action)
    multipliers, normal_rank = exact_particular_solution(
        normal,
        independent_right_hand_side,
    )
    normal_rank == length(independent_rows) ||
        error("exact minimum-norm normal matrix is singular")
    weights = transpose(independent_action) * multipliers
    action * weights == right_hand_side ||
        error("exact minimum-norm repair does not solve all kernel equations")
    action_rank = length(independent_rows)
    repair = zeros(BigRational, problem.variable_count)
    for (weight, basis_direction) in zip(weights, basis_directions)
        iszero(weight) && continue
        repair .+= weight .* basis_direction
    end
    trial = base .+ repair
    trial_matrix = GapRayPostprocess.exact_block_matrix(block, trial)
    all(iszero, trial_matrix * kernel) ||
        error("repaired target block does not have the declared exact kernel")
    all(iszero, exact_residuals(problem, trial)) ||
        error("kernel-repaired ray left an exact affine residual")
    improvement = objective_improvement(problem, trial)
    improvement > 0 || error("kernel-repaired ray lost objective improvement")

    removed_indices = last.(kernel_supports)
    target_proof = rigorous_psd_proof_with_exact_kernel(
        trial_matrix,
        kernel,
        removed_indices;
        precision=256,
    )
    println("kernel_action_rank\t", action_rank)
    println("nonzero_repair_weights\t", count(!iszero, weights))
    println(
        "target_reduced_float_minimum\t",
        minimum(eigvals(Symmetric(
            Float64.(trial_matrix[
                setdiff(1:block.dimension, removed_indices),
                setdiff(1:block.dimension, removed_indices),
            ]),
        ))),
    )
    println("target_reduced_proved\t", target_proof.proved)
    println("target_reduced_failed_pivot\t", target_proof.failed_pivot)
    println(
        "target_reduced_failed_pivot_upper\t",
        target_proof.failed_pivot_upper,
    )
    flush(stdout)
    target_proof.proved ||
        error("kernel-reduced target block is not rigorously positive definite")
    proofs = [
        candidate_index == block_index ?
        nothing :
        rigorous_psd_proof(candidate_block, trial; precision=256)
        for (candidate_index, candidate_block) in
            enumerate(problem.psd_blocks)
    ]
    all(
        proof -> isnothing(proof) || something(proof).proved,
        proofs,
    ) || error("kernel repair broke another PSD block")

    println("free_block_columns\t", length(block_columns))
    println(
        "kernel_supports\t",
        join((join(support, ',') for support in kernel_supports), ';'),
    )
    println("removed_kernel_coordinates\t", join(removed_indices, ','))
    println(
        "target_status\t",
        target_proof.status,
    )
    println(
        "target_minimum_reduced_pivot_lower\t",
        minimum(target_proof.pivot_lower_bounds),
    )
    for (candidate_index, proof) in enumerate(proofs)
        isnothing(proof) && continue
        println(
            "block\t",
            candidate_index,
            "\tstatus\t",
            something(proof).status,
            "\tproved\t",
            something(proof).proved,
        )
    end
    println("objective_improvement_exact\t", improvement)
    println("coefficient_tolerance\t", coefficient_tolerance)
    println("optimizer_invoked\tfalse")
    write_exact_ray(output_path, trial)
    println("repaired_exact_ray_path\t", abspath(output_path))
    return 0
end

exit(main())
