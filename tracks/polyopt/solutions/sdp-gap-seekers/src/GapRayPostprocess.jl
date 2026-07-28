module GapRayPostprocess

using LinearAlgebra
import JuMP

const MOI = JuMP.MOI
const MOIU = MOI.Utilities
const BigRational = Rational{BigInt}

export BigRational,
    ExactAffineRow,
    ExactRayProblem,
    PSDDirectionBlock,
    correct_with_private_pivots,
    deduplicate_affine_equalities!,
    equality_conditioning,
    exact_residuals,
    extract_exact_problem,
    float_psd_minima,
    rigorous_psd_proof,
    normalize_rational_ray,
    objective_improvement,
    private_pivot_rows,
    psd_block_scales,
    read_ray_values,
    write_exact_ray

struct ExactAffineRow
    terms::Vector{Pair{Int,BigRational}}
end

function vector_affine_function(variables::MOI.VectorOfVariables)
    return MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(
                output_index,
                MOI.ScalarAffineTerm(1.0, variable),
            ) for (output_index, variable) in enumerate(variables.variables)
        ],
        zeros(Float64, length(variables.variables)),
    )
end

vector_affine_function(affine::MOI.VectorAffineFunction) = affine

"""
    deduplicate_affine_equalities!(model)

Remove bit-for-bit duplicate scalar equalities and duplicate coordinates of a
vector-in-`Zeros` constraint. This preserves the affine feasible set exactly;
it does not apply a tolerance, drop merely dependent rows, or alter variables
and cone constraints. The returned counts make the transformation auditable.
"""
function deduplicate_affine_equalities!(model)
    variables = sort(MOI.get(model, MOI.ListOfVariableIndices()); by=index -> index.value)
    ordinals = Dict(variable => ordinal for (ordinal, variable) in enumerate(variables))
    scalar_seen = Dict{Any,Any}()
    scalar_removed = 0
    vector_removed = 0

    for (function_type, set_type) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
        indices = collect(MOI.get(
            model,
            MOI.ListOfConstraintIndices{function_type,set_type}(),
        ))
        if set_type <: MOI.EqualTo
            for constraint in indices
                moi_function = MOI.get(model, MOI.ConstraintFunction(), constraint)
                set = MOI.get(model, MOI.ConstraintSet(), constraint)
                row = scalar_row(ordinals, moi_function)
                offset = affine_constant(moi_function) - exact_rational(set.value)
                signature = (Tuple(row.terms), offset)
                if haskey(scalar_seen, signature)
                    MOI.delete(model, constraint)
                    scalar_removed += 1
                else
                    scalar_seen[signature] = constraint
                end
            end
        elseif set_type <: MOI.Zeros
            for constraint in indices
                original = MOI.get(model, MOI.ConstraintFunction(), constraint)
                affine = vector_affine_function(original)
                rows = vector_rows(ordinals, affine)
                constants = vector_constants(affine)
                seen = Set{Any}()
                retained = Int[]
                for output_index in eachindex(rows)
                    signature = (Tuple(rows[output_index].terms), constants[output_index])
                    if signature ∉ seen
                        push!(seen, signature)
                        push!(retained, output_index)
                    end
                end
                removed = length(rows) - length(retained)
                iszero(removed) && continue
                new_output = zeros(Int, length(rows))
                for (output_index, old_output) in enumerate(retained)
                    new_output[old_output] = output_index
                end
                new_terms = [
                    MOI.VectorAffineTerm(
                        new_output[term.output_index],
                        term.scalar_term,
                    ) for term in affine.terms if new_output[term.output_index] != 0
                ]
                new_function = MOI.VectorAffineFunction(
                    new_terms,
                    affine.constants[retained],
                )
                MOI.delete(model, constraint)
                MOI.add_constraint(model, new_function, MOI.Zeros(length(retained)))
                vector_removed += removed
            end
        end
    end
    return (
        scalar_equalities_removed=scalar_removed,
        vector_zero_coordinates_removed=vector_removed,
        total_removed=scalar_removed + vector_removed,
    )
end

struct PSDDirectionBlock
    dimension::Int
    coordinates::Vector{ExactAffineRow}
end

struct ExactRayProblem
    variable_count::Int
    equalities::Vector{ExactAffineRow}
    equality_offsets::Vector{BigRational}
    psd_blocks::Vector{PSDDirectionBlock}
    objective::ExactAffineRow
    objective_sense::MOI.OptimizationSense
end

affine_constant(::MOI.VariableIndex; coefficient_tolerance::Real=0) =
    zero(BigRational)

affine_constant(
    affine::MOI.ScalarAffineFunction;
    coefficient_tolerance::Real=0,
) = exact_rational(affine.constant; tolerance=coefficient_tolerance)

vector_constants(
    variables::MOI.VectorOfVariables;
    coefficient_tolerance::Real=0,
) = fill(zero(BigRational), length(variables.variables))

vector_constants(
    affine::MOI.VectorAffineFunction;
    coefficient_tolerance::Real=0,
) = [
    exact_rational(value; tolerance=coefficient_tolerance)
    for value in affine.constants
]

exact_rational(value::Real; tolerance::Real=0) =
    rationalize(BigInt, value; tol=tolerance)

function read_ray_values(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("ray file is empty: $path")
    lines[1] == "ordinal\tmoi_index\tname\tvalue" ||
        error("unexpected ray header")
    values = Float64[]
    moi_indices = Int[]
    for (expected_ordinal, line) in enumerate(lines[2:end])
        fields = split(line, '\t'; keepempty=true)
        length(fields) == 4 || error("malformed ray row $expected_ordinal")
        parse(Int, fields[1]) == expected_ordinal ||
            error("non-canonical ray ordinal")
        push!(moi_indices, parse(Int, fields[2]))
        push!(values, parse(Float64, fields[4]))
    end
    all(isfinite, values) || error("ray contains non-finite values")
    return moi_indices, values
end

function scalar_row(
    variable_ordinals::Dict{MOI.VariableIndex,Int},
    variable::MOI.VariableIndex,
    ;
    coefficient_tolerance::Real=0,
)
    return ExactAffineRow([variable_ordinals[variable] => one(BigRational)])
end

function scalar_row(
    variable_ordinals::Dict{MOI.VariableIndex,Int},
    affine::MOI.ScalarAffineFunction,
    ;
    coefficient_tolerance::Real=0,
)
    accumulated = Dict{Int,BigRational}()
    for term in affine.terms
        ordinal = variable_ordinals[term.variable]
        accumulated[ordinal] = get(accumulated, ordinal, zero(BigRational)) +
                               exact_rational(
                                   term.coefficient;
                                   tolerance=coefficient_tolerance,
                               )
    end
    terms = sort!(
        [ordinal => coefficient for (ordinal, coefficient) in accumulated if !iszero(coefficient)];
        by=first,
    )
    return ExactAffineRow(terms)
end

function vector_rows(
    variable_ordinals::Dict{MOI.VariableIndex,Int},
    variables::MOI.VectorOfVariables,
    ;
    coefficient_tolerance::Real=0,
)
    return [
        scalar_row(
            variable_ordinals,
            variable;
            coefficient_tolerance=coefficient_tolerance,
        ) for variable in variables.variables
    ]
end

function vector_rows(
    variable_ordinals::Dict{MOI.VariableIndex,Int},
    affine::MOI.VectorAffineFunction,
    ;
    coefficient_tolerance::Real=0,
)
    rows = [Dict{Int,BigRational}() for _ in 1:MOI.output_dimension(affine)]
    for term in affine.terms
        ordinal = variable_ordinals[term.scalar_term.variable]
        coefficient = exact_rational(
            term.scalar_term.coefficient;
            tolerance=coefficient_tolerance,
        )
        row = rows[term.output_index]
        row[ordinal] = get(row, ordinal, zero(BigRational)) + coefficient
    end
    return [
        ExactAffineRow(
            sort!(
                [ordinal => coefficient for (ordinal, coefficient) in row if !iszero(coefficient)];
                by=first,
            ),
        ) for row in rows
    ]
end

function extract_exact_problem(
    model_path::AbstractString;
    coefficient_tolerance::Real=0,
)
    coefficient_tolerance >= 0 ||
        throw(ArgumentError("coefficient_tolerance must be nonnegative"))
    model = MOI.FileFormats.Model(filename=String(model_path))
    MOI.read_from_file(model, String(model_path))
    variables = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=index -> index.value,
    )
    ordinals = Dict(variable => ordinal for (ordinal, variable) in enumerate(variables))
    equalities = ExactAffineRow[]
    equality_offsets = BigRational[]
    psd_blocks = PSDDirectionBlock[]

    for (function_type, set_type) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
        indices = MOI.get(
            model,
            MOI.ListOfConstraintIndices{function_type,set_type}(),
        )
        for constraint in indices
            moi_function = MOI.get(model, MOI.ConstraintFunction(), constraint)
            set = MOI.get(model, MOI.ConstraintSet(), constraint)
            if set isa MOI.EqualTo
                push!(
                    equalities,
                    scalar_row(
                        ordinals,
                        moi_function;
                        coefficient_tolerance=coefficient_tolerance,
                    ),
                )
                push!(
                    equality_offsets,
                    affine_constant(
                        moi_function;
                        coefficient_tolerance=coefficient_tolerance,
                    ) - exact_rational(
                        set.value;
                        tolerance=coefficient_tolerance,
                    ),
                )
            elseif set isa MOI.Zeros
                append!(
                    equalities,
                    vector_rows(
                        ordinals,
                        moi_function;
                        coefficient_tolerance=coefficient_tolerance,
                    ),
                )
                append!(
                    equality_offsets,
                    vector_constants(
                        moi_function;
                        coefficient_tolerance=coefficient_tolerance,
                    ),
                )
            elseif set isa MOI.PositiveSemidefiniteConeTriangle
                rows = vector_rows(
                    ordinals,
                    moi_function;
                    coefficient_tolerance=coefficient_tolerance,
                )
                length(rows) == div(set.side_dimension * (set.side_dimension + 1), 2) ||
                    error("PSD triangle dimension mismatch")
                push!(psd_blocks, PSDDirectionBlock(set.side_dimension, rows))
            elseif set isa MOI.Reals
                nothing
            else
                error("strict post-processing does not support $(typeof(set))")
            end
        end
    end

    objective_type = MOI.get(model, MOI.ObjectiveFunctionType())
    objective_function = MOI.get(model, MOI.ObjectiveFunction{objective_type}())
    objective = scalar_row(
        ordinals,
        objective_function;
        coefficient_tolerance=coefficient_tolerance,
    )
    sense = MOI.get(model, MOI.ObjectiveSense())
    sense in (MOI.MAX_SENSE, MOI.MIN_SENSE) ||
        error("strict post-processing requires a min/max objective")
    return ExactRayProblem(
        length(variables),
        equalities,
        equality_offsets,
        psd_blocks,
        objective,
        sense,
    )
end

function evaluate(row::ExactAffineRow, values::AbstractVector)
    return sum(coefficient * values[ordinal] for (ordinal, coefficient) in row.terms;
               init=zero(BigRational))
end

exact_residuals(problem::ExactRayProblem, values::AbstractVector) =
    [evaluate(row, values) for row in problem.equalities]

"""
    equality_conditioning(problem, values; precision=256)

Measure each homogeneous equality on a floating ray. `residual_high_precision`
is evaluated from the exact rational model coefficients and the exact binary
Float64 ray coordinates at `precision` bits. `backward_error` divides that
residual by the absolute term sum, while `row_scaled_residual` additionally
divides by the row infinity norm and the ray infinity norm. These metrics
separate a genuinely inaccurate ray from Float64 summation cancellation.
"""
function equality_conditioning(
    problem::ExactRayProblem,
    values::AbstractVector{<:AbstractFloat};
    precision::Integer=256,
)
    length(values) == problem.variable_count || error("model/ray variable count mismatch")
    precision >= 64 || throw(ArgumentError("precision must be at least 64 bits"))
    ray_scale = maximum(abs, values)
    ray_scale > 0 || error("ray is identically zero")
    return setprecision(precision) do
        [
            begin
                products = [Float64(coefficient) * values[ordinal] for
                            (ordinal, coefficient) in row.terms]
                residual_float = sum(products; init=0.0)
                residual_high = sum(
                    (BigFloat(coefficient) * BigFloat(values[ordinal]) for
                     (ordinal, coefficient) in row.terms);
                    init=BigFloat(0),
                )
                absolute_term_sum = sum(abs, products; init=0.0)
                coefficient_scale = maximum(
                    (abs(Float64(coefficient)) for (_, coefficient) in row.terms);
                    init=0.0,
                )
                high_float = Float64(residual_high)
                (
                    term_count=length(row.terms),
                    residual_float=residual_float,
                    residual_high_precision=high_float,
                    summation_difference=abs(residual_float - high_float),
                    absolute_term_sum=absolute_term_sum,
                    backward_error=absolute_term_sum == 0 ? 0.0 :
                                   abs(high_float) / absolute_term_sum,
                    coefficient_scale=coefficient_scale,
                    row_scaled_residual=coefficient_scale == 0 ? 0.0 :
                                        abs(high_float) /
                                        (coefficient_scale * ray_scale),
                )
            end for row in problem.equalities
        ]
    end
end

"""Return direction-coordinate scales for each PSD triangle block."""
function psd_block_scales(problem::ExactRayProblem, values::AbstractVector)
    return [
        begin
            coordinates = [abs(Float64(evaluate(row, values))) for row in block.coordinates]
            nonzero = filter(!iszero, coordinates)
            (
                dimension=block.dimension,
                coordinate_count=length(coordinates),
                nonzero_coordinates=length(nonzero),
                maximum=maximum(coordinates; init=0.0),
                minimum_nonzero=isempty(nonzero) ? 0.0 : minimum(nonzero),
            )
        end for block in problem.psd_blocks
    ]
end

function objective_improvement(problem::ExactRayProblem, values::AbstractVector)
    value = evaluate(problem.objective, values)
    return problem.objective_sense == MOI.MAX_SENSE ? value : -value
end

function normalize_rational_ray(
    values::AbstractVector{<:Real};
    rational_tolerance::Real=1e-12,
    zero_threshold::Real=0.0,
)
    rational_tolerance >= 0 || throw(ArgumentError("rational_tolerance must be nonnegative"))
    zero_threshold >= 0 || throw(ArgumentError("zero_threshold must be nonnegative"))
    scale, pivot = findmax(abs.(values))
    scale > 0 || error("ray is identically zero")
    normalized = Vector{BigRational}(undef, length(values))
    for index in eachindex(values)
        value = values[index] / scale
        normalized[index] = abs(value) <= zero_threshold ?
                            zero(BigRational) :
                            rationalize(BigInt, value; tol=rational_tolerance)
    end
    normalized[pivot] = signbit(values[pivot]) ? -one(BigRational) : one(BigRational)
    return normalized, Float64(scale), pivot
end

function write_exact_ray(path::AbstractString, values::AbstractVector{BigRational})
    open(path, "w") do io
        println(io, "ordinal\tnumerator\tdenominator")
        for (ordinal, value) in enumerate(values)
            println(io, ordinal, '\t', numerator(value), '\t', denominator(value))
        end
    end
    return path
end

function column_occurrences(problem::ExactRayProblem)
    occurrences = zeros(Int, problem.variable_count)
    for row in problem.equalities, (ordinal, _) in row.terms
        occurrences[ordinal] += 1
    end
    return occurrences
end

function private_pivot_rows(problem::ExactRayProblem; forbidden=Set{Int}())
    occurrences = column_occurrences(problem)
    pivots = Vector{Union{Nothing,Pair{Int,BigRational}}}(undef, length(problem.equalities))
    for (row_index, row) in enumerate(problem.equalities)
        candidates = [term for term in row.terms if occurrences[first(term)] == 1 && first(term) ∉ forbidden]
        pivots[row_index] = isempty(candidates) ? nothing :
                            first(sort(candidates; by=term -> (-abs(last(term)), first(term))))
    end
    return pivots
end

function correct_with_private_pivots(
    problem::ExactRayProblem,
    values::AbstractVector{BigRational};
    forbidden=Set{Int}(),
)
    corrected = copy(values)
    residuals = exact_residuals(problem, corrected)
    pivots = private_pivot_rows(problem; forbidden=forbidden)
    unresolved = Int[]
    corrections = Pair{Int,BigRational}[]
    for row_index in eachindex(problem.equalities)
        residual = residuals[row_index]
        iszero(residual) && continue
        pivot = pivots[row_index]
        if isnothing(pivot)
            push!(unresolved, row_index)
            continue
        end
        ordinal, coefficient = something(pivot)
        delta = -residual / coefficient
        corrected[ordinal] += delta
        push!(corrections, ordinal => delta)
    end
    return corrected, unresolved, corrections
end

function evaluate_block(block::PSDDirectionBlock, values::AbstractVector)
    matrix = zeros(Float64, block.dimension, block.dimension)
    for (index, row) in enumerate(block.coordinates)
        i, j = MOIU.inverse_trimap(index)
        value = Float64(evaluate(row, values))
        matrix[i, j] = value
        matrix[j, i] = value
    end
    return matrix
end

float_psd_minima(problem::ExactRayProblem, values::AbstractVector) =
    [eigmin(Symmetric(evaluate_block(block, values))) for block in problem.psd_blocks]

function exact_block_matrix(block::PSDDirectionBlock, values::AbstractVector)
    matrix = fill(zero(BigRational), block.dimension, block.dimension)
    for (index, row) in enumerate(block.coordinates)
        i, j = MOIU.inverse_trimap(index)
        value = evaluate(row, values)
        matrix[i, j] = value
        matrix[j, i] = value
    end
    return matrix
end

struct BFInterval
    lower::BigFloat
    upper::BigFloat
    function BFInterval(lower::BigFloat, upper::BigFloat)
        lower <= upper || error("invalid interval")
        new(lower, upper)
    end
end

down(f) = setrounding(f, BigFloat, RoundDown)
up(f) = setrounding(f, BigFloat, RoundUp)

interval(value::BigRational) = BFInterval(
    BigFloat(value, RoundDown),
    BigFloat(value, RoundUp),
)

Base.:+(left::BFInterval, right::BFInterval) = BFInterval(
    down(() -> left.lower + right.lower),
    up(() -> left.upper + right.upper),
)

Base.:-(left::BFInterval, right::BFInterval) = BFInterval(
    down(() -> left.lower - right.upper),
    up(() -> left.upper - right.lower),
)

function Base.:*(left::BFInterval, right::BFInterval)
    lower_products = BigFloat[
        down(() -> left.lower * right.lower),
        down(() -> left.lower * right.upper),
        down(() -> left.upper * right.lower),
        down(() -> left.upper * right.upper),
    ]
    upper_products = BigFloat[
        up(() -> left.lower * right.lower),
        up(() -> left.lower * right.upper),
        up(() -> left.upper * right.lower),
        up(() -> left.upper * right.upper),
    ]
    return BFInterval(minimum(lower_products), maximum(upper_products))
end

function Base.:/(left::BFInterval, right::BFInterval)
    right.lower > 0 || error("interval divisor is not strictly positive")
    lower_quotients = BigFloat[
        down(() -> left.lower / right.lower),
        down(() -> left.lower / right.upper),
        down(() -> left.upper / right.lower),
        down(() -> left.upper / right.upper),
    ]
    upper_quotients = BigFloat[
        up(() -> left.lower / right.lower),
        up(() -> left.lower / right.upper),
        up(() -> left.upper / right.lower),
        up(() -> left.upper / right.upper),
    ]
    return BFInterval(minimum(lower_quotients), maximum(upper_quotients))
end

function interval_ldlt_positive_definite(matrix::Matrix{BigRational}; precision=256)
    size(matrix, 1) == size(matrix, 2) || error("matrix is not square")
    n = size(matrix, 1)
    return setprecision(precision) do
        lower = Matrix{Union{Nothing,BFInterval}}(nothing, n, n)
        diagonal = Vector{BFInterval}(undef, n)
        pivot_lower_bounds = BigFloat[]
        for k in 1:n
            pivot = interval(matrix[k, k])
            for j in 1:k-1
                lkj = something(lower[k, j])
                pivot = pivot - lkj * lkj * diagonal[j]
            end
            if pivot.lower <= 0
                return (
                    proved=false,
                    failed_pivot=k,
                    pivot_lower_bounds=pivot_lower_bounds,
                )
            end
            diagonal[k] = pivot
            push!(pivot_lower_bounds, pivot.lower)
            for i in k+1:n
                numerator = interval(matrix[i, k])
                for j in 1:k-1
                    numerator = numerator -
                                something(lower[i, j]) *
                                something(lower[k, j]) *
                                diagonal[j]
                end
                lower[i, k] = numerator / pivot
            end
        end
        return (
            proved=true,
            failed_pivot=nothing,
            pivot_lower_bounds=pivot_lower_bounds,
        )
    end
end

function rigorous_psd_proof(block::PSDDirectionBlock, values::AbstractVector; precision=256)
    matrix = exact_block_matrix(block, values)
    zero_rows = Int[]
    active_rows = Int[]
    for row in axes(matrix, 1)
        if all(iszero, @view matrix[row, :])
            push!(zero_rows, row)
        else
            push!(active_rows, row)
        end
    end
    isempty(active_rows) && return (
        proved=true,
        status="exact_zero",
        dimension=size(matrix, 1),
        active_dimension=0,
        exact_zero_rows=Tuple(zero_rows),
        failed_pivot=nothing,
        minimum_pivot_lower=Inf,
    )
    active = matrix[active_rows, active_rows]
    result = interval_ldlt_positive_definite(active; precision=precision)
    return (
        proved=result.proved,
        status=result.proved ? "interval_ldlt_positive_semidefinite" : "interval_ldlt_failed",
        dimension=size(matrix, 1),
        active_dimension=length(active_rows),
        exact_zero_rows=Tuple(zero_rows),
        failed_pivot=result.failed_pivot,
        minimum_pivot_lower=isempty(result.pivot_lower_bounds) ?
                            -Inf :
                            minimum(result.pivot_lower_bounds),
    )
end

end
