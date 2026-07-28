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
    apply_variable_scaling!,
    backtransform_ray,
    column_equilibration,
    correct_with_private_pivots,
    deduplicate_affine_equalities!,
    equality_conditioning,
    equilibrate_rows!,
    exact_residuals,
    extract_exact_problem,
    float_psd_minima,
    rigorous_psd_proof,
    normalize_rational_ray,
    normalize_recession_problem!,
    objective_improvement,
    private_pivot_rows,
    psd_block_scales,
    ray_equilibration,
    ray_congruence_equilibration,
    read_ray_values,
    read_scale_map,
    transform_ray,
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

function update_maximum!(
    maxima::Dict{MOI.VariableIndex,Float64},
    variable::MOI.VariableIndex,
    coefficient::Real,
)
    value = abs(Float64(coefficient))
    isfinite(value) || error("model contains a non-finite coefficient")
    maxima[variable] = max(get(maxima, variable, 0.0), value)
    return maxima
end

function scan_coefficients!(maxima, function_value::MOI.ScalarAffineFunction)
    for term in function_value.terms
        update_maximum!(maxima, term.variable, term.coefficient)
    end
    return maxima
end

function scan_coefficients!(maxima, function_value::MOI.VectorAffineFunction)
    for term in function_value.terms
        update_maximum!(
            maxima,
            term.scalar_term.variable,
            term.scalar_term.coefficient,
        )
    end
    return maxima
end

function scan_coefficients!(maxima, variable::MOI.VariableIndex)
    update_maximum!(maxima, variable, 1.0)
end

function scan_coefficients!(maxima, variables::MOI.VectorOfVariables)
    for variable in variables.variables
        update_maximum!(maxima, variable, 1.0)
    end
    return maxima
end

function direct_psd_groups(model)
    groups = Vector{Vector{MOI.VariableIndex}}()
    owner = Dict{MOI.VariableIndex,Int}()
    for (function_type, set_type) in MOI.get(
        model,
        MOI.ListOfConstraintTypesPresent(),
    )
        set_type <: MOI.PositiveSemidefiniteConeTriangle || continue
        constraints = MOI.get(
            model,
            MOI.ListOfConstraintIndices{function_type,set_type}(),
        )
        for constraint in constraints
            function_value = MOI.get(
                model,
                MOI.ConstraintFunction(),
                constraint,
            )
            function_value isa MOI.VectorOfVariables || continue
            group = copy(function_value.variables)
            for variable in group
                haskey(owner, variable) &&
                    error("one variable belongs to multiple direct PSD blocks")
            end
            push!(groups, group)
            for variable in group
                owner[variable] = length(groups)
            end
        end
    end
    return groups
end

function enforce_uniform_group_scales!(
    scales,
    exponents,
    variables,
    groups,
    target_values,
    exponent_sign::Int,
)
    positions = Dict(
        variable => position for (position, variable) in enumerate(variables)
    )
    for group in groups
        group_positions = [positions[variable] for variable in group]
        target = maximum(target_values[position] for position in group_positions)
        if iszero(target)
            exponent = 0
        else
            exponent = exponent_sign * floor(Int, log2(target))
        end
        scale = ldexp(1.0, exponent)
        isfinite(scale) && scale > 0 ||
            error("direct PSD block scale is not finite and positive")
        for position in group_positions
            scales[position] = scale
            exponents[position] = exponent
        end
    end
    return nothing
end

"""
Choose an invertible power-of-two substitution `x_i = scale_i * z_i` so the
largest absolute affine/objective coefficient in each transformed column lies
in `[1,2)`. Variable-domain `Reals` constraints are scale invariant and do not
participate in the maxima.
"""
function column_equilibration(model)
    variables = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    maxima = Dict(variable => 0.0 for variable in variables)
    for (function_type, set_type) in MOI.get(
        model,
        MOI.ListOfConstraintTypesPresent(),
    )
        constraints = MOI.get(
            model,
            MOI.ListOfConstraintIndices{function_type,set_type}(),
        )
        for constraint in constraints
            function_value = MOI.get(
                model,
                MOI.ConstraintFunction(),
                constraint,
            )
            if set_type <: MOI.Reals &&
               (function_value isa MOI.VariableIndex ||
                function_value isa MOI.VectorOfVariables)
                continue
            end
            function_value isa Union{
                MOI.ScalarAffineFunction,
                MOI.VectorAffineFunction,
                MOI.VectorOfVariables,
            } || error(
                "column equilibration does not support $(typeof(function_value)) in $(set_type)",
            )
            scan_coefficients!(maxima, function_value)
        end
    end
    objective_type = MOI.get(model, MOI.ObjectiveFunctionType())
    objective = MOI.get(model, MOI.ObjectiveFunction{objective_type}())
    objective isa Union{MOI.VariableIndex,MOI.ScalarAffineFunction} ||
        error("column equilibration requires a scalar affine objective")
    scan_coefficients!(maxima, objective)

    scales = Float64[]
    exponents = Int[]
    for variable in variables
        maximum_value = maxima[variable]
        if iszero(maximum_value)
            push!(scales, 1.0)
            push!(exponents, 0)
            continue
        end
        exponent = -floor(Int, log2(maximum_value))
        scale = ldexp(1.0, exponent)
        isfinite(scale) && scale > 0 ||
            error("column equilibration scale is not finite and positive")
        transformed = maximum_value * scale
        1.0 <= transformed < 2.0 ||
            error("power-of-two column equilibration invariant failed")
        push!(scales, scale)
        push!(exponents, exponent)
    end
    groups = direct_psd_groups(model)
    # A direct PSD matrix variable may only be scaled uniformly by a positive
    # scalar. Choose that scalar from the largest column coefficient in the
    # whole block, then drop it from the cone constraint during rendering.
    enforce_uniform_group_scales!(
        scales,
        exponents,
        variables,
        groups,
        [maxima[variable] for variable in variables],
        -1,
    )
    return (
        variables=variables,
        maxima=[maxima[variable] for variable in variables],
        exponents=exponents,
        scales=scales,
        direct_psd_groups=groups,
    )
end

"""
Choose positive power-of-two scales from a reference ray. Direct PSD variables
receive one uniform scale per matrix block; all other variables use their own
nonzero magnitude. This is a formulation experiment only: the candidate must
be back-transformed and replayed against the original model.
"""
function ray_equilibration(model, values::AbstractVector{<:Real})
    variables = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    length(values) == length(variables) ||
        throw(ArgumentError("reference ray length differs from model"))
    all(isfinite, values) || error("reference ray contains non-finite values")
    magnitudes = Float64.(abs.(values))
    exponents = Int[]
    scales = Float64[]
    for magnitude in magnitudes
        exponent = iszero(magnitude) ? 0 : floor(Int, log2(magnitude))
        scale = ldexp(1.0, exponent)
        isfinite(scale) && scale > 0 ||
            error("reference-ray scale is not finite and positive")
        push!(exponents, exponent)
        push!(scales, scale)
    end
    groups = direct_psd_groups(model)
    enforce_uniform_group_scales!(
        scales,
        exponents,
        variables,
        groups,
        magnitudes,
        1,
    )
    return (
        variables=variables,
        maxima=magnitudes,
        exponents=exponents,
        scales=scales,
        direct_psd_groups=groups,
    )
end

function congruence_coordinate_scales(
    group::AbstractVector{MOI.VariableIndex},
    magnitudes::Dict{MOI.VariableIndex,Float64},
)
    dimension = MOIU.side_dimension_for_vectorized_dimension(length(group))
    div(dimension * (dimension + 1), 2) == length(group) ||
        error("invalid triangular PSD block dimension")
    diagonal_exponents = Int[]
    for index in 1:dimension
        variable = group[MOIU.trimap(index, index)]
        magnitude = magnitudes[variable]
        # x_ii = d_i^2 z_ii. Rounding half the binary exponent keeps each
        # nonzero reference diagonal within a factor two of one in z.
        exponent = iszero(magnitude) ? 0 :
                   round(Int, log2(magnitude) / 2)
        push!(diagonal_exponents, exponent)
    end
    exponents = Int[]
    scales = Float64[]
    for coordinate in eachindex(group)
        row, column = MOIU.inverse_trimap(coordinate)
        exponent = diagonal_exponents[row] + diagonal_exponents[column]
        scale = ldexp(1.0, exponent)
        isfinite(scale) && scale > 0 ||
            error("PSD congruence scale is not finite and positive")
        push!(exponents, exponent)
        push!(scales, scale)
    end
    return exponents, scales
end

"""
Choose an invertible diagonal-congruence substitution for each direct PSD
matrix from a reference ray:

    X = D * Z * D

Every diagonal entry of `D` is a positive (possibly square-root)
power-of-two, so each stored triangular coordinate still has a positive
power-of-two scale. Variables outside direct PSD blocks retain independent
reference-ray scaling. Unlike arbitrary coordinate scaling, diagonal
congruence preserves positive semidefiniteness in both directions.
"""
function ray_congruence_equilibration(
    model,
    values::AbstractVector{<:Real},
)
    variables = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    length(values) == length(variables) ||
        throw(ArgumentError("reference ray length differs from model"))
    all(isfinite, values) || error("reference ray contains non-finite values")
    magnitudes_vector = Float64.(abs.(values))
    magnitudes = Dict(
        variable => magnitudes_vector[index]
        for (index, variable) in enumerate(variables)
    )
    exponents = [
        iszero(magnitude) ? 0 : floor(Int, log2(magnitude))
        for magnitude in magnitudes_vector
    ]
    scales = [
        ldexp(1.0, exponent)
        for exponent in exponents
    ]
    all(scale -> isfinite(scale) && scale > 0, scales) ||
        error("reference-ray scale is not finite and positive")
    groups = direct_psd_groups(model)
    positions = Dict(
        variable => position for (position, variable) in enumerate(variables)
    )
    for group in groups
        group_exponents, group_scales =
            congruence_coordinate_scales(group, magnitudes)
        for (coordinate, variable) in enumerate(group)
            position = positions[variable]
            exponents[position] = group_exponents[coordinate]
            scales[position] = group_scales[coordinate]
        end
    end
    return (
        variables=variables,
        maxima=magnitudes_vector,
        exponents=exponents,
        scales=scales,
        direct_psd_groups=groups,
    )
end

function power_of_two_exponent(value::Float64)
    isfinite(value) && value > 0 ||
        error("variable scale is not finite and positive")
    significand(value) == 1.0 ||
        error("direct PSD coordinate scale is not a power of two")
    return exponent(value)
end

function validate_congruence_scales(block_scales::AbstractVector{Float64})
    dimension =
        MOIU.side_dimension_for_vectorized_dimension(length(block_scales))
    div(dimension * (dimension + 1), 2) == length(block_scales) ||
        error("invalid triangular PSD block dimension")
    exponents = power_of_two_exponent.(block_scales)
    diagonal_exponents = [
        exponents[MOIU.trimap(index, index)]
        for index in 1:dimension
    ]
    for coordinate in eachindex(exponents)
        row, column = MOIU.inverse_trimap(coordinate)
        2 * exponents[coordinate] ==
            diagonal_exponents[row] + diagonal_exponents[column] ||
            error(
                "direct PSD coordinate scales are not one diagonal congruence",
            )
    end
    return nothing
end

function zero_constant(
    function_value::MOI.ScalarAffineFunction{Float64},
)
    return MOI.ScalarAffineFunction(copy(function_value.terms), 0.0)
end

function zero_constants(
    function_value::MOI.VectorAffineFunction{Float64},
)
    return MOI.VectorAffineFunction(
        copy(function_value.terms),
        zeros(Float64, length(function_value.constants)),
    )
end

function scalar_affine_objective(
    objective::MOI.VariableIndex,
)
    return MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(1.0, objective)],
        0.0,
    )
end

function scalar_affine_objective(
    objective::MOI.ScalarAffineFunction{Float64},
)
    return objective
end

"""
Replace an affine conic optimization problem by its homogeneous recession
system and fix the improving objective direction to unit magnitude.

For maximization the added normalization is `c'x = 1`; for minimization it is
`c'x = -1`. Thus the returned feasibility problem has a solution exactly when
the original recession cone contains an improving direction. Affine
constants and equality right-hand sides are removed, cone directions are
preserved, and the objective becomes feasibility. Unsupported constraint
types fail closed.
"""
function normalize_recession_problem!(model)
    original_sense = MOI.get(model, MOI.ObjectiveSense())
    original_sense in (MOI.MAX_SENSE, MOI.MIN_SENSE) ||
        error("recession normalization requires min or max objective sense")
    objective_type = MOI.get(model, MOI.ObjectiveFunctionType())
    objective = scalar_affine_objective(
        MOI.get(model, MOI.ObjectiveFunction{objective_type}()),
    )
    isempty(objective.terms) &&
        error("cannot normalize a zero improving objective")
    all(term -> isfinite(term.coefficient), objective.terms) ||
        error("objective contains a non-finite coefficient")

    scalar_equalities = 0
    zero_coordinates = 0
    affine_psd_blocks = 0
    direct_psd_blocks = 0
    for (function_type, set_type) in MOI.get(
        model,
        MOI.ListOfConstraintTypesPresent(),
    )
        constraints = collect(MOI.get(
            model,
            MOI.ListOfConstraintIndices{function_type,set_type}(),
        ))
        for constraint in constraints
            function_value = MOI.get(
                model,
                MOI.ConstraintFunction(),
                constraint,
            )
            set = MOI.get(model, MOI.ConstraintSet(), constraint)
            if set isa MOI.EqualTo &&
               function_value isa MOI.ScalarAffineFunction{Float64}
                MOI.set(
                    model,
                    MOI.ConstraintFunction(),
                    constraint,
                    zero_constant(function_value),
                )
                MOI.set(
                    model,
                    MOI.ConstraintSet(),
                    constraint,
                    MOI.EqualTo(0.0),
                )
                scalar_equalities += 1
            elseif set isa MOI.Zeros &&
                   function_value isa MOI.VectorAffineFunction{Float64}
                MOI.set(
                    model,
                    MOI.ConstraintFunction(),
                    constraint,
                    zero_constants(function_value),
                )
                zero_coordinates += length(function_value.constants)
            elseif set isa MOI.Zeros &&
                   function_value isa MOI.VectorOfVariables
                zero_coordinates += length(function_value.variables)
            elseif set isa MOI.PositiveSemidefiniteConeTriangle &&
                   function_value isa MOI.VectorAffineFunction{Float64}
                MOI.set(
                    model,
                    MOI.ConstraintFunction(),
                    constraint,
                    zero_constants(function_value),
                )
                affine_psd_blocks += 1
            elseif set isa MOI.PositiveSemidefiniteConeTriangle &&
                   function_value isa MOI.VectorOfVariables
                direct_psd_blocks += 1
            elseif set isa MOI.Reals
                nothing
            else
                error(
                    "recession normalization does not support $(typeof(function_value)) in $(typeof(set))",
                )
            end
        end
    end

    orientation = original_sense == MOI.MAX_SENSE ? 1.0 : -1.0
    normalization = MOI.ScalarAffineFunction(
        [
            MOI.ScalarAffineTerm(
                orientation * term.coefficient,
                term.variable,
            )
            for term in objective.terms
        ],
        0.0,
    )
    normalization_constraint =
        MOI.add_constraint(model, normalization, MOI.EqualTo(1.0))
    zero_objective = MOI.ScalarAffineFunction{Float64}(
        MOI.ScalarAffineTerm{Float64}[],
        0.0,
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{typeof(zero_objective)}(),
        zero_objective,
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    return (
        original_sense=original_sense,
        objective_terms=length(objective.terms),
        scalar_equalities=scalar_equalities,
        zero_coordinates=zero_coordinates,
        affine_psd_blocks=affine_psd_blocks,
        direct_psd_blocks=direct_psd_blocks,
        normalization_constraint=normalization_constraint,
    )
end

function scaled_function(
    function_value::MOI.ScalarAffineFunction{Float64},
    scales::Dict{MOI.VariableIndex,Float64},
)
    return MOI.ScalarAffineFunction(
        [
            MOI.ScalarAffineTerm(
                term.coefficient * scales[term.variable],
                term.variable,
            )
            for term in function_value.terms
        ],
        function_value.constant,
    )
end

function scaled_function(
    function_value::MOI.VectorAffineFunction{Float64},
    scales::Dict{MOI.VariableIndex,Float64},
)
    return MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(
                term.output_index,
                MOI.ScalarAffineTerm(
                    term.scalar_term.coefficient *
                        scales[term.scalar_term.variable],
                    term.scalar_term.variable,
                ),
            )
            for term in function_value.terms
        ],
        copy(function_value.constants),
    )
end

"""
Apply `x_i = scale_i*z_i` in place to every affine constraint and scalar
objective. The caller must preserve the returned map and back-transform a
candidate ray before replay against the original model.
"""
function apply_variable_scaling!(model, equilibration=column_equilibration(model))
    scale_by_variable = Dict(
        variable => equilibration.scales[index]
        for (index, variable) in enumerate(equilibration.variables)
    )
    for (function_type, set_type) in MOI.get(
        model,
        MOI.ListOfConstraintTypesPresent(),
    )
        constraints = collect(MOI.get(
            model,
            MOI.ListOfConstraintIndices{function_type,set_type}(),
        ))
        for constraint in constraints
            function_value = MOI.get(
                model,
                MOI.ConstraintFunction(),
                constraint,
            )
            if set_type <: MOI.Reals &&
               (function_value isa MOI.VariableIndex ||
                function_value isa MOI.VectorOfVariables)
                continue
            end
            if set_type <: MOI.PositiveSemidefiniteConeTriangle &&
               function_value isa MOI.VectorOfVariables
                block_scales = [
                    scale_by_variable[variable]
                    for variable in function_value.variables
                ]
                validate_congruence_scales(block_scales)
                # X = D*Z*D with invertible positive diagonal D gives X PSD
                # iff Z PSD. The congruence factors are therefore omitted from
                # this cone constraint while they remain in every affine
                # occurrence of the variables.
                continue
            end
            function_value isa Union{
                MOI.ScalarAffineFunction{Float64},
                MOI.VectorAffineFunction{Float64},
            } || error(
                "cannot scale constraint function $(typeof(function_value))",
            )
            MOI.set(
                model,
                MOI.ConstraintFunction(),
                constraint,
                scaled_function(function_value, scale_by_variable),
            )
        end
    end
    objective_type = MOI.get(model, MOI.ObjectiveFunctionType())
    objective = MOI.get(model, MOI.ObjectiveFunction{objective_type}())
    if objective isa MOI.VariableIndex
        objective = MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(scale_by_variable[objective], objective)],
            0.0,
        )
    elseif objective isa MOI.ScalarAffineFunction{Float64}
        objective = scaled_function(objective, scale_by_variable)
    else
        error("cannot scale objective function $(typeof(objective))")
    end
    MOI.set(model, MOI.ObjectiveFunction{typeof(objective)}(), objective)
    return equilibration
end

function backtransform_ray(
    scaled_values::AbstractVector{<:Real},
    scales::AbstractVector{<:Real},
)
    length(scaled_values) == length(scales) ||
        throw(ArgumentError("scaled ray and scale map lengths differ"))
    values = Float64[
        Float64(scale) * Float64(value)
        for (scale, value) in zip(scales, scaled_values)
    ]
    all(isfinite, values) || error("back-transformed ray is non-finite")
    return values
end

function transform_ray(
    original_values::AbstractVector{<:Real},
    scales::AbstractVector{<:Real},
)
    length(original_values) == length(scales) ||
        throw(ArgumentError("original ray and scale map lengths differ"))
    values = Float64[
        Float64(value) / Float64(scale)
        for (scale, value) in zip(scales, original_values)
    ]
    all(isfinite, values) || error("transformed ray is non-finite")
    return values
end

function reciprocal_power_of_two(maximum_value::Real)
    value = abs(Float64(maximum_value))
    isfinite(value) || error("row contains a non-finite magnitude")
    iszero(value) && return 1.0, 0
    exponent = -floor(Int, log2(value))
    factor = ldexp(1.0, exponent)
    transformed = value * factor
    isfinite(factor) && factor > 0 ||
        error("row equilibration factor is not finite and positive")
    1.0 <= transformed < 2.0 ||
        error("row equilibration invariant failed")
    return factor, exponent
end

function scalar_maximum(function_value::MOI.ScalarAffineFunction, set)
    return maximum(
        (
            abs(function_value.constant),
            abs(set.value),
            (abs(term.coefficient) for term in function_value.terms)...,
        );
        init=0.0,
    )
end

function scale_scalar_function(function_value, factor)
    return MOI.ScalarAffineFunction(
        [
            MOI.ScalarAffineTerm(
                factor * term.coefficient,
                term.variable,
            )
            for term in function_value.terms
        ],
        factor * function_value.constant,
    )
end

function coordinate_factors(function_value::MOI.VectorAffineFunction)
    maxima = abs.(function_value.constants)
    for term in function_value.terms
        maxima[term.output_index] = max(
            maxima[term.output_index],
            abs(term.scalar_term.coefficient),
        )
    end
    factors = Float64[]
    exponents = Int[]
    for maximum_value in maxima
        factor, exponent = reciprocal_power_of_two(maximum_value)
        push!(factors, factor)
        push!(exponents, exponent)
    end
    return factors, exponents
end

function scale_vector_coordinates(function_value, factors)
    return MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(
                term.output_index,
                MOI.ScalarAffineTerm(
                    factors[term.output_index] *
                        term.scalar_term.coefficient,
                    term.scalar_term.variable,
                ),
            )
            for term in function_value.terms
        ],
        [
            factors[index] * function_value.constants[index]
            for index in eachindex(function_value.constants)
        ],
    )
end

"""
Apply positive power-of-two row scaling after a variable substitution:

- each scalar equality and both sides by one factor;
- each coordinate of a vector-in-`Zeros` equality by its own factor;
- each affine PSD matrix by one common factor;
- the scalar objective by one positive factor.

These operations preserve the feasible set and min/max improvement direction.
"""
function equilibrate_rows!(model)
    scalar_equalities = 0
    zero_coordinates = 0
    affine_psd_blocks = 0
    exponents = Int[]
    for (function_type, set_type) in MOI.get(
        model,
        MOI.ListOfConstraintTypesPresent(),
    )
        constraints = collect(MOI.get(
            model,
            MOI.ListOfConstraintIndices{function_type,set_type}(),
        ))
        for constraint in constraints
            function_value = MOI.get(
                model,
                MOI.ConstraintFunction(),
                constraint,
            )
            set = MOI.get(model, MOI.ConstraintSet(), constraint)
            if set isa MOI.EqualTo &&
               function_value isa MOI.ScalarAffineFunction{Float64}
                factor, exponent =
                    reciprocal_power_of_two(scalar_maximum(function_value, set))
                MOI.set(
                    model,
                    MOI.ConstraintFunction(),
                    constraint,
                    scale_scalar_function(function_value, factor),
                )
                MOI.set(
                    model,
                    MOI.ConstraintSet(),
                    constraint,
                    MOI.EqualTo(factor * set.value),
                )
                push!(exponents, exponent)
                scalar_equalities += 1
            elseif set isa MOI.Zeros &&
                   function_value isa MOI.VectorAffineFunction{Float64}
                factors, coordinate_exponents =
                    coordinate_factors(function_value)
                MOI.set(
                    model,
                    MOI.ConstraintFunction(),
                    constraint,
                    scale_vector_coordinates(function_value, factors),
                )
                append!(exponents, coordinate_exponents)
                zero_coordinates += length(factors)
            elseif set isa MOI.PositiveSemidefiniteConeTriangle &&
                   function_value isa MOI.VectorAffineFunction{Float64}
                maximum_value = maximum(
                    (
                        maximum(abs, function_value.constants; init=0.0),
                        maximum(
                            (
                                abs(term.scalar_term.coefficient)
                                for term in function_value.terms
                            );
                            init=0.0,
                        ),
                    ),
                )
                factor, exponent =
                    reciprocal_power_of_two(maximum_value)
                factors = fill(factor, length(function_value.constants))
                MOI.set(
                    model,
                    MOI.ConstraintFunction(),
                    constraint,
                    scale_vector_coordinates(function_value, factors),
                )
                push!(exponents, exponent)
                affine_psd_blocks += 1
            elseif set isa MOI.PositiveSemidefiniteConeTriangle &&
                   function_value isa MOI.VectorOfVariables
                nothing
            elseif set isa MOI.Reals
                nothing
            else
                error(
                    "row equilibration does not support $(typeof(function_value)) in $(typeof(set))",
                )
            end
        end
    end

    objective_type = MOI.get(model, MOI.ObjectiveFunctionType())
    objective = MOI.get(model, MOI.ObjectiveFunction{objective_type}())
    objective isa MOI.ScalarAffineFunction{Float64} ||
        error("row equilibration requires a scalar affine objective")
    objective_maximum = maximum(
        (
            abs(objective.constant),
            (abs(term.coefficient) for term in objective.terms)...,
        );
        init=0.0,
    )
    objective_factor, objective_exponent =
        reciprocal_power_of_two(objective_maximum)
    MOI.set(
        model,
        MOI.ObjectiveFunction{typeof(objective)}(),
        scale_scalar_function(objective, objective_factor),
    )
    return (
        scalar_equalities=scalar_equalities,
        zero_coordinates=zero_coordinates,
        affine_psd_blocks=affine_psd_blocks,
        minimum_exponent=isempty(exponents) ? 0 : minimum(exponents),
        maximum_exponent=isempty(exponents) ? 0 : maximum(exponents),
        objective_exponent=objective_exponent,
        objective_factor=objective_factor,
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

function read_scale_map(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("scale map is empty")
    lines[1] ==
        "ordinal\tmoi_index\tname\tcolumn_max\texponent\tx_per_z" ||
        error("unexpected scale-map header")
    scales = Float64[]
    indices = Int[]
    names = String[]
    for (ordinal, line) in enumerate(lines[2:end])
        fields = split(line, '\t'; keepempty=true)
        length(fields) == 6 || error("malformed scale-map row")
        parse(Int, fields[1]) == ordinal ||
            error("noncanonical scale-map ordinal")
        push!(indices, parse(Int, fields[2]))
        push!(names, fields[3])
        push!(scales, parse(Float64, fields[6]))
    end
    all(isfinite, scales) && all(>(0), scales) ||
        error("scale map contains an invalid scale")
    return indices, names, scales
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
