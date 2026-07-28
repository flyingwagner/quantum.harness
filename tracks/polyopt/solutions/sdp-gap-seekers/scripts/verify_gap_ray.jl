#!/usr/bin/env julia

module GapRayVerifier

using LinearAlgebra
import JuMP
const MOI = JuMP.MOI
const MOIU = MOI.Utilities

function read_ray(path::AbstractString)
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
    return moi_indices, values
end

function triangle_matrix(values::AbstractVector, side_dimension::Integer)
    length(values) == div(side_dimension * (side_dimension + 1), 2) ||
        error("PSD triangle dimension mismatch")
    matrix = zeros(Float64, side_dimension, side_dimension)
    for (index, value) in enumerate(values)
        row, column = MOIU.inverse_trimap(index)
        matrix[row, column] = value
        matrix[column, row] = value
    end
    return matrix
end

function direction_value(value_function, _, variable::MOI.VariableIndex)
    return value_function(variable)
end

function direction_value(value_function, _, variables::MOI.VectorOfVariables)
    return value_function.(variables.variables)
end

function direction_value(value_function, _, affine::MOI.ScalarAffineFunction)
    # Do not evaluate the affine function and subtract its value at the origin:
    # that loses low-order bits when a large ray is combined with a nonzero
    # constant. A recession direction contains the linear part only.
    return sum(
        term.coefficient * value_function(term.variable)
        for term in affine.terms;
        init=0.0,
    )
end

function direction_value(value_function, _, affine::MOI.VectorAffineFunction)
    values = zeros(Float64, MOI.output_dimension(affine))
    for term in affine.terms
        values[term.output_index] +=
            term.scalar_term.coefficient *
            value_function(term.scalar_term.variable)
    end
    return values
end


function direction_value(_, _, moi_function)
    error("unsupported non-affine ray function: $(typeof(moi_function))")
end

function verify(
    model_path::AbstractString,
    ray_path::AbstractString;
    absolute_tolerance::Real=1e-12,
    relative_tolerance::Real=1e-12,
)
    isfinite(absolute_tolerance) && absolute_tolerance >= 0 ||
        throw(ArgumentError("absolute_tolerance must be finite and nonnegative"))
    isfinite(relative_tolerance) && relative_tolerance >= 0 ||
        throw(ArgumentError("relative_tolerance must be finite and nonnegative"))
    model = MOI.FileFormats.Model(filename=String(model_path))
    MOI.read_from_file(model, String(model_path))
    variables = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=index -> index.value,
    )
    moi_indices, ray_values = read_ray(ray_path)
    length(variables) == length(ray_values) ||
        error("model/ray variable count mismatch")
    all(isfinite, ray_values) || error("ray contains non-finite values")
    for (ordinal, variable) in enumerate(variables)
        variable.value == moi_indices[ordinal] ||
            error("model/ray MOI index mismatch at ordinal $ordinal")
        MOI.get(model, MOI.VariableName(), variable) == "C$ordinal" ||
            error("model/ray generic-name mismatch at ordinal $ordinal")
    end
    ray = Dict(variable => ray_values[index] for (index, variable) in enumerate(variables))
    value_function(variable) = ray[variable]
    scale = maximum(abs, ray_values)
    scale > 0.0 || error("ray is identically zero")
    normalized_tolerance =
        Float64(absolute_tolerance) / scale + Float64(relative_tolerance)

    equality_count = 0
    psd_count = 0
    scalar_cone_count = 0
    max_equality_residual = 0.0
    max_scalar_cone_violation = 0.0
    min_psd_eigenvalue = Inf
    unsupported_sets = String[]

    for (function_type, set_type) in
        MOI.get(model, MOI.ListOfConstraintTypesPresent())
        indices = MOI.get(
            model,
            MOI.ListOfConstraintIndices{function_type, set_type}(),
        )
        for constraint in indices
            function_value = direction_value(
                value_function,
                model,
                MOI.get(model, MOI.ConstraintFunction(), constraint),
            )
            set = MOI.get(model, MOI.ConstraintSet(), constraint)
            if set isa MOI.EqualTo
                equality_count += 1
                max_equality_residual = max(
                    max_equality_residual,
                    abs(function_value),
                )
            elseif set isa MOI.Zeros
                equality_count += length(function_value)
                max_equality_residual = max(
                    max_equality_residual,
                    maximum(abs, function_value; init=0.0),
                )
            elseif set isa MOI.GreaterThan
                scalar_cone_count += 1
                max_scalar_cone_violation = max(
                    max_scalar_cone_violation,
                    max(0.0, -function_value),
                )
            elseif set isa MOI.LessThan
                scalar_cone_count += 1
                max_scalar_cone_violation = max(
                    max_scalar_cone_violation,
                    max(0.0, function_value),
                )
            elseif set isa MOI.Interval
                scalar_cone_count += 1
                max_scalar_cone_violation = max(
                    max_scalar_cone_violation,
                    abs(function_value),
                )
            elseif set isa MOI.Nonnegatives
                scalar_cone_count += length(function_value)
                max_scalar_cone_violation = max(
                    max_scalar_cone_violation,
                    max(0.0, -minimum(function_value; init=0.0)),
                )
            elseif set isa MOI.Nonpositives
                scalar_cone_count += length(function_value)
                max_scalar_cone_violation = max(
                    max_scalar_cone_violation,
                    max(0.0, maximum(function_value; init=0.0)),
                )
            elseif set isa MOI.PositiveSemidefiniteConeTriangle
                psd_count += 1
                matrix = triangle_matrix(function_value, set.side_dimension)
                min_psd_eigenvalue = min(
                    min_psd_eigenvalue,
                    eigmin(Symmetric(matrix)),
                )
            elseif set isa MOI.Reals
                scalar_cone_count += length(function_value)
            else
                push!(unsupported_sets, string(typeof(set)))
            end
        end
    end
    isempty(unsupported_sets) ||
        error("unsupported constraint sets: $(join(unique(unsupported_sets), ", "))")

    objective_type = MOI.get(model, MOI.ObjectiveFunctionType())
    objective = MOI.get(model, MOI.ObjectiveFunction{objective_type}())
    objective_value = direction_value(value_function, model, objective)
    objective_sense = MOI.get(model, MOI.ObjectiveSense())
    objective_sense in (MOI.MAX_SENSE, MOI.MIN_SENSE) ||
        error("ray verification requires a maximization or minimization objective")
    improving_objective =
        objective_sense == MOI.MAX_SENSE ?
        objective_value :
        -objective_value
    min_psd_eigenvalue = psd_count == 0 ? NaN : min_psd_eigenvalue
    max_cone_violation = max(
        max_scalar_cone_violation,
        psd_count == 0 ? 0.0 : max(0.0, -min_psd_eigenvalue),
    )
    max_equality_residual_relative = max_equality_residual / scale
    max_cone_violation_relative = max_cone_violation / scale
    improving_objective_relative = improving_objective / scale
    equality_pass = max_equality_residual_relative <= normalized_tolerance
    cone_pass = max_cone_violation_relative <= normalized_tolerance
    objective_pass = improving_objective_relative > normalized_tolerance
    verified = equality_pass && cone_pass && objective_pass
    rejection_reasons = String[]
    equality_pass || push!(rejection_reasons, "equality_residual")
    cone_pass || push!(rejection_reasons, "cone_violation")
    objective_pass || push!(rejection_reasons, "objective_direction")
    return (
        verified=verified,
        verdict=verified ? "accepted_floating_point_ray" : "rejected",
        rigor="floating_point_replay",
        equality_pass=equality_pass,
        cone_pass=cone_pass,
        objective_pass=objective_pass,
        rejection_reasons=Tuple(rejection_reasons),
        variable_count=length(variables),
        equality_count=equality_count,
        psd_count=psd_count,
        scalar_cone_count=scalar_cone_count,
        variable_scale=scale,
        absolute_tolerance=Float64(absolute_tolerance),
        relative_tolerance=Float64(relative_tolerance),
        normalized_tolerance=normalized_tolerance,
        max_equality_residual=max_equality_residual,
        max_equality_residual_relative=max_equality_residual_relative,
        min_psd_eigenvalue=min_psd_eigenvalue,
        min_psd_eigenvalue_relative=min_psd_eigenvalue / scale,
        max_scalar_cone_violation=max_scalar_cone_violation,
        max_cone_violation=max_cone_violation,
        max_cone_violation_relative=max_cone_violation_relative,
        objective_sense=string(objective_sense),
        objective_ray=objective_value,
        improving_objective=improving_objective,
        improving_objective_relative=improving_objective_relative,
    )
end

function main(args=ARGS)
    length(args) == 2 ||
        error("usage: verify_gap_ray.jl MODEL.mof.json[.gz] RAY.variables.tsv")
    result = verify(args[1], args[2])
    for field in propertynames(result)
        println(field, '\t', getproperty(result, field))
    end
    return result.verified ? 0 : 1
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(GapRayVerifier.main())
end
