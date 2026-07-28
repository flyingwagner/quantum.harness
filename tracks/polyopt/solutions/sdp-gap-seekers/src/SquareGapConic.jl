module SquareGapConic

import JuMP
using ..SquareJ1J2Prototype: PauliWord
import ..GenericGapModel
using ..GenericGapModel:
    GapProblem,
    NoStateSymmetry,
    StateMonomial
import ..CoreMGK
using ..CoreMGK:
    CoreMGKPlan,
    GaussianRational,
    ScalarMoment,
    a_gamma_coefficients,
    core_mgk_pair,
    core_mgk_plan
import ..SharedCoreWire
using ..SharedCoreWire:
    canonical_scalar_row,
    encode_value,
    row_id

const MOI = JuMP.MOI
const MOIU = MOI.Utilities
const BigRational = Rational{BigInt}

export ExactAffineConstraint,
    RealPSDPlan,
    RealPSDTerm,
    SquareConicPlan,
    audit_rendered_mof,
    build_square_conic_plan,
    realify_hermitian_coefficient,
    render_mof

struct ExactAffineConstraint
    kind::Symbol
    source_id::String
    terms::Vector{Pair{Int,BigRational}}
    rhs::BigRational
end

struct RealPSDTerm
    output_index::Int
    row_index::Int
    coefficient::BigRational
end

struct RealPSDPlan
    role::Symbol
    complex_dimension::Int
    real_dimension::Int
    terms::Vector{RealPSDTerm}
end

struct SquareConicPlan
    source::CoreMGKPlan
    gamma::BigRational
    rows::Vector{ScalarMoment}
    row_ids::Vector{String}
    normalization::ExactAffineConstraint
    stationarity::Vector{ExactAffineConstraint}
    stationarity_selector_entries::Int
    stationarity_exact_duplicates_removed::Int
    psd_blocks::Vector{RealPSDPlan}
    objective_sense::Symbol
    objective_terms::Vector{Pair{Int,BigRational}}
end

struct RawPSDTerm
    output_index::Int
    row::ScalarMoment
    coefficient::BigRational
end

function exact_rational(value::Integer)
    return BigInt(value) // BigInt(1)
end

function exact_rational(value::Rational)
    return BigInt(numerator(value)) // BigInt(denominator(value))
end

exact_rational(value::AbstractFloat) =
    throw(ArgumentError("Square conic rendering requires exact rational gamma"))

function byte_vector_less(left, right)
    for index in 1:min(length(left), length(right))
        left[index] == right[index] || return left[index] < right[index]
    end
    return length(left) < length(right)
end

function canonical_row_order(left::ScalarMoment, right::ScalarMoment)
    left_bytes = encode_value(canonical_scalar_row(left))
    right_bytes = encode_value(canonical_scalar_row(right))
    left_bytes == right_bytes && return row_id(left) < row_id(right)
    return byte_vector_less(left_bytes, right_bytes)
end

"""
Map one upper-triangle complex Hermitian coefficient to the upper triangle of
`[Re(A) -Im(A); Im(A) Re(A)]`.

Returned triples are `(row, column, real_coefficient)` with `row <= column`.
No trace-inner-product factor of two is applied here: these are primal matrix
entries, not dual `Tr(AQ)` packed coefficients.
"""
function realify_hermitian_coefficient(
    coefficient::GaussianRational,
    j::Int,
    k::Int,
    dimension::Int,
)
    1 <= j <= k <= dimension ||
        throw(ArgumentError("complex matrix pair must satisfy 1 <= j <= k <= n"))
    real_part = real(coefficient)
    imag_part = imag(coefficient)
    result = Tuple{Int,Int,BigRational}[]
    if j == k
        iszero(imag_part) ||
            throw(ArgumentError("Hermitian diagonal coefficient is not real"))
        iszero(real_part) || begin
            push!(result, (j, j, real_part))
            push!(result, (j + dimension, j + dimension, real_part))
        end
        return result
    end
    if !iszero(real_part)
        push!(result, (j, k, real_part))
        push!(
            result,
            (j + dimension, k + dimension, real_part),
        )
    end
    if !iszero(imag_part)
        push!(result, (j, k + dimension, -imag_part))
        push!(result, (k, j + dimension, imag_part))
    end
    return result
end

function append_realified!(
    terms::Vector{RawPSDTerm},
    row::ScalarMoment,
    coefficient::GaussianRational,
    j::Int,
    k::Int,
    dimension::Int,
)
    for (matrix_row, matrix_column, real_coefficient) in
        realify_hermitian_coefficient(coefficient, j, k, dimension)
        push!(
            terms,
            RawPSDTerm(
                MOIU.trimap(matrix_row, matrix_column),
                row,
                real_coefficient,
            ),
        )
    end
    return terms
end

function positive_raw_terms(source::CoreMGKPlan)
    dimension = length(source.positive_basis.entries)
    terms = RawPSDTerm[]
    for j in 1:dimension, k in j:dimension
        wiring = core_mgk_pair(source, :positive, j, k)
        component = only(wiring.component_records)
        component.component == :M || error("positive core component is not M")
        for record in component.coefficients
            append_realified!(
                terms,
                record.row,
                record.coefficient,
                j,
                k,
                dimension,
            )
        end
    end
    return terms
end

function gap_raw_terms(source::CoreMGKPlan, gamma::BigRational)
    dimension = length(source.gap_basis.entries)
    terms = RawPSDTerm[]
    for j in 1:dimension, k in j:dimension
        wiring = core_mgk_pair(source, :gap, j, k)
        for (row, affine) in a_gamma_coefficients(wiring)
            coefficient = affine.constant + gamma * affine.gamma
            append_realified!(terms, row, coefficient, j, k, dimension)
        end
    end
    return terms
end

function add_gaussian!(
    row::Dict{ScalarMoment,GaussianRational},
    moment::ScalarMoment,
    coefficient::GaussianRational,
)
    iszero(coefficient) && return row
    value = get(row, moment, zero(GaussianRational)) + coefficient
    iszero(value) ? delete!(row, moment) : (row[moment] = value)
    return row
end

function stationarity_raw_rows(
    problem::GapProblem,
    source::CoreMGKPlan,
)
    problem.symmetry isa NoStateSymmetry ||
        throw(ArgumentError("stationarity renderer supports no unapplied symmetry"))
    maximum_degree = 2problem.d - 2
    site_ids = sort!(copy(problem.patch.inner_ids))
    entries = GenericGapModel.one_symbol_entries(site_ids, maximum_degree)
    candidates = Tuple{Symbol,String,Dict{ScalarMoment,BigRational}}[]
    for (entry_index, entry) in enumerate(entries)
        gaussian = Dict{ScalarMoment,GaussianRational}()
        commutator = CoreMGK.commutator_polynomial(
            source.hamiltonian_terms,
            entry.operator_word,
        )
        for (word, coefficient) in commutator
            moment = CoreMGK.scalarize(entry.state_symbols, word)
            add_gaussian!(gaussian, moment, coefficient)
        end
        real_row = Dict{ScalarMoment,BigRational}(
            moment => real(coefficient)
            for (moment, coefficient) in gaussian
            if !iszero(real(coefficient))
        )
        imag_row = Dict{ScalarMoment,BigRational}(
            moment => imag(coefficient)
            for (moment, coefficient) in gaussian
            if !iszero(imag(coefficient))
        )
        isempty(real_row) || push!(
            candidates,
            (:stationarity_real, "q$(entry_index):real", real_row),
        )
        isempty(imag_row) || push!(
            candidates,
            (:stationarity_imag, "q$(entry_index):imag", imag_row),
        )
    end

    # Remove only exact duplicates. The first source entry/component in the
    # deterministic stationarity selector owns the retained row.
    seen = Set{Any}()
    result = Tuple{Symbol,String,Dict{ScalarMoment,BigRational}}[]
    for candidate in candidates
        coefficients = last(candidate)
        signature = Tuple(sort!(
            collect(coefficients);
            by=pair -> row_id(first(pair)),
        ))
        signature in seen && continue
        push!(seen, signature)
        push!(result, candidate)
    end
    return result, length(entries), length(candidates) - length(result)
end

function canonicalize_psd_terms(
    raw_terms::Vector{RawPSDTerm},
    row_indices::Dict{ScalarMoment,Int},
)
    accumulated = Dict{Tuple{Int,Int},BigRational}()
    for term in raw_terms
        key = (term.output_index, row_indices[term.row])
        value = get(accumulated, key, zero(BigRational)) + term.coefficient
        iszero(value) ? delete!(accumulated, key) : (accumulated[key] = value)
    end
    terms = [
        RealPSDTerm(output_index, row_index, coefficient)
        for ((output_index, row_index), coefficient) in accumulated
    ]
    sort!(terms; by=term -> (term.output_index, term.row_index))
    return terms
end

function affine_constraint(
    kind::Symbol,
    source_id::AbstractString,
    raw_terms::Dict{ScalarMoment,BigRational},
    rhs::BigRational,
    row_indices::Dict{ScalarMoment,Int},
)
    terms = sort!(
        [
            row_indices[row] => coefficient
            for (row, coefficient) in raw_terms
            if !iszero(coefficient)
        ];
        by=first,
    )
    return ExactAffineConstraint(kind, String(source_id), terms, rhs)
end

"""
Build a solver-free exact primal feasibility plan for one structured Square
problem. Numeric gamma is evaluated only here; the underlying core M/K/G tensor
remains gamma-independent.
"""
function build_square_conic_plan(problem::GapProblem)
    gamma = exact_rational(problem.gamma)
    source = core_mgk_plan(problem)
    positive_terms = positive_raw_terms(source)
    gap_terms = gap_raw_terms(source, gamma)
    stationarity_raw, stationarity_entry_count, duplicate_count =
        stationarity_raw_rows(problem, source)
    identity = ScalarMoment(PauliWord[])

    row_set = Set{ScalarMoment}([identity])
    for term in Iterators.flatten((positive_terms, gap_terms))
        push!(row_set, term.row)
    end
    for (_, _, coefficients) in stationarity_raw
        union!(row_set, keys(coefficients))
    end
    rows = collect(row_set)
    sort!(rows; lt=canonical_row_order)
    row_ids = row_id.(rows)
    length(unique(row_ids)) == length(row_ids) ||
        error("canonical scalar-row content ID collision")
    row_indices = Dict(row => index for (index, row) in enumerate(rows))

    normalization = affine_constraint(
        :normalization,
        "L(1)=1",
        Dict(identity => one(BigRational)),
        one(BigRational),
        row_indices,
    )
    stationarity = [
        affine_constraint(kind, source_id, coefficients, zero(BigRational), row_indices)
        for (kind, source_id, coefficients) in stationarity_raw
    ]
    positive_dimension = length(source.positive_basis.entries)
    gap_dimension = length(source.gap_basis.entries)
    psd_blocks = [
        RealPSDPlan(
            :positive,
            positive_dimension,
            2positive_dimension,
            canonicalize_psd_terms(positive_terms, row_indices),
        ),
        RealPSDPlan(
            :gap,
            gap_dimension,
            2gap_dimension,
            canonicalize_psd_terms(gap_terms, row_indices),
        ),
    ]

    # Store selector diagnostics in the source plan's state-class string only
    # as display is intentionally avoided; enforce them as construction gates.
    stationarity_entry_count > 0 || error("stationarity selector is empty")
    duplicate_count >= 0 || error("invalid stationarity duplicate count")
    return SquareConicPlan(
        source,
        gamma,
        rows,
        row_ids,
        normalization,
        stationarity,
        stationarity_entry_count,
        duplicate_count,
        psd_blocks,
        :feasibility,
        Pair{Int,BigRational}[],
    )
end

function checked_float(value::BigRational)
    rendered = Float64(value)
    isfinite(rendered) || error("exact conic coefficient overflows Float64")
    rationalize(BigInt, rendered; tol=1e-14) == value ||
        error("exact conic coefficient is not recoverable at tolerance 1e-14")
    return rendered
end

function scalar_function(
    variables,
    terms::Vector{Pair{Int,BigRational}},
)
    return MOI.ScalarAffineFunction(
        [
            MOI.ScalarAffineTerm(
                checked_float(coefficient),
                variables[row_index],
            )
            for (row_index, coefficient) in terms
        ],
        0.0,
    )
end

"""Render the exact plan to MOF without constructing or invoking an optimizer."""
function render_mof(plan::SquareConicPlan, output_path::AbstractString)
    model = MOIU.Model{Float64}()
    variables = MOI.add_variables(model, length(plan.rows))
    for (index, variable) in enumerate(variables)
        MOI.set(model, MOI.VariableName(), variable, "y:$index:$(plan.row_ids[index])")
    end
    for constraint in [plan.normalization; plan.stationarity]
        MOI.add_constraint(
            model,
            scalar_function(variables, constraint.terms),
            MOI.EqualTo(checked_float(constraint.rhs)),
        )
    end
    for block in plan.psd_blocks
        output_dimension = div(
            block.real_dimension * (block.real_dimension + 1),
            2,
        )
        terms = [
            MOI.VectorAffineTerm(
                term.output_index,
                MOI.ScalarAffineTerm(
                    checked_float(term.coefficient),
                    variables[term.row_index],
                ),
            )
            for term in block.terms
        ]
        function_value = MOI.VectorAffineFunction(
            terms,
            zeros(Float64, output_dimension),
        )
        MOI.add_constraint(
            model,
            function_value,
            MOI.PositiveSemidefiniteConeTriangle(block.real_dimension),
        )
    end
    MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    MOI.write_to_file(model, String(output_path))
    return (
        variables=length(variables),
        affine_equalities=1 + length(plan.stationarity),
        psd_dimensions=[block.real_dimension for block in plan.psd_blocks],
        objective_sense=:feasibility,
        optimizer_invoked=false,
    )
end

function accumulated_scalar_terms(function_value, variable_positions)
    result = Dict{Int,Float64}()
    for term in function_value.terms
        position = variable_positions[term.variable]
        value = get(result, position, 0.0) + term.coefficient
        iszero(value) ? delete!(result, position) : (result[position] = value)
    end
    return result
end

function accumulated_vector_terms(function_value, variable_positions)
    result = Dict{Tuple{Int,Int},Float64}()
    for term in function_value.terms
        position = variable_positions[term.scalar_term.variable]
        key = (term.output_index, position)
        value = get(result, key, 0.0) + term.scalar_term.coefficient
        iszero(value) ? delete!(result, key) : (result[key] = value)
    end
    return result
end

"""
Independently parse a rendered MOF and compare every affine and PSD coefficient
to the exact source plan after the declared checked Float64 conversion.
No optimizer is constructed.
"""
function audit_rendered_mof(plan::SquareConicPlan, model_path::AbstractString)
    model = MOI.FileFormats.Model(filename=String(model_path))
    MOI.read_from_file(model, String(model_path))
    variables = sort!(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    length(variables) == length(plan.rows) ||
        error("MOF variable count differs from the exact conic plan")
    variable_positions = Dict(
        variable => position for (position, variable) in enumerate(variables)
    )
    for (position, variable) in enumerate(variables)
        expected_name = "y:$position:$(plan.row_ids[position])"
        MOI.get(model, MOI.VariableName(), variable) == expected_name ||
            error("MOF variable name/order differs from the exact conic plan")
    end

    equality_indices = sort!(
        MOI.get(
            model,
            MOI.ListOfConstraintIndices{
                MOI.ScalarAffineFunction{Float64},
                MOI.EqualTo{Float64},
            }(),
        );
        by=constraint -> constraint.value,
    )
    exact_equalities = [plan.normalization; plan.stationarity]
    length(equality_indices) == length(exact_equalities) ||
        error("MOF affine-equality count differs from the exact conic plan")
    for (constraint_index, exact_constraint) in
        zip(equality_indices, exact_equalities)
        function_value = MOI.get(
            model,
            MOI.ConstraintFunction(),
            constraint_index,
        )
        set_value = MOI.get(model, MOI.ConstraintSet(), constraint_index)
        iszero(function_value.constant) ||
            error("MOF affine equality has a nonzero function constant")
        set_value.value == checked_float(exact_constraint.rhs) ||
            error("MOF affine right-hand side differs from the exact plan")
        expected_terms = Dict(
            row_index => checked_float(coefficient)
            for (row_index, coefficient) in exact_constraint.terms
        )
        accumulated_scalar_terms(function_value, variable_positions) ==
            expected_terms ||
            error("MOF affine coefficients differ from the exact plan")
    end

    psd_indices = sort!(
        MOI.get(
            model,
            MOI.ListOfConstraintIndices{
                MOI.VectorAffineFunction{Float64},
                MOI.PositiveSemidefiniteConeTriangle,
            }(),
        );
        by=constraint -> constraint.value,
    )
    length(psd_indices) == length(plan.psd_blocks) ||
        error("MOF PSD-block count differs from the exact conic plan")
    for (constraint_index, block) in zip(psd_indices, plan.psd_blocks)
        function_value = MOI.get(
            model,
            MOI.ConstraintFunction(),
            constraint_index,
        )
        set_value = MOI.get(model, MOI.ConstraintSet(), constraint_index)
        set_value.side_dimension == block.real_dimension ||
            error("MOF PSD dimension differs from the exact conic plan")
        all(iszero, function_value.constants) ||
            error("MOF PSD function contains a nonzero constant")
        expected_terms = Dict(
            (term.output_index, term.row_index) =>
                checked_float(term.coefficient)
            for term in block.terms
        )
        accumulated_vector_terms(function_value, variable_positions) ==
            expected_terms ||
            error("MOF PSD coefficients differ from the exact conic plan")
    end
    MOI.get(model, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE ||
        error("MOF objective sense is not feasibility")
    isempty(plan.objective_terms) ||
        error("exact conic plan unexpectedly has objective terms")
    return (
        variables=length(variables),
        affine_equalities=length(equality_indices),
        psd_dimensions=[
            MOI.get(model, MOI.ConstraintSet(), constraint).side_dimension
            for constraint in psd_indices
        ],
        affine_coefficients=sum(length, (
            constraint.terms for constraint in exact_equalities
        )),
        psd_coefficients=sum(length, (
            block.terms for block in plan.psd_blocks
        )),
        exact_coefficient_match=true,
        objective_sense=:feasibility,
        optimizer_invoked=false,
    )
end

end
