module KagomeSourceAudit

using SHA

include(joinpath(@__DIR__, "GapRayPostprocess.jl"))
using .GapRayPostprocess

const MOI = GapRayPostprocess.MOI
const MOIU = MOI.Utilities
const BigRational = GapRayPostprocess.BigRational

export audit_kagome_source_assembly, pauli_real_action

const SOURCE_COMMIT = "a1171c906ff2cc2901e58c2426397a2f68c32bb7"
const PATCH_SHA256 =
    "332c0931ac810289aa3713af0948f259c01189270706af58b262d60d994d4abd"
const SOURCE_SHA256 = Dict(
    "SpectralGap.jl" =>
        "940cd72b9c4bea39b6daaefd2b9797c54df450cb0afbfb4d9322b2df1b3838bb",
    "basicfunction.jl" =>
        "2095cf7401355f37e9d17915b3ab29d44712d8e40f750eb8449f8c294229b03a",
    "sdp.jl" =>
        "4ea362723bd7601e67db3bc27f21a4a11506791ce2b9b82cb7e4a60b0f6bae10",
    "strengthening.jl" =>
        "de56b12b17049f81f689d4caef193b9dfd3bf50061fc78b1bf1547a748f7c57b",
)

file_sha256(path::AbstractString) =
    open(path, "r") do io
        bytes2hex(sha256(io))
    end

function source_gate(repo_root::AbstractString, spectralgap_root::AbstractString)
    patch = joinpath(
        repo_root,
        "tracks",
        "polyopt",
        "solutions",
        "sdp-gap-seekers",
        "spectralgap_a1171c9.patch",
    )
    patch_hash = file_sha256(patch)
    patch_hash == PATCH_SHA256 || error("SpectralGap patch SHA-256 mismatch")
    commit = readchomp(`git -C $spectralgap_root rev-parse HEAD`)
    commit == SOURCE_COMMIT || error("SpectralGap base commit mismatch")
    hashes = Dict(
        filename => file_sha256(joinpath(spectralgap_root, "src", filename))
        for filename in keys(SOURCE_SHA256)
    )
    hashes == SOURCE_SHA256 ||
        error("patched SpectralGap source SHA-256 mismatch")
    return (commit=commit, patch_sha256=patch_hash, files=hashes)
end

function exact_coefficient(value)
    isreal(value) ||
        error("source realification left a complex coefficient: $value")
    return rationalize(BigInt, Float64(real(value)); tol=1e-14)
end

function add_coefficient!(
    rows::Vector{Dict{Int,BigRational}},
    row_index::Union{Nothing,Integer},
    ordinal::Integer,
    coefficient,
)
    isnothing(row_index) && error("source support was not present in tsupp")
    rational = exact_coefficient(coefficient)
    row = rows[row_index]
    row[ordinal] = get(row, ordinal, zero(BigRational)) + rational
    iszero(row[ordinal]) && delete!(row, ordinal)
    return nothing
end

canonical_row(row::Dict{Int,BigRational}) =
    sort!(collect(row); by=first)

canonical_row(row::ExactAffineRow) = row.terms

function row_mismatches(expected, actual)
    length(expected) == length(actual) ||
        return collect(1:max(length(expected), length(actual)))
    return [
        index for index in eachindex(expected) if
        canonical_row(expected[index]) != canonical_row(actual[index])
    ]
end

function row_digest(rows)
    lines = String["kagome-source-affine-rows-v1"]
    for (row_index, row) in enumerate(rows)
        terms = join(
            (
                string(
                    ordinal,
                    ":",
                    numerator(coefficient),
                    "/",
                    denominator(coefficient),
                ) for (ordinal, coefficient) in canonical_row(row)
            ),
            ",",
        )
        push!(lines, string(row_index, "|", terms))
    end
    return bytes2hex(sha256(codeunits(join(lines, "\n") * "\n")))
end

function block_ordinal(block::PSDDirectionBlock, row::Int, column::Int)
    coordinate = block.coordinates[MOIU.trimap(row, column)]
    length(coordinate.terms) == 1 ||
        error("PSD coordinate is not a single source variable")
    ordinal, coefficient = only(coordinate.terms)
    coefficient == 1 || error("PSD coordinate coefficient is not one")
    return ordinal
end

function decode_pauli_string(code::Integer, site_count::Integer)
    code >= 0 || error("Pauli code must be nonnegative")
    labels = zeros(Int, site_count)
    remaining = code
    for site in site_count:-1:1
        labels[site] = remaining % 4
        remaining ÷= 4
    end
    iszero(remaining) || error("Pauli code exceeds the requested site count")
    return labels
end

"""
    pauli_real_action(labels, column_state)

Return the zero-based row state and exact real matrix entry for a tensor
product of `I/X/Y/Z`, encoded by `0/1/2/3`, acting on a zero-based
computational-basis column. This helper accepts only an even number of `Y`
factors, for which every nonzero matrix entry is real.
"""
function pauli_real_action(
    labels::AbstractVector{<:Integer},
    column_state::Integer,
)
    site_count = length(labels)
    0 <= column_state < (1 << site_count) ||
        error("computational-basis column is out of range")
    y_count = count(==(2), labels)
    iseven(y_count) || error("Pauli action is not real")
    row_state = Int(column_state)
    sign_exponent = y_count ÷ 2
    for (site, label) in enumerate(labels)
        label in 0:3 || error("invalid Pauli label")
        bit_position = site_count - site
        bit = (column_state >> bit_position) & 1
        if label in (1, 2)
            row_state ⊻= 1 << bit_position
        end
        if label in (2, 3)
            sign_exponent += bit
        end
    end
    return row_state, iseven(sign_exponent) ? 1 : -1
end

function strengthening_blocks(site_count::Integer)
    site_count == 9 ||
        error("the audited Kagome strengthening block uses nine sites")
    return [
        [
            state + 1 for state in 0:((1 << site_count) - 1)
            if count_ones(state) == weight
        ]
        for weight in 0:4
    ]
end

function kagome_source_objects(spectralgap)
    N = 13
    d = 3
    lso = 5
    gamma = Float64(159 // 125)
    triples = [
        [1, 2, 3],
        [1, 4, 5],
        [2, 6, 7],
        [3, 8, 9],
        [4, 10, 11],
        [5, 12, 13],
    ]
    edges = Vector{Vector{Int}}()
    inner_triples = [[1, 2, 3], [1, 4, 5]]
    inner_edges = Vector{Vector{Int}}()
    supports = Vector{Int}[]
    for triangle in triples
        for (left, right) in (
            (triangle[1], triangle[2]),
            (triangle[1], triangle[3]),
            (triangle[2], triangle[3]),
        )
            for axis in 0:2
                push!(supports, [3left - 2 + axis, 3right - 2 + axis])
            end
        end
    end
    coefficients = fill(0.25, length(supports))
    H = spectralgap.ncpoly(supports, coefficients)
    basis = [
        spectralgap.get_kagome_basis(
            N,
            triples,
            edges,
            d;
            label=label,
        )
        for label in (1, 2)
    ]
    gbasis = [
        spectralgap.get_kagome_bulkbasis(
            N,
            inner_triples,
            inner_edges,
            d - 1;
            label=label,
        )
        for label in (1, 2)
    ]
    lb = length.(basis)
    lgb = length.(gbasis)

    tsupp = Vector{Vector{Int}}[]
    for block in eachindex(basis), j in 1:lb[block], k in j:lb[block]
        word, coefficient = spectralgap.reduce!(
            [basis[block][j][1]; basis[block][k][1]],
            N;
            model="kagome",
        )
        iszero(coefficient) && continue
        state_factors = [basis[block][j][2]; basis[block][k][2]]
        push!(
            tsupp,
            isempty(word) ? sort(state_factors) :
            sort([state_factors; [word]]),
        )
    end
    for block in eachindex(gbasis), i in 1:lgb[block], j in i:lgb[block]
        words = spectralgap.PSDstate_entry(
            gbasis[block][i][1],
            gbasis[block][j][1],
            H,
            N;
            model="kagome",
        )[1]
        state_factors = [gbasis[block][i][2]; gbasis[block][j][2]]
        for word in words
            push!(
                tsupp,
                isempty(word) ? sort(state_factors) :
                sort([state_factors; [word]]),
            )
        end
        if !spectralgap.isz(gbasis[block][i][1]; model="kagome") &&
           !spectralgap.isz(gbasis[block][j][1]; model="kagome")
            rotated = [
                spectralgap.reduce_perm(gbasis[block][i][1]),
                spectralgap.reduce_perm(gbasis[block][j][1]),
            ]
            push!(tsupp, sort([state_factors; rotated]))
        end
    end

    strengthening_strings = Vector{Tuple{Vector{Int},Vector{Vector{Int}}}}()
    for code in 0:(4^9 - 1)
        labels = decode_pauli_string(code, 9)
        all(axis -> iseven(count(==(axis), labels)), 1:3) || continue
        sites = findall(!iszero, labels)
        word = 3 .* (sites .- 1) .+ labels[sites]
        key = [spectralgap.reduce_perm(word)]
        push!(tsupp, key)
        push!(strengthening_strings, (labels, key))
    end
    sort!(tsupp)
    unique!(tsupp)

    monomials = spectralgap.generate_mons_kagome(5, lso)
    monomials = spectralgap.filter_mons(
        monomials,
        tsupp,
        H,
        N;
        model="kagome",
    )
    return (
        N=N,
        d=d,
        lso=lso,
        gamma=gamma,
        triples=triples,
        edges=edges,
        inner_triples=inner_triples,
        inner_edges=inner_edges,
        supports=supports,
        coefficients=coefficients,
        H=H,
        basis=basis,
        gbasis=gbasis,
        lb=lb,
        lgb=lgb,
        tsupp=tsupp,
        monomials=monomials,
        strengthening_strings=strengthening_strings,
        strengthening_blocks=strengthening_blocks(9),
    )
end

function add_basis_rows!(rows, source, problem, spectralgap)
    for block in eachindex(source.basis)
        for j in 1:source.lb[block], k in j:source.lb[block]
            word, coefficient = spectralgap.reduce!(
                [
                    source.basis[block][j][1];
                    source.basis[block][k][1]
                ],
                source.N;
                realify=true,
                model="kagome",
            )
            iszero(coefficient) && continue
            state_factors = [
                source.basis[block][j][2];
                source.basis[block][k][2]
            ]
            key = isempty(word) ? sort(state_factors) :
                  sort([state_factors; [word]])
            add_coefficient!(
                rows,
                spectralgap.bfind(source.tsupp, key),
                block_ordinal(problem.psd_blocks[block], j, k),
                (j == k ? 1 : 2) * coefficient,
            )
        end
    end
    return nothing
end

function add_gap_rows!(rows, source, problem, spectralgap)
    gap_offset = length(source.basis)
    for block in eachindex(source.gbasis)
        for i in 1:source.lgb[block], j in i:source.lgb[block]
            words, coefficients = spectralgap.PSDstate_entry(
                source.gbasis[block][i][1],
                source.gbasis[block][j][1],
                source.H,
                source.N;
                realify=true,
                model="kagome",
            )
            state_factors = [
                source.gbasis[block][i][2];
                source.gbasis[block][j][2]
            ]
            ordinal = block_ordinal(
                problem.psd_blocks[gap_offset + block],
                i,
                j,
            )
            factor = i == j ? 1 : 2
            for (word, coefficient) in zip(words, coefficients)
                key = isempty(word) ? sort(state_factors) :
                      sort([state_factors; [word]])
                add_coefficient!(
                    rows,
                    spectralgap.bfind(source.tsupp, key),
                    ordinal,
                    factor * coefficient,
                )
            end

            word, coefficient = spectralgap.reduce!(
                [
                    source.gbasis[block][i][1];
                    source.gbasis[block][j][1]
                ],
                source.N;
                realify=true,
                model="kagome",
            )
            if !iszero(coefficient)
                key = isempty(word) ? sort(state_factors) :
                      sort([state_factors; [word]])
                add_coefficient!(
                    rows,
                    spectralgap.bfind(source.tsupp, key),
                    ordinal,
                    -factor * coefficient * source.gamma,
                )
            end

            if !spectralgap.isz(
                source.gbasis[block][i][1];
                model="kagome",
            ) && !spectralgap.isz(
                source.gbasis[block][j][1];
                model="kagome",
            )
                rotated = [
                    spectralgap.reduce_perm(source.gbasis[block][i][1]),
                    spectralgap.reduce_perm(source.gbasis[block][j][1]),
                ]
                add_coefficient!(
                    rows,
                    spectralgap.bfind(
                        source.tsupp,
                        sort([state_factors; rotated]),
                    ),
                    ordinal,
                    factor * source.gamma,
                )
            end
        end
    end
    return nothing
end

function add_stationarity_rows!(
    rows,
    source,
    scalar_ordinals,
    spectralgap,
)
    length(scalar_ordinals) == length(source.monomials) + 1 ||
        error("MOF scalar-variable count does not match stationarity plus lambda")
    for (monomial_index, monomial) in enumerate(source.monomials)
        ordinal = scalar_ordinals[monomial_index]
        for (term_index, term) in enumerate(source.H.supp)
            word, coefficient = spectralgap.reduce!(
                [term; monomial],
                source.N;
                model="kagome",
            )
            iszero(imag(coefficient)) && continue
            add_coefficient!(
                rows,
                spectralgap.bfind(source.tsupp, [word]),
                ordinal,
                source.H.coe[term_index] * imag(coefficient),
            )
        end
    end
    return nothing
end

function add_strengthening_rows!(rows, source, problem, spectralgap)
    block_offset = length(source.basis) + length(source.gbasis)
    block_lookup = zeros(Int, 1 << 9)
    local_lookup = zeros(Int, 1 << 9)
    for (block_index, states) in enumerate(source.strengthening_blocks)
        for (local_index, global_index) in enumerate(states)
            block_lookup[global_index] = block_index
            local_lookup[global_index] = local_index
        end
    end
    for (string_index, (labels, key)) in
        enumerate(source.strengthening_strings)
        row_index = spectralgap.bfind(source.tsupp, key)
        for column_state in 0:((1 << 9) - 1)
            row_state, coefficient =
                pauli_real_action(labels, column_state)
            column_global = column_state + 1
            row_global = row_state + 1
            block_index = block_lookup[column_global]
            iszero(block_index) && continue
            block_lookup[row_global] == block_index || continue
            local_column = local_lookup[column_global]
            local_row = local_lookup[row_global]
            ordinal = block_ordinal(
                problem.psd_blocks[block_offset + block_index],
                min(local_row, local_column),
                max(local_row, local_column),
            )
            add_coefficient!(
                rows,
                row_index,
                ordinal,
                coefficient,
            )
        end
        if string_index % 5000 == 0
            println(
                "progress\tstrengthening_strings\t",
                string_index,
                '/',
                length(source.strengthening_strings),
            )
            flush(stdout)
        end
    end
    return nothing
end

function expected_rows(source, problem::ExactRayProblem, spectralgap)
    dimensions = [
        source.lb;
        source.lgb;
        length.(source.strengthening_blocks)
    ]
    [block.dimension for block in problem.psd_blocks] == dimensions ||
        error("MOF PSD blocks do not match source basis dimensions")
    length(problem.equalities) == length(source.tsupp) ||
        error("MOF equality count does not match source support count")
    rows = [Dict{Int,BigRational}() for _ in source.tsupp]
    add_basis_rows!(rows, source, problem, spectralgap)
    add_gap_rows!(rows, source, problem, spectralgap)

    psd_ordinals = Set(
        ordinal
        for block in problem.psd_blocks
        for coordinate in block.coordinates
        for (ordinal, _) in coordinate.terms
    )
    scalar_ordinals = sort!(setdiff(
        collect(1:problem.variable_count),
        collect(psd_ordinals),
    ))
    add_stationarity_rows!(rows, source, scalar_ordinals, spectralgap)
    add_strengthening_rows!(rows, source, problem, spectralgap)
    lambda_ordinal = last(scalar_ordinals)
    add_coefficient!(rows, 1, lambda_ordinal, 1)
    return rows, scalar_ordinals, lambda_ordinal
end

"""
    audit_kagome_source_assembly(
        spectralgap, repo_root, spectralgap_root, model_path,
    )

Recompute every affine coefficient of the locked Kagome N=13, d=3, lso=5,
gamma=159/125 sign-symmetric source formulation and compare it exactly with
the intended rational reconstruction of the exported MOF. The audit includes
all five `posepsd9!` strengthening blocks. No optimizer is constructed or
invoked.
"""
function audit_kagome_source_assembly(
    spectralgap,
    repo_root::AbstractString,
    spectralgap_root::AbstractString,
    model_path::AbstractString,
)
    gate = source_gate(repo_root, spectralgap_root)
    problem = extract_exact_problem(
        model_path;
        coefficient_tolerance=1e-14,
    )
    source = kagome_source_objects(spectralgap)
    expected, scalar_ordinals, lambda_ordinal =
        expected_rows(source, problem, spectralgap)
    mismatches = row_mismatches(expected, problem.equalities)
    all(iszero, problem.equality_offsets) ||
        error("exported MOF has a nonzero affine right-hand side")
    problem.objective_sense == MOI.MAX_SENSE ||
        error("exported MOF objective is not maximize")
    problem.objective.terms == [lambda_ordinal => one(BigRational)] ||
        error("exported MOF objective is not +lambda")
    isempty(mismatches) ||
        error(
            "source/MOF affine coefficient mismatch in " *
            "$(length(mismatches)) rows; first=$(first(mismatches))",
        )
    expected_coefficients = sort(unique(
        coefficient for row in expected for (_, coefficient) in row
    ))
    actual_coefficients = sort(unique(
        coefficient
        for row in problem.equalities
        for (_, coefficient) in row.terms
    ))
    expected_coefficients == actual_coefficients ||
        error("source/MOF coefficient inventories differ")
    return (
        source_gate=gate,
        hamiltonian_supports=source.supports,
        hamiltonian_coefficients=source.coefficients,
        pauli_encoding="3i-2=X_i;3i-1=Y_i;3i=Z_i",
        patch=(
            sites=source.N,
            triangles=source.triples,
            edges=source.edges,
            inner_triangles=source.inner_triples,
            inner_edges=source.inner_edges,
        ),
        state_class=(
            sign_symmetry="odd parity of every X/Y/Z axis removed by kagome isz",
            spin_symmetry="cyclic X-to-Y-to-Z representatives via reduce_perm",
            semantics="symmetry-restricted infinite-volume KMS ground-state state-polynomial relaxation",
        ),
        matrix_orientation=(
            source_entry="upper triangle",
            mof_storage="MOI PositiveSemidefiniteConeTriangle trimap",
            off_diagonal_factor="symmetric source entries accumulated twice",
        ),
        gamma=159 // 125,
        basis_dimensions=source.lb,
        gap_basis_dimensions=source.lgb,
        strengthening_dimensions=length.(source.strengthening_blocks),
        stationarity_monomials=length(source.monomials),
        scalar_variable_ordinals=scalar_ordinals,
        lambda_ordinal=lambda_ordinal,
        variable_count=problem.variable_count,
        equality_count=length(problem.equalities),
        equality_offsets_zero=true,
        psd_dimensions=[block.dimension for block in problem.psd_blocks],
        objective_sense=problem.objective_sense,
        objective=problem.objective.terms,
        coefficient_inventory=actual_coefficients,
        affine_row_sha256=row_digest(expected),
        rows_compared=length(expected),
        coefficient_mismatches=length(mismatches),
        optimizer_invoked=false,
    )
end

end
