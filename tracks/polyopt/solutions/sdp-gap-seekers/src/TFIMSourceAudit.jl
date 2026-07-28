module TFIMSourceAudit

using SHA

include(joinpath(@__DIR__, "GapRayPostprocess.jl"))
using .GapRayPostprocess

const MOI = GapRayPostprocess.MOI
const MOIU = MOI.Utilities
const BigRational = GapRayPostprocess.BigRational

export audit_tfim_source_assembly, row_mismatches

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

function exact_coefficient(value)
    isreal(value) || error("source realification left a complex coefficient: $value")
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
    lines = String["tfim-source-affine-rows-v1"]
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
    hashes == SOURCE_SHA256 || error("patched SpectralGap source SHA-256 mismatch")
    return (commit=commit, patch_sha256=patch_hash, files=hashes)
end

function tfim_source_objects(spectralgap)
    N = 9
    d = 2
    lso = 6
    gamma = Float64(201 // 800)
    supports = [[3i, 3(i + 1)] for i in 1:(N - 1)]
    append!(supports, [[3i - 2] for i in 1:N])
    coefficients = [fill(-1.0, N - 1); fill(0.5, N)]
    H = spectralgap.ncpoly(supports, coefficients)
    basis = [spectralgap.get_basis(N, d; label=i) for i in (1, 2)]
    gbasis = [spectralgap.get_bulkbasis(N, d - 1; label=i) for i in (1, 2)]
    lb = length.(basis)
    lgb = length.(gbasis)

    tsupp = Vector{Vector{Int}}[]
    for block in eachindex(basis), j in 1:lb[block], k in j:lb[block]
        bi, coefficient = spectralgap.reduce!(
            [basis[block][j][1]; basis[block][k][1]],
            N,
        )
        iszero(coefficient) && continue
        state_factors = [basis[block][j][2]; basis[block][k][2]]
        push!(tsupp, isempty(bi) ? sort(state_factors) :
              sort([state_factors; [bi]]))
    end
    for block in eachindex(gbasis), i in 1:lgb[block], j in i:lgb[block]
        words = spectralgap.PSDstate_entry(
            gbasis[block][i][1],
            gbasis[block][j][1],
            H,
            N,
        )[1]
        state_factors = [gbasis[block][i][2]; gbasis[block][j][2]]
        for word in words
            push!(
                tsupp,
                isempty(word) ? sort(state_factors) :
                sort([state_factors; [word]]),
            )
        end
        if !spectralgap.isz(gbasis[block][i][1]; model="Ising") &&
           !spectralgap.isz(gbasis[block][j][1]; model="Ising")
            mirrored = [
                spectralgap.reduce_mirror(gbasis[block][i][1], N),
                spectralgap.reduce_mirror(gbasis[block][j][1], N),
            ]
            push!(tsupp, sort([state_factors; mirrored]))
        end
    end
    sort!(tsupp)
    unique!(tsupp)

    monomials = spectralgap.generate_mons(N, lso)
    monomials = spectralgap.filter_mons(monomials, tsupp, H, N)
    return (
        N=N,
        d=d,
        lso=lso,
        gamma=gamma,
        supports=supports,
        coefficients=coefficients,
        H=H,
        basis=basis,
        gbasis=gbasis,
        lb=lb,
        lgb=lgb,
        tsupp=tsupp,
        monomials=monomials,
    )
end

function expected_rows(source, problem::ExactRayProblem, spectralgap)
    dimensions = [source.lb; source.lgb]
    [block.dimension for block in problem.psd_blocks] == dimensions ||
        error("MOF PSD blocks do not match source basis dimensions")
    rows = [Dict{Int,BigRational}() for _ in source.tsupp]

    for block in eachindex(source.basis)
        for j in 1:source.lb[block], k in j:source.lb[block]
            word, coefficient = spectralgap.reduce!(
                [source.basis[block][j][1]; source.basis[block][k][1]],
                source.N;
                realify=true,
            )
            iszero(coefficient) && continue
            state_factors = [
                source.basis[block][j][2];
                source.basis[block][k][2]
            ]
            key = isempty(word) ? sort(state_factors) :
                  sort([state_factors; [word]])
            row_index = spectralgap.bfind(source.tsupp, key)
            ordinal = block_ordinal(problem.psd_blocks[block], j, k)
            add_coefficient!(
                rows,
                row_index,
                ordinal,
                j == k ? coefficient : 2coefficient,
            )
        end
    end

    gap_offset = length(source.basis)
    for block in eachindex(source.gbasis)
        for i in 1:source.lgb[block], j in i:source.lgb[block]
            words, coefficients = spectralgap.PSDstate_entry(
                source.gbasis[block][i][1],
                source.gbasis[block][j][1],
                source.H,
                source.N;
                realify=true,
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
                model="Ising",
            ) && !spectralgap.isz(
                source.gbasis[block][j][1];
                model="Ising",
            )
                mirrored = [
                    spectralgap.reduce_mirror(
                        source.gbasis[block][i][1],
                        source.N,
                    ),
                    spectralgap.reduce_mirror(
                        source.gbasis[block][j][1],
                        source.N,
                    ),
                ]
                add_coefficient!(
                    rows,
                    spectralgap.bfind(
                        source.tsupp,
                        sort([state_factors; mirrored]),
                    ),
                    ordinal,
                    factor * source.gamma,
                )
            end
        end
    end

    psd_ordinals = Set(
        ordinal for block in problem.psd_blocks for coordinate in block.coordinates
        for (ordinal, _) in coordinate.terms
    )
    scalar_ordinals = setdiff(
        collect(1:problem.variable_count),
        collect(psd_ordinals),
    )
    sort!(scalar_ordinals)
    length(scalar_ordinals) == length(source.monomials) + 1 ||
        error("MOF scalar-variable count does not match stationarity plus lambda")
    for (monomial_index, monomial) in enumerate(source.monomials)
        ordinal = scalar_ordinals[monomial_index]
        for (term_index, term) in enumerate(source.H.supp)
            word, coefficient = spectralgap.reduce!(
                [term; monomial],
                source.N,
            )
            iszero(imag(coefficient)) && continue
            row_index = spectralgap.bfind(source.tsupp, [word])
            add_coefficient!(
                rows,
                row_index,
                ordinal,
                source.H.coe[term_index] * imag(coefficient),
            )
        end
    end
    lambda_ordinal = last(scalar_ordinals)
    add_coefficient!(rows, 1, lambda_ordinal, 1)
    return rows, scalar_ordinals, lambda_ordinal
end

"""
    audit_tfim_source_assembly(spectralgap, repo_root, spectralgap_root, model_path)

Recompute every affine coefficient of the locked TFIM N=9, g=1/2, d=2,
lso=6, gamma=201/800 sign- and reflection-symmetric source formulation and
compare it exactly with the rational reconstruction of the exported MOF.
No optimizer is constructed or invoked.
"""
function audit_tfim_source_assembly(
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
    source = tfim_source_objects(spectralgap)
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
        coefficient for row in problem.equalities for (_, coefficient) in row.terms
    ))
    expected_coefficients == actual_coefficients ||
        error("source/MOF coefficient inventories differ")
    return (
        source_gate=gate,
        hamiltonian_supports=source.supports,
        hamiltonian_coefficients=source.coefficients,
        pauli_encoding="3i-2=X_i;3i-1=Y_i;3i=Z_i",
        boundary="open-chain-N9",
        state_class=(
            sign_symmetry="global-X spin flip; odd Y/Z parity removed by isz",
            spatial_symmetry="open-chain reflection via reduce_mirror",
            semantics="symmetry-restricted infinite-volume KMS ground-state state-polynomial relaxation",
        ),
        matrix_orientation=(
            source_entry="upper triangle j<=k",
            mof_storage="MOI PositiveSemidefiniteConeTriangle trimap",
            off_diagonal_factor="2 applied once during source affine assembly",
        ),
        gamma=201 // 800,
        basis_dimensions=source.lb,
        gap_basis_dimensions=source.lgb,
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
