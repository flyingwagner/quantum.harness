#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "TFIMSourceAudit.jl"))
using .TFIMSourceAudit
import SpectralGap

function render(value)
    if value isa AbstractVector
        return join(value, ",")
    elseif value isa Pair
        return string(first(value), "=>", last(value))
    elseif value isa Tuple || value isa NamedTuple
        return repr(value)
    end
    return string(value)
end

function main(args=ARGS)
    length(args) in (1, 2) || error(
        "usage: audit_tfim_source_assembly.jl TFIM.mof.json[.gz] [SPECTRALGAP_ROOT]",
    )
    repo_root = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
    spectralgap_root =
        length(args) == 2 ? abspath(args[2]) :
        joinpath(repo_root, ".external", "SpectralGap")
    result = audit_tfim_source_assembly(
        SpectralGap,
        repo_root,
        spectralgap_root,
        args[1],
    )
    println("source_commit\t", result.source_gate.commit)
    println("patch_sha256\t", result.source_gate.patch_sha256)
    println("hamiltonian_supports\t", repr(result.hamiltonian_supports))
    println("hamiltonian_coefficients\t", join(result.hamiltonian_coefficients, ','))
    println("pauli_encoding\t", result.pauli_encoding)
    println("boundary\t", result.boundary)
    println("state_class\t", repr(result.state_class))
    println("matrix_orientation\t", repr(result.matrix_orientation))
    println("gamma\t", result.gamma)
    println("basis_dimensions\t", render(result.basis_dimensions))
    println("gap_basis_dimensions\t", render(result.gap_basis_dimensions))
    println("stationarity_monomials\t", result.stationarity_monomials)
    println("lambda_ordinal\t", result.lambda_ordinal)
    println("variable_count\t", result.variable_count)
    println("equality_count\t", result.equality_count)
    println("equality_offsets_zero\t", result.equality_offsets_zero)
    println("psd_dimensions\t", render(result.psd_dimensions))
    println("objective_sense\t", result.objective_sense)
    println("objective\t", join(render.(result.objective), ','))
    println("coefficient_inventory\t", render(result.coefficient_inventory))
    println("affine_row_sha256\t", result.affine_row_sha256)
    println("rows_compared\t", result.rows_compared)
    println("coefficient_mismatches\t", result.coefficient_mismatches)
    println("optimizer_invoked\t", result.optimizer_invoked)
    println("source_assembly_equal\ttrue")
    return 0
end

exit(main())
