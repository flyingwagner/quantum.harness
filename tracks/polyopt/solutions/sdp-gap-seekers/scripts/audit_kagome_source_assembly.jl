#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "KagomeSourceAudit.jl"))
using .KagomeSourceAudit
import SpectralGap

render(value) = value isa AbstractVector ? join(value, ",") : string(value)

function main(args=ARGS)
    length(args) in (1, 2) || error(
        "usage: audit_kagome_source_assembly.jl " *
        "KAGOME.mof.json[.gz] [SPECTRALGAP_ROOT]",
    )
    repo_root = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
    spectralgap_root =
        length(args) == 2 ? abspath(args[2]) :
        joinpath(repo_root, ".external", "SpectralGap")
    result = audit_kagome_source_assembly(
        SpectralGap,
        repo_root,
        spectralgap_root,
        args[1],
    )
    println("source_commit\t", result.source_gate.commit)
    println("patch_sha256\t", result.source_gate.patch_sha256)
    println("hamiltonian_supports\t", repr(result.hamiltonian_supports))
    println(
        "hamiltonian_coefficients\t",
        join(result.hamiltonian_coefficients, ','),
    )
    println("pauli_encoding\t", result.pauli_encoding)
    println("patch\t", repr(result.patch))
    println("state_class\t", repr(result.state_class))
    println("matrix_orientation\t", repr(result.matrix_orientation))
    println("gamma\t", result.gamma)
    println("basis_dimensions\t", render(result.basis_dimensions))
    println("gap_basis_dimensions\t", render(result.gap_basis_dimensions))
    println(
        "strengthening_dimensions\t",
        render(result.strengthening_dimensions),
    )
    println("stationarity_monomials\t", result.stationarity_monomials)
    println("lambda_ordinal\t", result.lambda_ordinal)
    println("variable_count\t", result.variable_count)
    println("equality_count\t", result.equality_count)
    println("equality_offsets_zero\t", result.equality_offsets_zero)
    println("psd_dimensions\t", render(result.psd_dimensions))
    println("objective_sense\t", result.objective_sense)
    println("objective\t", join(string.(result.objective), ','))
    println(
        "coefficient_inventory\t",
        render(result.coefficient_inventory),
    )
    println("affine_row_sha256\t", result.affine_row_sha256)
    println("rows_compared\t", result.rows_compared)
    println("coefficient_mismatches\t", result.coefficient_mismatches)
    println("optimizer_invoked\t", result.optimizer_invoked)
    println("source_assembly_equal\ttrue")
    return 0
end

exit(main())
