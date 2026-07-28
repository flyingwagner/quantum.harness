#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "SquareJ1J2Prototype.jl"))
using .SquareJ1J2Prototype
include(joinpath(@__DIR__, "..", "src", "GenericGapModel.jl"))
using .GenericGapModel
include(joinpath(@__DIR__, "..", "src", "CoreMGK.jl"))
using .CoreMGK
include(joinpath(@__DIR__, "..", "src", "SharedCoreWire.jl"))
using .SharedCoreWire
include(joinpath(@__DIR__, "..", "src", "SquareGapConic.jl"))
using .SquareGapConic

function main(args=ARGS)
    length(args) == 1 || error(
        "usage: emit_square_conic_mof.jl OUTPUT.mof.json[.gz]",
    )
    problem = GapProblem(
        square_patch_geometry(1),
        square_j1j2_model(1 // 2),
        1 // 10,
        2;
        basis_mode=:structured,
        basis_spec=StructuredBasisSpec(:one_symbol_lift, 1),
    )
    plan = build_square_conic_plan(problem)
    result = render_mof(plan, args[1])
    println("problem_sha256\t", plan.source.source_plan.problem_sha256)
    println("positive_basis_sha256\t", plan.source.positive_basis.sha256)
    println("gap_basis_sha256\t", plan.source.gap_basis.sha256)
    println("state_class\t", plan.source.state_class)
    println("gamma\t", plan.gamma)
    println("scalar_rows\t", length(plan.rows))
    println("normalization_rhs\t", plan.normalization.rhs)
    println(
        "stationarity_selector_entries\t",
        plan.stationarity_selector_entries,
    )
    println("stationarity_equalities\t", length(plan.stationarity))
    println(
        "stationarity_exact_duplicates_removed\t",
        plan.stationarity_exact_duplicates_removed,
    )
    println("psd_complex_dimensions\t", join(
        (block.complex_dimension for block in plan.psd_blocks),
        ',',
    ))
    println("psd_real_dimensions\t", join(result.psd_dimensions, ','))
    println("objective_sense\t", result.objective_sense)
    println("optimizer_invoked\t", result.optimizer_invoked)
    println("output\t", abspath(args[1]))
    return 0
end

exit(main())
