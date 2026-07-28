#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "SquareJ1J2Prototype.jl"))
using .SquareJ1J2Prototype
include(joinpath(@__DIR__, "..", "src", "GenericGapModel.jl"))
using .GenericGapModel
include(joinpath(@__DIR__, "..", "src", "CoreMGK.jl"))
include(joinpath(@__DIR__, "..", "src", "SharedCoreWire.jl"))
include(joinpath(@__DIR__, "..", "src", "SquareGapConic.jl"))
using .SquareGapConic

function main(args=ARGS)
    length(args) == 1 || error(
        "usage: audit_square_conic_mof.jl MODEL.mof.json[.gz]",
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
    result = audit_rendered_mof(plan, args[1])
    println("problem_sha256\t", plan.source.source_plan.problem_sha256)
    println("positive_basis_sha256\t", plan.source.positive_basis.sha256)
    println("gap_basis_sha256\t", plan.source.gap_basis.sha256)
    for key in propertynames(result)
        println(key, '\t', getproperty(result, key))
    end
    return 0
end

exit(main())
