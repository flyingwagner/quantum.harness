#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "SquareJ1J2Prototype.jl"))
using .SquareJ1J2Prototype
include(joinpath(@__DIR__, "..", "src", "GenericGapModel.jl"))
using .GenericGapModel
include(joinpath(@__DIR__, "..", "src", "CoreMGK.jl"))
using .CoreMGK
include(joinpath(@__DIR__, "..", "src", "SharedCoreWire.jl"))
using .SharedCoreWire
include(joinpath(@__DIR__, "..", "src", "SquareCoreInventory.jl"))
using .SquareCoreInventory

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: emit_square_core_inventory.jl MATH.aicore ENVELOPE.aicoreenv",
    )
    problem = GapProblem(
        square_patch_geometry(1),
        square_j1j2_model(1 // 2),
        1 // 10,
        2;
        basis_mode=:structured,
        basis_spec=StructuredBasisSpec(:one_symbol_lift, 1),
    )
    result = write_square_core_inventory(problem, args[1], args[2])
    println("math_scope\tcore_mgk")
    println("gamma_in_math_bytes\tfalse")
    println("relaxation_id\t", result.relaxation_id)
    println("math_byte_count\t", length(result.math_bytes))
    println("math_sha256\t", result.math_sha256)
    println("envelope_byte_count\t", length(result.envelope_bytes))
    println("envelope_sha256\t", result.envelope_sha256)
    println("blocks\t2")
    println("scalar_rows\t", result.scalar_rows)
    println("wiring_pairs\t", result.wiring_pairs)
    println("component_records\t", result.component_records)
    println("nonzero_coefficients\t", result.nonzero_coefficients)
    println("optimizer_invoked\tfalse")
    println("math_output\t", abspath(args[1]))
    println("envelope_output\t", abspath(args[2]))
    return 0
end

exit(main())
