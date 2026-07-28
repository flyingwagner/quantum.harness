#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "SquareJ1J2Prototype.jl"))
include(joinpath(@__DIR__, "..", "src", "GenericGapModel.jl"))
include(joinpath(@__DIR__, "..", "src", "CoreMGK.jl"))
include(joinpath(@__DIR__, "..", "src", "SharedCoreWire.jl"))
include(joinpath(@__DIR__, "..", "src", "SquareCoreInventory.jl"))
using .SquareCoreInventory

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: validate_square_core_inventory.jl MATH.aicore ENVELOPE.aicoreenv",
    )
    result = validate_square_core_inventory(args[1], args[2])
    for key in propertynames(result)
        println(key, '\t', getproperty(result, key))
    end
    return 0
end

exit(main())
