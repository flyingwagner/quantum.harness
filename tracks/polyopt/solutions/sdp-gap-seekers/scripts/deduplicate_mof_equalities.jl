#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess
import JuMP

const MOI = JuMP.MOI

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: deduplicate_mof_equalities.jl INPUT.mof.json[.gz] OUTPUT.mof.json[.gz]",
    )
    abspath(args[1]) == abspath(args[2]) && error("input and output paths must differ")
    isfile(args[2]) && error("refusing to overwrite existing output: $(args[2])")
    model = MOI.FileFormats.Model(filename=String(args[1]))
    MOI.read_from_file(model, String(args[1]))
    variables_before = length(MOI.get(model, MOI.ListOfVariableIndices()))
    summary = deduplicate_affine_equalities!(model)
    variables_after = length(MOI.get(model, MOI.ListOfVariableIndices()))
    variables_before == variables_after || error("variable inventory changed")
    MOI.write_to_file(model, String(args[2]))
    println("variables\t", variables_after)
    println("scalar_equalities_removed\t", summary.scalar_equalities_removed)
    println("vector_zero_coordinates_removed\t", summary.vector_zero_coordinates_removed)
    println("total_equalities_removed\t", summary.total_removed)
    println("output\t", abspath(args[2]))
    return 0
end

exit(main())
