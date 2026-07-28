#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess
import JuMP

const MOI = JuMP.MOI

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: equilibrate_mof_rows.jl INPUT.mof.json[.gz] OUTPUT.mof.json[.gz]",
    )
    input_path, output_path = args
    abspath(input_path) == abspath(output_path) &&
        error("input and output paths must differ")
    ispath(output_path) && error("refusing to overwrite output: $output_path")
    model = MOI.FileFormats.Model(filename=String(input_path))
    MOI.read_from_file(model, String(input_path))
    variables_before = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    summary = equilibrate_rows!(model)
    variables_after = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    variables_after == variables_before || error("variable inventory changed")
    MOI.write_to_file(model, String(output_path))

    println("variables\t", length(variables_after))
    println("row_scalar_equalities\t", summary.scalar_equalities)
    println("row_zero_coordinates\t", summary.zero_coordinates)
    println("row_affine_psd_blocks\t", summary.affine_psd_blocks)
    println("row_minimum_exponent\t", summary.minimum_exponent)
    println("row_maximum_exponent\t", summary.maximum_exponent)
    println("objective_exponent\t", summary.objective_exponent)
    println("optimization_invoked\tfalse")
    println("output\t", abspath(output_path))
    return 0
end

exit(main())
