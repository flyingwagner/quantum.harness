#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess
import JuMP

const MOI = JuMP.MOI

function main(args=ARGS)
    length(args) == 2 || error(
        "usage: normalize_mof_recession.jl INPUT.mof.json[.gz] OUTPUT.mof.json[.gz]",
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
    summary = normalize_recession_problem!(model)
    variables_after = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    variables_after == variables_before || error("variable inventory changed")
    MOI.write_to_file(model, String(output_path))

    println("variables\t", length(variables_after))
    println("original_objective_sense\t", summary.original_sense)
    println("objective_terms\t", summary.objective_terms)
    println("homogeneous_scalar_equalities\t", summary.scalar_equalities)
    println("homogeneous_zero_coordinates\t", summary.zero_coordinates)
    println("homogeneous_affine_psd_blocks\t", summary.affine_psd_blocks)
    println("direct_psd_blocks\t", summary.direct_psd_blocks)
    println("normalization_equalities\t1")
    println("new_objective_sense\t", MOI.get(model, MOI.ObjectiveSense()))
    println("optimization_invoked\tfalse")
    println("output\t", abspath(output_path))
    return 0
end

exit(main())
