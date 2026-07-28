#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess
import JuMP

const MOI = JuMP.MOI

function main(args=ARGS)
    length(args) in (3, 4) || error(
        "usage: equilibrate_mof_columns.jl INPUT.mof.json[.gz] OUTPUT.mof.json[.gz] SCALE.tsv [REFERENCE_RAY.tsv]",
    )
    input_path, output_path, scale_path = args
    abspath(input_path) == abspath(output_path) &&
        error("input and output paths must differ")
    ispath(output_path) && error("refusing to overwrite output: $output_path")
    ispath(scale_path) && error("refusing to overwrite scale map: $scale_path")
    model = MOI.FileFormats.Model(filename=String(input_path))
    MOI.read_from_file(model, String(input_path))
    variables_before = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    if length(args) == 4
        ray_indices, ray_values = read_ray_values(args[4])
        ray_indices == [variable.value for variable in variables_before] ||
            error("reference-ray MOI indices differ from the model")
        equilibration = ray_equilibration(model, ray_values)
        mode = "reference-ray-power-of-two"
    else
        equilibration = column_equilibration(model)
        mode = "column-power-of-two"
    end
    apply_variable_scaling!(model, equilibration)
    row_equilibration = equilibrate_rows!(model)
    variables_after = sort(
        MOI.get(model, MOI.ListOfVariableIndices());
        by=variable -> variable.value,
    )
    variables_after == variables_before || error("variable inventory changed")
    MOI.write_to_file(model, String(output_path))
    open(scale_path, "w") do io
        println(io, "ordinal\tmoi_index\tname\tcolumn_max\texponent\tx_per_z")
        for (ordinal, variable) in enumerate(equilibration.variables)
            name = MOI.get(model, MOI.VariableName(), variable)
            println(
                io,
                ordinal,
                '\t',
                variable.value,
                '\t',
                name,
                '\t',
                repr(equilibration.maxima[ordinal]),
                '\t',
                equilibration.exponents[ordinal],
                '\t',
                repr(equilibration.scales[ordinal]),
            )
        end
    end
    final_columns = column_equilibration(model)
    transformed_maxima = filter(!iszero, final_columns.maxima)
    println("variables\t", length(variables_after))
    println("mode\t", mode)
    println("direct_psd_blocks\t", length(equilibration.direct_psd_groups))
    println("minimum_scale\t", minimum(equilibration.scales))
    println("maximum_scale\t", maximum(equilibration.scales))
    println("minimum_exponent\t", minimum(equilibration.exponents))
    println("maximum_exponent\t", maximum(equilibration.exponents))
    println("row_scalar_equalities\t", row_equilibration.scalar_equalities)
    println("row_zero_coordinates\t", row_equilibration.zero_coordinates)
    println("row_affine_psd_blocks\t", row_equilibration.affine_psd_blocks)
    println("row_minimum_exponent\t", row_equilibration.minimum_exponent)
    println("row_maximum_exponent\t", row_equilibration.maximum_exponent)
    println("objective_exponent\t", row_equilibration.objective_exponent)
    println("minimum_nonzero_transformed_column_max\t", minimum(transformed_maxima))
    println("maximum_transformed_column_max\t", maximum(transformed_maxima))
    println("output\t", abspath(output_path))
    println("scale_map\t", abspath(scale_path))
    return 0
end

exit(main())
