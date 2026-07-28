#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess

function main(args=ARGS)
    length(args) == 3 || error(
        "usage: backtransform_mof_ray.jl SCALE.tsv SCALED.variables.tsv ORIGINAL.variables.tsv",
    )
    scale_indices, scale_names, scales = read_scale_map(args[1])
    ray_indices, scaled_values = read_ray_values(args[2])
    ray_indices == scale_indices || error("ray/scale MOI indices differ")
    ray_lines = readlines(args[2])
    ray_names = [split(line, '\t'; keepempty=true)[3] for line in ray_lines[2:end]]
    ray_names == scale_names || error("ray/scale variable names differ")
    values = backtransform_ray(scaled_values, scales)
    ispath(args[3]) && error("refusing to overwrite output: $(args[3])")
    open(args[3], "w") do io
        println(io, "ordinal\tmoi_index\tname\tvalue")
        for ordinal in eachindex(values)
            println(
                io,
                ordinal,
                '\t',
                ray_indices[ordinal],
                '\t',
                ray_names[ordinal],
                '\t',
                repr(values[ordinal]),
            )
        end
    end
    println("variables\t", length(values))
    println("scaled_ray_scale\t", maximum(abs, scaled_values))
    println("original_ray_scale\t", maximum(abs, values))
    println("output\t", abspath(args[3]))
    return 0
end

exit(main())
