#!/usr/bin/env julia

import JuMP
import MosekTools

const MOI = JuMP.MOI

function inventory(model)
    constraints = [
        (
            function_type=string(function_type),
            set_type=string(set_type),
            count=length(MOI.get(
                model,
                MOI.ListOfConstraintIndices{function_type,set_type}(),
            )),
        ) for (function_type, set_type) in
        MOI.get(model, MOI.ListOfConstraintTypesPresent())
    ]
    return (
        variables=length(MOI.get(model, MOI.ListOfVariableIndices())),
        constraints=constraints,
    )
end

function load_model(path)
    model = MOI.FileFormats.Model(filename=String(path))
    MOI.read_from_file(model, String(path))
    return model
end

function print_inventory(model)
    model_inventory = inventory(model)
    println("variables\t", model_inventory.variables)
    for row in model_inventory.constraints
        println("constraint_type\t", row.function_type, '\t', row.set_type, '\t', row.count)
    end
end

function safe_status(optimizer, attribute, fallback)
    return try
        MOI.get(optimizer, attribute)
    catch
        fallback
    end
end

function write_result(path, fields)
    open(path, "w") do io
        println(io, "key\tvalue")
        for (key, value) in pairs(fields)
            rendered = replace(string(value), '\n' => "\\n", '\t' => "\\t")
            println(io, key, '\t', rendered)
        end
    end
end

function solve(model_path, output_prefix)
    model = load_model(model_path)
    variables = sort(MOI.get(model, MOI.ListOfVariableIndices()); by=index -> index.value)
    optimizer = MOI.instantiate(MosekTools.Optimizer; with_bridge_type=Float64)
    MOI.set(optimizer, MOI.Silent(), false)
    threads = parse(Int, get(ENV, "MOF_SOLVE_THREADS", "1"))
    time_limit = parse(Float64, get(ENV, "MOF_SOLVE_MAX_TIME", "3600"))
    threads > 0 || error("MOF_SOLVE_THREADS must be positive")
    time_limit > 0 || error("MOF_SOLVE_MAX_TIME must be positive")
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS"),
        threads,
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME"),
        time_limit,
    )
    index_map = MOI.copy_to(optimizer, model)
    elapsed = @elapsed MOI.optimize!(optimizer)

    termination = safe_status(optimizer, MOI.TerminationStatus(), MOI.OTHER_ERROR)
    primal = safe_status(optimizer, MOI.PrimalStatus(), MOI.NO_SOLUTION)
    dual = safe_status(optimizer, MOI.DualStatus(), MOI.NO_SOLUTION)
    objective = safe_status(optimizer, MOI.ObjectiveValue(), NaN)
    raw_status = safe_status(optimizer, MOI.RawStatusString(), "unavailable")
    result_count = safe_status(optimizer, MOI.ResultCount(), 0)

    write_result(
        output_prefix * ".result.tsv",
        (
            model=abspath(model_path),
            mosektools_version=Base.pkgversion(MosekTools),
            threads=threads,
            time_limit_seconds=time_limit,
            walltime_seconds=elapsed,
            termination=termination,
            primal=primal,
            dual=dual,
            raw_status=raw_status,
            result_count=result_count,
            objective=objective,
        ),
    )
    open(output_prefix * ".variables.tsv", "w") do io
        println(io, "ordinal\tmoi_index\tname\tvalue")
        for (ordinal, variable) in enumerate(variables)
            value = try
                MOI.get(optimizer, MOI.VariablePrimal(), index_map[variable])
            catch
                NaN
            end
            name = try
                MOI.get(model, MOI.VariableName(), variable)
            catch
                ""
            end
            println(io, ordinal, '\t', variable.value, '\t', name, '\t', repr(value))
        end
    end
    println("result\t", abspath(output_prefix * ".result.tsv"))
    println("variables\t", abspath(output_prefix * ".variables.tsv"))
    return 0
end

function main(args=ARGS)
    if length(args) == 2 && args[1] == "--preflight"
        model = load_model(args[2])
        print_inventory(model)
        println("mosektools_version\t", Base.pkgversion(MosekTools))
        println("optimization_invoked\tfalse")
        return 0
    end
    length(args) == 2 || error(
        "usage: solve_exported_mof.jl [--preflight] MODEL.mof.json[.gz] OUTPUT_PREFIX",
    )
    return solve(args[1], args[2])
end

exit(main())
