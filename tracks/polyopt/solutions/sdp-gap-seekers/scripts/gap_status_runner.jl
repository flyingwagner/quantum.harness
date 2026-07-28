#!/usr/bin/env julia

include(joinpath(@__DIR__, "gap_status_runner_lib.jl"))

try
    exit(GapStatusRunner.main())
catch err
    bt = catch_backtrace()
    showerror(stderr, err, bt)
    println(stderr)
    exit(64)
end
