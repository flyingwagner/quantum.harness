#!/usr/bin/env julia

using SHA

include(joinpath(@__DIR__, "..", "src", "SquareJ1J2Prototype.jl"))
using .SquareJ1J2Prototype
include(joinpath(@__DIR__, "..", "src", "GenericGapModel.jl"))
using .GenericGapModel
include(joinpath(@__DIR__, "..", "src", "CoreMGK.jl"))
include(joinpath(@__DIR__, "..", "src", "SharedCoreWire.jl"))
include(joinpath(@__DIR__, "..", "src", "SquareCoreInventory.jl"))
using .SquareCoreInventory
include(joinpath(@__DIR__, "..", "src", "SquareGapConic.jl"))
using .SquareGapConic
include(joinpath(@__DIR__, "..", "src", "SquareStatusEnvelope.jl"))
using .SquareStatusEnvelope

file_sha256(path::AbstractString) =
    open(path, "r") do io
        bytes2hex(sha256(io))
    end

function progress(phase::AbstractString)
    println("progress\t", phase)
    flush(stdout)
    return nothing
end

function core_progress(event)
    println(
        "progress\t",
        event.phase,
        '\t',
        event.completed,
        '/',
        event.total,
    )
    flush(stdout)
    return nothing
end

function git_object(
    repo::AbstractString,
    git_dir::AbstractString,
    expression::AbstractString,
)
    value = readchomp(
        `git --git-dir=$git_dir --work-tree=$repo rev-parse $expression`,
    )
    occursin(r"^[0-9a-f]{40}$", value) ||
        error("git returned a noncanonical object ID")
    return value
end

function require_clean_source(repo::AbstractString, git_dir::AbstractString)
    changes = readchomp(
        `git --git-dir=$git_dir --work-tree=$repo status --porcelain --untracked-files=no`,
    )
    isempty(changes) ||
        error("tracked source differs from the source commit:\n$changes")
    return nothing
end

function main(args=ARGS)
    length(args) == 8 || error(
        "usage: emit_square_status_envelope.jl CORE.aicore CORE.aicoreenv MODEL.mof.json[.gz] REPO_ROOT GIT_DIR PROJECT.toml MANIFEST.toml OUTPUT.aisqstatus",
    )
    core_path, core_envelope_path, mof_path, repo_root, git_dir,
        project_path, manifest_path, output_path = args
    isdir(repo_root) || error("repository root does not exist")
    isdir(git_dir) || error("Git directory does not exist")
    require_clean_source(repo_root, git_dir)
    source_commit = git_object(repo_root, git_dir, "HEAD")
    source_tree = git_object(repo_root, git_dir, "HEAD^{tree}")

    problem = GapProblem(
        square_patch_geometry(1),
        square_j1j2_model(1 // 2),
        1 // 10,
        2;
        basis_mode=:structured,
        basis_spec=StructuredBasisSpec(:one_symbol_lift, 1),
    )
    progress("core_source_rebuild_start")
    expected_core = build_square_core_inventory(
        problem;
        progress_every=50_000,
        progress_callback=core_progress,
    )
    progress("core_source_rebuild_complete")
    progress("core_artifact_validation_start")
    parsed_core =
        validate_square_core_inventory(core_path, core_envelope_path)
    progress("core_artifact_validation_complete")
    source_rebuild_match =
        expected_core.math_sha256 == parsed_core.math_sha256 &&
        expected_core.envelope_sha256 == parsed_core.envelope_sha256
    source_rebuild_match ||
        error("supplied core artifacts differ from the current source rebuild")
    core_validation = merge(
        parsed_core,
        (source_rebuild_match=true,),
    )
    expected_core = nothing
    GC.gc()

    progress("conic_source_rebuild_start")
    plan = build_square_conic_plan(problem)
    progress("conic_source_rebuild_complete")
    progress("mof_exact_replay_start")
    mof_audit = audit_rendered_mof(plan, mof_path)
    progress("mof_exact_replay_complete")
    identity = (
        L=1,
        d=2,
        g=BigInt(1) // BigInt(2),
        gamma=plan.gamma,
        state_class=plan.source.state_class,
        problem_sha256=plan.source.source_plan.problem_sha256,
        positive_basis_sha256=plan.source.positive_basis.sha256,
        gap_basis_sha256=plan.source.gap_basis.sha256,
        scalar_rows=length(plan.rows),
        affine_equalities=1 + length(plan.stationarity),
        psd_real_dimensions=[
            block.real_dimension
            for block in plan.psd_blocks
        ],
    )
    result = build_square_status_envelope(
        identity;
        source_commit=source_commit,
        source_tree=source_tree,
        project_sha256=file_sha256(project_path),
        manifest_sha256=file_sha256(manifest_path),
        core_validation=core_validation,
        mof_sha256=file_sha256(mof_path),
        mof_audit=mof_audit,
    )
    progress("status_envelope_write")
    write_square_status_envelope(output_path, result)
    println("source_commit\t", source_commit)
    println("source_tree\t", source_tree)
    println("core_math_sha256\t", core_validation.math_sha256)
    println("core_envelope_sha256\t", core_validation.envelope_sha256)
    println("core_source_rebuild_match\t", core_validation.source_rebuild_match)
    println("mof_sha256\t", file_sha256(mof_path))
    println("exact_coefficient_match\t", mof_audit.exact_coefficient_match)
    println("status\tunsolved")
    println("optimizer_invoked\tfalse")
    println("envelope_sha256\t", result.sha256)
    println("output\t", abspath(output_path))
    return 0
end

exit(main())
