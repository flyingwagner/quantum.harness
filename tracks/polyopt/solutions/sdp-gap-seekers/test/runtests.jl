using Test

include(joinpath(@__DIR__, "..", "src", "SquareJ1J2Prototype.jl"))
using .SquareJ1J2Prototype

@testset "square patch geometry" begin
    for L in 1:4
        patch = square_patch(L; g=1//2)
        side = 2L + 1
        @test length(patch.sites) == side^2
        @test length(patch.inner_ids) == (side - 2)^2
        @test count(b -> b.kind == :J1, patch.bonds) == 2side * (side - 1)
        @test count(b -> b.kind == :J2, patch.bonds) == 2(side - 1)^2
        @test validate_inner_buffer(patch)
        @test length(unique((b.kind, b.i, b.j) for b in patch.bonds)) ==
              length(patch.bonds)
    end
end

@testset "status runner static safety gates" begin
    solution_dir = normpath(joinpath(@__DIR__, ".."))
    patch_text = read(joinpath(solution_dir, "spectralgap_a1171c9.patch"), String)
    prepare_text = read(joinpath(solution_dir, "scripts", "prepare_gap_status_run.sh"), String)
    array_text = read(joinpath(solution_dir, "scripts", "gap_status_array.sbatch"), String)
    finalize_text = read(joinpath(solution_dir, "scripts", "finalize_gap_status_run.sh"), String)

    @test length(findall("_objv = try", patch_text)) == 4
    @test length(findall("_objv_available = !isnothing(_objv)", patch_text)) == 4
    @test occursin("_lambda_value = try", patch_text)
    @test occursin(
        "_lambda_value_available = !isnothing(_lambda_value)",
        patch_text,
    )
    lambda_assignment = findfirst("_lambda_value = try", patch_text)
    lambda_availability =
        findfirst("_lambda_value_available = !isnothing(_lambda_value)", patch_text)
    @test !isnothing(lambda_assignment)
    @test !isnothing(lambda_availability)
    @test first(something(lambda_assignment)) < first(something(lambda_availability))
    @test !occursin("+    _objv = objective_value(model)", patch_text)

    @test occursin(raw"mkdir \"$RUN_DIR\" ||", prepare_text)
    @test occursin(
        raw"cp -a \"$SPECTRALGAP_PATH\" \"$FROZEN_SPECTRALGAP\"",
        prepare_text,
    )
    @test occursin("--startup-file=no --history-file=no", prepare_text)
    @test occursin(raw"--project=\"$FREEZE_DIR\"", prepare_text)

    @test occursin("gap-status-array-contract-v1", array_text)
    @test occursin(
        raw"--run-spec \"$RUN_SPEC\" --cell-index \"$TASK_ID\"",
        array_text,
    )
    @test occursin(raw"--project=\"$FREEZE_DIR\"", array_text)
    @test occursin("GAP_STATUS_EXPECTED_JULIA_SHA256", array_text)
    @test !occursin("LD_LIBRARY_PATH", array_text)
    @test !occursin(raw"--model \"$MODEL\" --gamma \"$GAMMA\"", array_text)

    @test occursin("sacct-cell-outcomes.tsv", finalize_text)
    @test occursin("FINAL_OUTCOME=success", finalize_text)
    @test occursin("FINAL_OUTCOME=failed", finalize_text)
    @test occursin("cell directory set mismatch", finalize_text)
    @test occursin("source gate absent", finalize_text)
end

@testset "Pauli canonicalization" begin
    one, identity = pauli_word([(1, :X), (1, :X)])
    @test one == 1
    @test isempty(identity.ops)

    phase_xy, xy = pauli_word([(1, :X), (1, :Y)])
    phase_yx, yx = pauli_word([(1, :Y), (1, :X)])
    _, z = pauli_word([(1, :Z)])
    @test phase_xy == im
    @test phase_yx == -im
    @test xy == z == yx

    phase_12, word_12 = pauli_word([(1, :X), (2, :Y)])
    phase_21, word_21 = pauli_word([(2, :Y), (1, :X)])
    @test phase_12 == phase_21 == 1
    @test word_12 == word_21

    _, x = pauli_word([(1, :X)])
    _, y = pauli_word([(1, :Y)])
    phase_left, left = multiply_words(x, y)
    phase_right, right = multiply_words(y, x)
    @test phase_left == im
    @test phase_right == -im
    @test left == right == z
end

@testset "bare Pauli basis counts" begin
    for nsites in 0:6, d in 0:4
        words = enumerate_pauli_words(nsites, d)
        @test length(words) == operator_word_count(nsites, d)
        @test length(unique(words)) == length(words)
    end
    @test operator_word_count(9, 2) == 352
    @test operator_word_count(25, 2) == 2776
end

@testset "full state-polynomial formal counts" begin
    @test full_state_basis_count_by_degree(1, 0) == BigInt[1]
    @test full_state_basis_count_by_degree(1, 1) == BigInt[1, 6]
    @test full_state_basis_count_by_degree(1, 2) == BigInt[1, 6, 15]
    @test full_state_basis_count(9, 2) == 1810
    @test one_symbol_lift_count(9, 2) == 703

    for nsites in 1:8
        counts = [full_state_basis_count(nsites, d) for d in 0:4]
        @test issorted(counts)
    end
end

@testset "storage estimates" begin
    @test dense_complex_matrix_bytes(10) == 1600
    @test real_embedding_matrix_bytes(10) == 3200
end

include(joinpath(@__DIR__, "..", "src", "LocalSpinIdentities.jl"))
using .LocalSpinIdentities

@testset "exact local spin identities" begin
    checks = local_identity_checks()
    for (name, result) in checks
        if result isa Bool
            @test result
        end
    end
    @test checks["bond_projector_traces"] == (1, 3)
    @test checks["triangle_projector_traces"] == (4, 4)
    @test checks["plaquette_projector_traces"] == (2, 9, 5)
    @test checks["joint_projector_traces"] == (1, 3, 3, 1, 3, 5)
end

include(joinpath(@__DIR__, "..", "src", "GenericGapModel.jl"))
using .GenericGapModel

@testset "generic solver-free problem adapter" begin
    for (L, expected_j1_bonds, expected_j2_bonds) in (
        (1, 12, 8),
        (2, 40, 32),
        (3, 84, 72),
    )
        patch = square_patch_geometry(L)
        model = square_j1j2_model(1//2)
        @test validate_model_buffer(model, patch)
        terms = instantiate_terms(model, patch)
        @test count(term -> term.tag == :J1, terms) == 3expected_j1_bonds
        @test count(term -> term.tag == :J2, terms) == 3expected_j2_bonds
        @test all(iszero ∘ imag ∘ (term -> term.coefficient), terms)
        @test all(
            term -> real(term.coefficient) == 1//4,
            filter(term -> term.tag == :J1, terms),
        )
        @test all(
            term -> real(term.coefficient) == 1//8,
            filter(term -> term.tag == :J2, terms),
        )
    end

    patch = square_patch_geometry(1)
    model = square_j1j2_model(1//2)
    integer_model = square_j1j2_model(0)
    @test all(
        term -> term.coefficient isa ComplexF64,
        instantiate_terms(integer_model, patch),
    )
    structured_problem = GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:one_symbol,
    )
    structured_plan = assembly_plan(structured_problem)
    @test structured_plan.local_terms == 60
    @test structured_plan.positive_basis_dimension == 703
    @test structured_plan.gap_basis_dimension == 7
    @test !structured_plan.symmetry_declared
    @test structured_plan.problem_sha256 ==
          assembly_plan(structured_problem).problem_sha256

    complete_problem = GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:full_count_only,
    )
    complete_plan = assembly_plan(complete_problem)
    @test complete_plan.positive_basis_dimension == 1810
    @test complete_plan.gap_basis_dimension == 7
    @test complete_plan.problem_sha256 != structured_plan.problem_sha256

    symmetric_problem = GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:one_symbol,
        symmetry=ExplicitStateSymmetry("D4", ["C4", "mirror"]),
    )
    symmetric_plan = assembly_plan(symmetric_problem)
    @test symmetric_plan.symmetry_declared
    @test symmetric_plan.problem_sha256 != structured_plan.problem_sha256

    supports, coefficients = legacy_ncpoly_data(structured_problem)
    @test length(supports) == length(coefficients) == 60
    @test all(length(support) == 2 for support in supports)
    @test count(==(Float64(1//4)), coefficients) == 36
    @test count(==(Float64(1//8)), coefficients) == 24

    changed_model = square_j1j2_model(107//200)
    changed_patch = square_patch_geometry(1)
    changed_problem = GapProblem(changed_patch, changed_model, 1//10, 2)
    @test assembly_plan(changed_problem).problem_sha256 !=
          structured_plan.problem_sha256

    bad_sites = [Site(0, 0)]
    bad_patch = LocalPatch("bad-unbuffered", 0, bad_sites, Dict(Site(0, 0) => 1), [1])
    @test !validate_model_buffer(model, bad_patch)
    @test_throws ArgumentError GapProblem(bad_patch, model, 0//1, 2)
end

include(joinpath(@__DIR__, "..", "src", "SmallEDOracle.jl"))
using .SmallEDOracle

@testset "small finite-patch ED construction oracle" begin
    comparison = compare_hamiltonian_builders(1; g=1//2)
    @test comparison.max_builder_difference == 0
    @test comparison.hermiticity_error == 0
    @test comparison.trace == 0
end

include(joinpath(@__DIR__, "..", "scripts", "gap_status_runner_lib.jl"))
using .GapStatusRunner

@testset "solver-free status runner contract" begin
    points = allowed_points()
    @test length(points) == 6
    @test point_id.(points) == [
        "tfim-n9-g0p5-d2-lso6-gamma0p25",
        "tfim-n9-g0p5-d2-lso6-gamma0p26",
        "kagome-n13-d3-lso5-gamma1",
        "kagome-n13-d3-lso5-gamma1p2",
        "kagome-n13-d3-lso5-gamma1p26",
        "kagome-n13-d3-lso5-gamma1p28",
    ]

    tfim = parse_point("tfim", "0.25")
    kagome = parse_point("kagome", "1.28")
    @test_throws ArgumentError parse_point("tfim", "0.30")
    @test_throws ArgumentError parse_point("tfim", "0.25"; N_text="10")
    @test_throws ArgumentError parse_point("kagome", "1.28"; d_text="4")
    @test_throws ArgumentError parse_point("kagome", "1.28"; g_text="0.5")

    tfim_supports, tfim_coefficients = build_hamiltonian_data(tfim)
    kagome_supports, kagome_coefficients = build_hamiltonian_data(kagome)
    @test length(tfim_supports) == length(tfim_coefficients) == 17
    @test count(==(-1 // 1), tfim_coefficients) == 8
    @test count(==(1 // 2), tfim_coefficients) == 9
    @test length(kagome_supports) == length(kagome_coefficients) == 54
    @test all(==(1 // 4), kagome_coefficients)
    @test hamiltonian_metadata(tfim).fingerprint ==
          "6b9818740fdc9299abbcd592540ad4630b572a94e6d2345cfd36120e3016a221"
    @test hamiltonian_metadata(kagome).fingerprint ==
          "d8276e59e709c7356c3a2ae25eae6de0bbee8920dfe9413fc7ff178c2f96b4ab"
    @test basis_metadata(tfim).fingerprint ==
          "47f03764c25af510833463339fa90acaf41f6a79bb15729a1b73eaed57bc4e54"
    @test basis_metadata(kagome).fingerprint ==
          "adfbe4ef28077e71039ea1054bab956331297f14942b0b0b444aaba96905060e"

    named = adapt_certify_result(
        (
            flag=0,
            termination=:SLOW_PROGRESS,
            primal=:INFEASIBILITY_CERTIFICATE,
            dual=:NO_SOLUTION,
            objective=1.25,
        ),
    )
    @test named.adapter == "namedtuple-v1"
    @test isequal(named.flag, 0)
    @test named.termination == "SLOW_PROGRESS"
    @test named.primal == "INFEASIBILITY_CERTIFICATE"
    @test named.dual == "NO_SOLUTION"
    @test named.objective == 1.25
    @test named.objective_availability == "available"

    objective_missing = adapt_certify_result(
        (
            flag=0,
            termination=:INFEASIBLE,
            primal=:NO_SOLUTION,
            dual=:INFEASIBILITY_CERTIFICATE,
            objective=nothing,
            objective_available=false,
        ),
    )
    @test objective_missing.termination == "INFEASIBLE"
    @test objective_missing.dual == "INFEASIBILITY_CERTIFICATE"
    @test isnothing(objective_missing.objective)
    @test objective_missing.objective_availability == "unavailable"

    legacy = adapt_certify_result(0)
    @test legacy.adapter == "legacy-int-explicit-unavailable-status"
    @test isequal(legacy.flag, 0)
    @test legacy.termination == "UNAVAILABLE_LEGACY_INT"
    @test legacy.primal == "UNAVAILABLE_LEGACY_INT"
    @test legacy.dual == "UNAVAILABLE_LEGACY_INT"
    @test isnothing(legacy.objective)
    @test legacy.objective_availability == "unavailable-legacy-api"
    @test_throws ArgumentError adapt_certify_result((flag=1, termination=:OPTIMAL))

    fixture_record = execute_point(
        tfim,
        _ -> (
            flag=1,
            termination=:OPTIMAL,
            primal=:FEASIBLE_POINT,
            dual=:FEASIBLE_POINT,
            objective=0.0,
        );
        source=(mode="fixture",),
    )
    @test isnothing(fixture_record.exception)
    @test isequal(fixture_record.flag, 1)
    @test fixture_record.termination == "OPTIMAL"
    @test fixture_record.residual.availability == "unavailable"
    @test fixture_record.witness.availability == "unavailable"

    exception_record = execute_point(
        kagome,
        _ -> error("fixture failure");
        source=(mode="fixture",),
    )
    @test !isnothing(exception_record.exception)
    @test isnothing(exception_record.flag)
    @test exception_record.termination == "UNAVAILABLE_EXCEPTION"
    @test exception_record.exception.message == "fixture failure"
    @test occursin("fixture failure", exception_record.exception.stacktrace)

    encoded = json_encode(fixture_record)
    @test startswith(encoded, "{")
    @test occursin("\"termination\":\"OPTIMAL\"", encoded)
    @test occursin("\"witness\":{\"availability\":\"unavailable\"", encoded)

    dry = build_dry_run_record(tfim, normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..")), joinpath(@__DIR__, "missing-env"))
    @test dry.solver_called == false
    @test dry.environment.lock_state == "unlocked-manifest-missing"
    @test !isdefined(Main, :SpectralGap)

    mktempdir() do directory
        spec_path = joinpath(directory, "run_spec.json")
        run_id = "fixture-run"
        open(spec_path, "w") do io
            print(
                io,
                """
                {
                  "schema_version": "gap-status-run-spec-v1",
                  "run_id": "$run_id",
                  "cells": [
                    {"cell_id": "01-tfim-gamma-0p25", "params": {"model": "tfim", "gamma": "0.25"}},
                    {"cell_id": "02-tfim-gamma-0p26", "params": {"model": "tfim", "gamma": "0.26"}},
                    {"cell_id": "03-kagome-gamma-1", "params": {"model": "kagome", "gamma": "1"}},
                    {"cell_id": "04-kagome-gamma-1p2", "params": {"model": "kagome", "gamma": "1.2"}},
                    {"cell_id": "05-kagome-gamma-1p26", "params": {"model": "kagome", "gamma": "1.26"}},
                    {"cell_id": "06-kagome-gamma-1p28", "params": {"model": "kagome", "gamma": "1.28"}}
                  ]
                }
                """,
            )
        end
        selected, selected_run = parse_run_spec_point(spec_path, 6)
        @test point_id(selected) == "kagome-n13-d3-lso5-gamma1p28"
        @test selected_run == run_id
        @test_throws ArgumentError parse_run_spec_point(spec_path, 7)

        tampered = replace(read(spec_path, String), "\"1.28\"" => "\"1.29\"")
        open(spec_path, "w") do io
            write(io, tampered)
        end
        @test_throws ArgumentError parse_run_spec_point(spec_path, 6)
    end
end
