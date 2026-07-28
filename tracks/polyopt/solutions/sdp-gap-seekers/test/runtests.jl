using Test
using SHA

include(joinpath(@__DIR__, "..", "scripts", "verify_gap_ray.jl"))
using .GapRayVerifier
include(joinpath(@__DIR__, "..", "src", "GapRayPostprocess.jl"))
using .GapRayPostprocess
include(joinpath(@__DIR__, "..", "src", "TFIMSourceAudit.jl"))
using .TFIMSourceAudit
include(joinpath(@__DIR__, "gap_ray_verifier_tests.jl"))

include(joinpath(@__DIR__, "..", "src", "SquareJ1J2Prototype.jl"))
using .SquareJ1J2Prototype

@testset "TFIM source-audit fail-closed row comparison" begin
    one = BigInt(1) // BigInt(1)
    half = BigInt(1) // BigInt(2)
    expected = [
        Dict{Int,Rational{BigInt}}(1 => one, 3 => -half),
        Dict{Int,Rational{BigInt}}(),
    ]
    same = [
        Dict{Int,Rational{BigInt}}(3 => -half, 1 => one),
        Dict{Int,Rational{BigInt}}(),
    ]
    changed_sign = [
        Dict{Int,Rational{BigInt}}(1 => one, 3 => half),
        Dict{Int,Rational{BigInt}}(),
    ]
    missing_row = same[1:1]
    @test isempty(row_mismatches(expected, same))
    @test row_mismatches(expected, changed_sign) == [1]
    @test row_mismatches(expected, missing_row) == [1, 2]

    audit_text = read(
        joinpath(@__DIR__, "..", "scripts", "audit_tfim_source_assembly.jl"),
        String,
    )
    @test occursin("source_assembly_equal\\ttrue", audit_text)
    @test occursin("optimizer_invoked", audit_text)
    @test !occursin("optimize!", audit_text)
end

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
    @test !occursin("+using Clarabel", patch_text)
    @test !occursin("_select_optimizer", patch_text)
    @test bytes2hex(sha256(codeunits(patch_text))) ==
          "332c0931ac810289aa3713af0948f259c01189270706af58b262d60d994d4abd"
    @test occursin("MSK_IPAR_PTF_WRITE_SOLUTIONS", patch_text)
    @test occursin("getsolutioninfo", patch_text)
    @test occursin("equality_residual_relative", patch_text)
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
    baseline_problem = GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:one_symbol,
    )
    baseline_plan = assembly_plan(baseline_problem)
    @test baseline_plan.local_terms == 60
    @test baseline_plan.positive_basis_dimension == 703
    @test baseline_plan.gap_basis_dimension == 7
    @test !baseline_plan.is_complete
    @test baseline_plan.positive_basis_sha256 === nothing
    @test baseline_plan.gap_basis_sha256 === nothing
    @test !baseline_plan.symmetry_declared
    @test baseline_plan.problem_sha256 ==
          assembly_plan(baseline_problem).problem_sha256

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
    @test complete_plan.is_complete
    @test complete_plan.positive_basis_sha256 === nothing
    @test complete_plan.gap_basis_sha256 === nothing
    @test complete_plan.problem_sha256 != baseline_plan.problem_sha256

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
    @test symmetric_plan.problem_sha256 != baseline_plan.problem_sha256
    joined_generator_problem = GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:one_symbol,
        symmetry=ExplicitStateSymmetry("D4", ["C4|mirror"]),
    )
    @test assembly_plan(joined_generator_problem).problem_sha256 !=
          symmetric_plan.problem_sha256

    supports, coefficients = legacy_ncpoly_data(baseline_problem)
    @test length(supports) == length(coefficients) == 60
    @test all(length(support) == 2 for support in supports)
    @test count(==(Float64(1//4)), coefficients) == 36
    @test count(==(Float64(1//8)), coefficients) == 24

    changed_model = square_j1j2_model(107//200)
    changed_patch = square_patch_geometry(1)
    changed_problem = GapProblem(changed_patch, changed_model, 1//10, 2)
    @test assembly_plan(changed_problem).problem_sha256 !=
          baseline_plan.problem_sha256

    bad_sites = [Site(0, 0)]
    bad_patch = LocalPatch("bad-unbuffered", 0, bad_sites, Dict(Site(0, 0) => 1), [1])
    @test !validate_model_buffer(model, bad_patch)
    @test_throws ArgumentError GapProblem(bad_patch, model, 0//1, 2)
end

@testset "structured basis manifests" begin
    patch = square_patch_geometry(1)
    model = square_j1j2_model(1//2)
    spec = StructuredBasisSpec(:one_symbol_lift, 1)

    @test_throws ArgumentError GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:structured,
    )
    @test_throws ArgumentError GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:one_symbol,
        basis_spec=spec,
    )

    problem = GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:structured,
        basis_spec=spec,
    )
    positive = basis_manifest(problem, :positive)
    gap = basis_manifest(problem, :gap)
    plan = assembly_plan(problem)

    @test positive.role == :positive
    @test gap.role == :gap
    @test positive.family == gap.family == :one_symbol_lift
    @test positive.family_version == gap.family_version == 1
    @test !positive.is_complete
    @test gap.is_complete
    @test positive.max_degree == 2
    @test gap.max_degree == 1
    @test positive.site_ids == collect(eachindex(patch.sites))
    @test gap.site_ids == patch.inner_ids
    @test length(positive.entries) == 703
    @test length(gap.entries) == 7
    @test length(unique(positive.entries)) == length(positive.entries)
    @test length(unique(gap.entries)) == length(gap.entries)
    @test all(entry -> state_monomial_degree(entry) <= 2, positive.entries)
    @test all(entry -> state_monomial_degree(entry) <= 1, gap.entries)
    @test all(
        entry -> all(
            factor -> factor[1] in patch.inner_ids,
            Iterators.flatten((
                entry.operator_word.ops,
                (word.ops for word in entry.state_symbols)...,
            )),
        ),
        gap.entries,
    )
    @test state_monomial_string(first(positive.entries)) == "zeta=[];op=I"
    @test issorted(state_monomial_degree.(positive.entries))
    @test issubset(Set(gap.entries), Set(positive.entries))
    @test positive.sha256 == basis_manifest(problem, :positive).sha256
    @test gap.sha256 == basis_manifest(problem, :gap).sha256
    @test positive.sha256 ==
          "83befe24c09bccdc7d228fc60c606d301dd76c10688121e1e466d43a583d5c13"
    @test gap.sha256 ==
          "5be3d2db7be104d1bc431898496e8e34116787a7f14a30886fa6933924bea169"
    @test positive.sha256 != gap.sha256
    @test validate_basis_manifest(positive)
    @test validate_basis_manifest(gap)
    @test validate_basis_manifest(positive, problem, :positive)
    @test validate_basis_manifest(gap, problem, :gap)
    @test !validate_basis_manifest(positive, problem, :gap)
    @test !validate_basis_manifest(gap, problem, :positive)

    forge_manifest = function(role, site_ids, max_degree)
        entries, is_complete, selection_rule =
            GenericGapModel.structured_basis_contents(spec, site_ids, max_degree)
        fingerprint = GenericGapModel.manifest_fingerprint(
            spec,
            role,
            site_ids,
            max_degree,
            is_complete,
            selection_rule,
            entries,
        )
        return BasisManifest(
            role,
            spec.family,
            spec.version,
            site_ids,
            max_degree,
            entries,
            is_complete,
            selection_rule,
            fingerprint,
        )
    end

    role_flipped = forge_manifest(:gap, positive.site_ids, positive.max_degree)
    @test validate_basis_manifest(role_flipped)
    @test !validate_basis_manifest(role_flipped, problem, :gap)

    wrong_gap_sites = [only(gap.site_ids) + 1]
    wrong_site_manifest = forge_manifest(:gap, wrong_gap_sites, gap.max_degree)
    @test validate_basis_manifest(wrong_site_manifest)
    @test !validate_basis_manifest(wrong_site_manifest, problem, :gap)

    wrong_degree_manifest =
        forge_manifest(:gap, gap.site_ids, gap.max_degree + 1)
    @test validate_basis_manifest(wrong_degree_manifest)
    @test !validate_basis_manifest(wrong_degree_manifest, problem, :gap)

    wrong_hash_manifest = BasisManifest(
        gap.role,
        gap.family,
        gap.family_version,
        gap.site_ids,
        gap.max_degree,
        gap.entries,
        gap.is_complete,
        gap.selection_rule,
        repeat("0", 64),
    )
    @test !validate_basis_manifest(wrong_hash_manifest)
    @test !validate_basis_manifest(wrong_hash_manifest, problem, :gap)

    input_state_word = PauliWord([(1, UInt8(1))])
    input_operator_word = PauliWord([(2, UInt8(2))])
    owned_monomial = StateMonomial([input_state_word], input_operator_word)
    push!(input_state_word.ops, (3, UInt8(3)))
    push!(input_operator_word.ops, (4, UInt8(1)))
    @test state_monomial_string(owned_monomial) ==
          "zeta=[1X];op=2Y"

    tampered = deepcopy(positive)
    push!(tampered.entries, first(tampered.entries))
    @test !validate_basis_manifest(tampered)
    nested_tamper = deepcopy(gap)
    push!(nested_tamper.entries[2].operator_word.ops, (99, UInt8(1)))
    @test !validate_basis_manifest(nested_tamper)

    constructor_sites = copy(gap.site_ids)
    constructor_entries = deepcopy(gap.entries)
    owned_manifest = BasisManifest(
        gap.role,
        gap.family,
        gap.family_version,
        constructor_sites,
        gap.max_degree,
        constructor_entries,
        gap.is_complete,
        gap.selection_rule,
        gap.sha256,
    )
    push!(constructor_sites, 99)
    push!(constructor_entries, first(constructor_entries))
    push!(constructor_entries[2].operator_word.ops, (99, UInt8(1)))
    @test owned_manifest.site_ids == gap.site_ids
    @test validate_basis_manifest(owned_manifest)

    truncated_entries = positive.entries[1:1]
    truncated_sha = GenericGapModel.manifest_fingerprint(
        spec,
        positive.role,
        positive.site_ids,
        positive.max_degree,
        positive.is_complete,
        positive.selection_rule,
        truncated_entries,
    )
    truncated = BasisManifest(
        positive.role,
        positive.family,
        positive.family_version,
        positive.site_ids,
        positive.max_degree,
        truncated_entries,
        positive.is_complete,
        positive.selection_rule,
        truncated_sha,
    )
    @test !validate_basis_manifest(truncated)
    @test plan.positive_basis_dimension == length(positive.entries)
    @test plan.gap_basis_dimension == length(gap.entries)
    @test !plan.is_complete
    @test plan.positive_basis_sha256 == positive.sha256
    @test plan.gap_basis_sha256 == gap.sha256
    @test plan.problem_sha256 ==
          "f6f7cd7a0cc2e053e40ecd82f52a24438536869e3340b959cd7f68cab4467f4e"

    for nsites in (1, 2, 9), max_degree in 0:2
        site_ids = collect(1:nsites)
        entries, is_complete, _ =
            GenericGapModel.structured_basis_contents(spec, site_ids, max_degree)
        @test BigInt(length(entries)) ==
              one_symbol_lift_count(nsites, max_degree)
        @test is_complete ==
              (max_degree <= 1)
        @test is_complete ==
              (BigInt(length(entries)) ==
               full_state_basis_count(nsites, max_degree))
    end

    higher_problem = GapProblem(
        patch,
        model,
        1//10,
        3;
        basis_mode=:structured,
        basis_spec=spec,
    )
    higher_positive = basis_manifest(higher_problem, :positive)
    @test higher_positive.entries[1:length(positive.entries)] == positive.entries
    higher_gap = basis_manifest(higher_problem, :gap)
    @test higher_gap.entries[1:length(gap.entries)] == gap.entries
    @test !higher_gap.is_complete

    changed_model_problem = GapProblem(
        patch,
        square_j1j2_model(107//200),
        1//5,
        2;
        basis_mode=:structured,
        basis_spec=spec,
    )
    @test basis_manifest(changed_model_problem, :positive).sha256 == positive.sha256
    @test assembly_plan(changed_model_problem).problem_sha256 != plan.problem_sha256

    changed_gamma_problem = GapProblem(
        patch,
        model,
        1//5,
        2;
        basis_mode=:structured,
        basis_spec=spec,
    )
    @test basis_manifest(changed_gamma_problem, :positive).sha256 == positive.sha256
    @test assembly_plan(changed_gamma_problem).problem_sha256 != plan.problem_sha256

    wider_patch = square_patch_geometry(2)
    wider_problem = GapProblem(
        wider_patch,
        model,
        1//10,
        2;
        basis_mode=:structured,
        basis_spec=spec,
    )
    wider_gap = basis_manifest(wider_problem, :gap)
    @test wider_gap.site_ids == wider_patch.inner_ids
    @test wider_gap.site_ids != collect(1:length(wider_gap.site_ids))
    @test length(wider_gap.entries) == 55
    @test validate_basis_manifest(wider_gap)

    permuted_inner_ids = reverse(wider_patch.inner_ids)
    permuted_patch = LocalPatch(
        wider_patch.name,
        wider_patch.level,
        copy(wider_patch.sites),
        copy(wider_patch.site_to_id),
        copy(permuted_inner_ids),
    )
    permuted_problem = GapProblem(
        permuted_patch,
        model,
        1//10,
        2;
        basis_mode=:structured,
        basis_spec=spec,
    )
    permuted_gap = basis_manifest(permuted_problem, :gap)
    @test permuted_gap.site_ids == wider_gap.site_ids
    @test permuted_gap.entries == wider_gap.entries
    @test permuted_gap.sha256 == wider_gap.sha256
    @test permuted_problem.patch.inner_ids == permuted_inner_ids
    @test validate_basis_manifest(permuted_gap, permuted_problem, :gap)
    wider_plan = assembly_plan(wider_problem)
    permuted_plan = assembly_plan(permuted_problem)
    @test permuted_plan.problem_sha256 == wider_plan.problem_sha256
    @test all(
        getproperty(permuted_plan, field) == getproperty(wider_plan, field)
        for field in propertynames(wider_plan)
    )

    duplicate_patch = LocalPatch(
        wider_patch.name,
        wider_patch.level,
        copy(wider_patch.sites),
        copy(wider_patch.site_to_id),
        [wider_patch.inner_ids; first(wider_patch.inner_ids)],
    )
    duplicate_problem = GapProblem(
        duplicate_patch,
        model,
        1//10,
        2;
        basis_mode=:structured,
        basis_spec=spec,
    )
    @test_throws ArgumentError basis_manifest(duplicate_problem, :gap)

    invalid_patch = LocalPatch(
        wider_patch.name,
        wider_patch.level,
        copy(wider_patch.sites),
        copy(wider_patch.site_to_id),
        copy(wider_patch.inner_ids),
    )
    invalid_problem = GapProblem(
        invalid_patch,
        model,
        1//10,
        2;
        basis_mode=:structured,
        basis_spec=spec,
    )
    push!(invalid_problem.patch.inner_ids, length(invalid_patch.sites) + 1)
    @test_throws ArgumentError basis_manifest(invalid_problem, :gap)
end

include(joinpath(@__DIR__, "..", "src", "CoreMGK.jl"))
using .CoreMGK
include(joinpath(@__DIR__, "..", "src", "SharedCoreWire.jl"))
using .SharedCoreWire

@testset "shared core canonical wire grammar" begin
    grammar_hex =
        "4149434f5245310a4f343a53343a666c616742313b53353a6c6162656c" *
        "53343a783a790a53343a6e6f6e654e3b53363a76616c7565734c333a49" *
        "303b492d323b51312f323b0a"
    grammar_bytes = hex2bytes(grammar_hex)
    @test length(grammar_bytes) == 70
    @test bytes2hex(sha256(grammar_bytes)) ==
          "b75e0d4cab1d35247d2f654bd9c65566cac050be7202baadbf82c42120958975"
    grammar_value = decode_math_bytes(grammar_bytes)
    @test canonical_math_bytes(grammar_value) == grammar_bytes
    envelope_bytes = canonical_framed_bytes("AICOREENV1", grammar_value)
    @test canonical_framed_bytes(
        "AICOREENV1",
        decode_framed_bytes(envelope_bytes, "AICOREENV1"),
    ) == envelope_bytes
    @test_throws ArgumentError canonical_framed_bytes("not-a-frame", grammar_value)

    hz_hex =
        "4149434f5245310a4f383a5331323a675f70726f647563745f7879432d31" *
        "2f312c302f313b53363a675f7a5f787943302f312c312f313b53363a6b5f" *
        "695f787943302f312c2d322f313b53363a6d5f7a5f787943302f312c312f" *
        "313b5331373a7061636b5f675f70726f647563745f7265512d322f313b53" *
        "31313a7061636b5f675f7a5f696d51322f313b5331313a7061636b5f6b5f" *
        "695f696d512d342f313b5331313a7061636b5f6d5f7a5f696d51322f313b" *
        "0a"
    hz_bytes = hex2bytes(hz_hex)
    @test length(hz_bytes) == 181
    @test bytes2hex(sha256(hz_bytes)) ==
          "4c8f8281c27eee74f10dac79b0861c256f46bece9b49f6204bbfd01503970ad4"
    @test canonical_math_bytes(decode_math_bytes(hz_bytes)) == hz_bytes

    @test_throws ErrorException decode_math_bytes(
        Vector{UInt8}(codeunits("AICORE1\nO2:S1:bI1;S1:aI2;\n")),
    )
    @test_throws ErrorException decode_math_bytes(
        Vector{UInt8}(codeunits("AICORE1\nQ2/2;\n")),
    )
    @test_throws ArgumentError canonical_math_bytes(Dict("x" => 0.5))
    @test_throws ArgumentError canonical_math_bytes(("not", "a", "list"))

    _, x_word = pauli_word([(2, :X)])
    _, y_word = pauli_word([(1, :Y)])
    monomial = StateMonomial([x_word], y_word)
    scalar_row = ScalarMoment([x_word, y_word])
    @test startswith(term_id(x_word), "h:")
    @test startswith(entry_id(monomial), "be:")
    @test startswith(row_id(scalar_row), "row:")
    @test term_id(x_word) != term_id(x_word; site_namespace="legacy")
    reversed_row = ScalarMoment([y_word, x_word])
    @test row_id(scalar_row) == row_id(reversed_row)
end

include(joinpath(@__DIR__, "..", "src", "SquareGapConic.jl"))
using .SquareGapConic
include(joinpath(@__DIR__, "..", "src", "SquareCoreInventory.jl"))
using .SquareCoreInventory
import JuMP

@testset "complex Hermitian to real PSD rendering" begin
    one = BigInt(1) // BigInt(1)
    two = BigInt(2) // BigInt(1)
    @test realify_hermitian_coefficient(
        GaussianRational(one, two),
        1,
        2,
        2,
    ) == [
        (1, 2, one),
        (3, 4, one),
        (1, 4, -two),
        (2, 3, two),
    ]
    @test realify_hermitian_coefficient(
        GaussianRational(one, zero(one)),
        1,
        1,
        2,
    ) == [(1, 1, one), (3, 3, one)]
    @test isempty(realify_hermitian_coefficient(
        zero(GaussianRational),
        1,
        2,
        2,
    ))
    @test_throws ArgumentError realify_hermitian_coefficient(
        GaussianRational(one, one),
        1,
        1,
        2,
    )
    @test_throws ArgumentError realify_hermitian_coefficient(
        GaussianRational(one, zero(one)),
        2,
        1,
        2,
    )
end

@testset "Square shared-core inventory declarations" begin
    problem = GapProblem(
        square_patch_geometry(1),
        square_j1j2_model(1 // 2),
        1 // 10,
        2;
        basis_mode=:structured,
        basis_spec=StructuredBasisSpec(:one_symbol_lift, 1),
    )
    source = core_mgk_plan(problem)
    positive_selector =
        SquareCoreInventory.selector_spec(source.positive_basis)
    gap_selector = SquareCoreInventory.selector_spec(source.gap_basis)
    @test SquareCoreInventory.field(
        positive_selector,
        "basis_selector_sha256",
    ) != SquareCoreInventory.field(gap_selector, "basis_selector_sha256")
    @test SquareCoreInventory.field(
        positive_selector,
        "basis_selector_id",
    ) == "structured-one-symbol-lift"
    action, reference = SquareCoreInventory.unrestricted_action(
        [source.positive_basis, source.gap_basis],
    )
    @test SquareCoreInventory.field(action, "mode") == "unrestricted"
    @test SquareCoreInventory.field(reference, "action_sha256") ==
          SquareCoreInventory.field(action, "action_sha256")
    config = SquareCoreInventory.exact_model_config(problem)
    parameters = SquareCoreInventory.field(config, "parameters")
    @test SquareCoreInventory.field(parameters, "j1") == 1 // 1
    @test SquareCoreInventory.field(parameters, "j2_over_j1") == 1 // 2
    @test SquareCoreInventory.field(parameters, "spin_operator_scale") == 1 // 2
    inventory_source = read(
        joinpath(@__DIR__, "..", "src", "SquareCoreInventory.jl"),
        String,
    )
    @test !occursin("optimize!", inventory_source)
    @test !occursin("Mosek", inventory_source)
end

@testset "solver-free Square conic render" begin
    problem = GapProblem(
        square_patch_geometry(1),
        square_j1j2_model(1 // 2),
        1 // 10,
        2;
        basis_mode=:structured,
        basis_spec=StructuredBasisSpec(:one_symbol_lift, 1),
    )
    source = core_mgk_plan(problem)
    one = BigInt(1) // BigInt(1)
    identity = ScalarMoment(PauliWord[])
    plan = SquareConicPlan(
        source,
        BigInt(1) // BigInt(10),
        [identity],
        [row_id(identity)],
        ExactAffineConstraint(:normalization, "L(1)=1", [1 => one], one),
        ExactAffineConstraint[],
        7,
        0,
        [
            RealPSDPlan(:positive, 1, 2, [RealPSDTerm(1, 1, one)]),
            RealPSDPlan(:gap, 1, 2, [RealPSDTerm(1, 1, one)]),
        ],
        :feasibility,
        Pair{Int,Rational{BigInt}}[],
    )
    @test plan.gamma == BigInt(1) // BigInt(10)
    @test plan.normalization.rhs == BigInt(1) // BigInt(1)
    @test plan.stationarity_selector_entries == 7
    @test isempty(plan.stationarity)
    @test plan.stationarity_exact_duplicates_removed == 0
    @test plan.objective_sense == :feasibility
    @test isempty(plan.objective_terms)
    @test [block.complex_dimension for block in plan.psd_blocks] == [1, 1]
    @test [block.real_dimension for block in plan.psd_blocks] == [2, 2]

    mktempdir() do directory
        path = joinpath(directory, "tiny-square.mof.json")
        rendered = render_mof(plan, path)
        @test rendered.optimizer_invoked == false
        model = JuMP.MOI.FileFormats.Model(filename=path)
        JuMP.MOI.read_from_file(model, path)
        @test JuMP.MOI.get(model, JuMP.MOI.ObjectiveSense()) ==
              JuMP.MOI.FEASIBILITY_SENSE
        @test length(JuMP.MOI.get(
            model,
            JuMP.MOI.ListOfVariableIndices(),
        )) == length(plan.rows)
        equalities = JuMP.MOI.get(
            model,
            JuMP.MOI.ListOfConstraintIndices{
                JuMP.MOI.ScalarAffineFunction{Float64},
                JuMP.MOI.EqualTo{Float64},
            }(),
        )
        @test length(equalities) == 1
        @test JuMP.MOI.get(
            model,
            JuMP.MOI.ConstraintSet(),
            only(equalities),
        ).value == 1.0
        psd_constraints = JuMP.MOI.get(
            model,
            JuMP.MOI.ListOfConstraintIndices{
                JuMP.MOI.VectorAffineFunction{Float64},
                JuMP.MOI.PositiveSemidefiniteConeTriangle,
            }(),
        )
        @test [
            JuMP.MOI.get(model, JuMP.MOI.ConstraintSet(), constraint).side_dimension
            for constraint in psd_constraints
        ] == [2, 2]
        audit = audit_rendered_mof(plan, path)
        @test audit.variables == 1
        @test audit.affine_equalities == 1
        @test audit.psd_dimensions == [2, 2]
        @test audit.exact_coefficient_match
        @test audit.optimizer_invoked == false
    end
end

function exact_component_map(components, name)
    component = only(filter(record -> record.component == name, components))
    return Dict(
        coefficient.row => coefficient.coefficient
        for coefficient in component.coefficients
    )
end

function conjugate_maps_equal(left, right)
    return keys(left) == keys(right) &&
           all(right[row] == conj(value) for (row, value) in left)
end

@testset "exact core M/G/K pair algebra" begin
    _, x = pauli_word([(1, :X)])
    _, y = pauli_word([(1, :Y)])
    _, z = pauli_word([(1, :Z)])
    identity = PauliWord()
    x_entry = StateMonomial(PauliWord[], x)
    y_entry = StateMonomial(PauliWord[], y)
    h_z = LocalPauliTerm[
        LocalPauliTerm(Complex(1//1, 0//1), z, :hand_golden, Site(0, 0)),
    ]
    row_i = ScalarMoment(PauliWord[])
    row_x = ScalarMoment([x])
    row_y = ScalarMoment([y])
    row_z = ScalarMoment([z])
    row_xy = ScalarMoment([x, y])

    m_xy = exact_component_map(positive_pair_components(x_entry, y_entry), :M)
    @test m_xy == Dict(row_z => GaussianRational(0, 1))

    gap_xy = gap_pair_components(h_z, x_entry, y_entry)
    @test exact_component_map(gap_xy, :K) ==
          Dict(row_i => GaussianRational(0, -2))
    @test exact_component_map(gap_xy, :G_moment) ==
          Dict(row_z => GaussianRational(0, 1))
    @test exact_component_map(gap_xy, :G_product) ==
          Dict(row_xy => GaussianRational(-1, 0))

    gap_yx = gap_pair_components(h_z, y_entry, x_entry)
    for component in (:K, :G_moment, :G_product)
        @test conjugate_maps_equal(
            exact_component_map(gap_xy, component),
            exact_component_map(gap_yx, component),
        )
    end
    @test conjugate_maps_equal(
        m_xy,
        exact_component_map(positive_pair_components(y_entry, x_entry), :M),
    )

    gap_xx = gap_pair_components(h_z, x_entry, x_entry)
    @test exact_component_map(gap_xx, :K) ==
          Dict(row_z => GaussianRational(-2, 0))
    @test exact_component_map(gap_xx, :G_moment) ==
          Dict(row_i => GaussianRational(1, 0))
    @test exact_component_map(gap_xx, :G_product) ==
          Dict(ScalarMoment([x, x]) => GaussianRational(-1, 0))

    wiring_xy = ExactPairWiring(:gap, 1, 2, x_entry, y_entry, gap_xy)
    a_gamma = a_gamma_coefficients(wiring_xy)
    @test a_gamma[row_i] == GammaAffineCoefficient(
        GaussianRational(0, -2),
        GaussianRational(0, 0),
    )
    @test a_gamma[row_z] == GammaAffineCoefficient(
        GaussianRational(0, 0),
        GaussianRational(0, -1),
    )
    @test a_gamma[row_xy] == GammaAffineCoefficient(
        GaussianRational(0, 0),
        GaussianRational(1, 0),
    )

    @test pack_upper_coefficient(GaussianRational(0, 1), 1, 2) ==
          (real=0//1, imag=2//1)
    @test pack_upper_coefficient(GaussianRational(0, -2), 1, 2) ==
          (real=0//1, imag=-4//1)
    @test pack_upper_coefficient(GaussianRational(-1, 0), 1, 2) ==
          (real=-2//1, imag=0//1)
    @test pack_upper_coefficient(GaussianRational(-2, 0), 1, 1) ==
          (real=-2//1, imag=0//1)
    @test_throws ArgumentError pack_upper_coefficient(
        GaussianRational(0, 1),
        1,
        1,
    )

    scalar_x = StateMonomial([x], identity)
    zero_k = only(filter(
        record -> record.component == :K,
        gap_pair_components(h_z, scalar_x, scalar_x),
    ))
    @test zero_k.status == :computed_exact_zero
    @test zero_k.zero_reason == :algebraic
    @test isempty(zero_k.coefficients)

    inexact_h = LocalPauliTerm[
        LocalPauliTerm(ComplexF64(1.0, 0.0), z, :inexact, Site(0, 0)),
    ]
    @test_throws ArgumentError gap_pair_components(inexact_h, x_entry, y_entry)
end

@testset "Square J1-J2 core M/G/K source gate" begin
    patch = square_patch_geometry(1)
    model = square_j1j2_model(1//2)
    spec = StructuredBasisSpec(:one_symbol_lift, 1)
    problem = GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:structured,
        basis_spec=spec,
    )
    plan = core_mgk_plan(problem)
    @test plan.state_class ==
          "unrestricted infinite-volume KMS ground states; flat basis; no symmetry quotient"
    @test length(plan.hamiltonian_terms) == 60
    @test length(plan.positive_basis.entries) == 703
    @test length(plan.gap_basis.entries) == 7
    @test plan.positive_basis.sha256 ==
          "83befe24c09bccdc7d228fc60c606d301dd76c10688121e1e466d43a583d5c13"
    @test plan.gap_basis.sha256 ==
          "5be3d2db7be104d1bc431898496e8e34116787a7f14a30886fa6933924bea169"

    identity_wiring = core_mgk_pair(plan, :positive, 1, 1)
    @test exact_component_map(identity_wiring.component_records, :M) ==
          Dict(ScalarMoment(PauliWord[]) => GaussianRational(1, 0))
    for j in 1:length(plan.gap_basis.entries), k in j:length(plan.gap_basis.entries)
        wiring = core_mgk_pair(plan, :gap, j, k)
        @test getproperty.(wiring.component_records, :component) ==
              [:K, :G_moment, :G_product]
        reverse_components = gap_pair_components(
            plan.hamiltonian_terms,
            plan.gap_basis.entries[k],
            plan.gap_basis.entries[j],
        )
        for component in (:K, :G_moment, :G_product)
            @test conjugate_maps_equal(
                exact_component_map(wiring.component_records, component),
                exact_component_map(reverse_components, component),
            )
        end
    end

    for (j, k) in ((1, 2), (2, 7), (100, 353), (353, 703))
        forward = positive_pair_components(
            plan.positive_basis.entries[j],
            plan.positive_basis.entries[k],
        )
        reverse = positive_pair_components(
            plan.positive_basis.entries[k],
            plan.positive_basis.entries[j],
        )
        @test conjugate_maps_equal(
            exact_component_map(forward, :M),
            exact_component_map(reverse, :M),
        )
    end

    @test_throws ArgumentError core_mgk_pair(plan, :gap, 2, 1)
    restricted = GapProblem(
        patch,
        model,
        1//10,
        2;
        basis_mode=:structured,
        basis_spec=spec,
        symmetry=ExplicitStateSymmetry("D4", ["C4", "mirror"]),
    )
    @test_throws ArgumentError core_mgk_plan(restricted)
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
    tfim_scan = parse_point(
        "tfim",
        "0.258";
        allow_experimental_gamma=true,
    )
    kagome_scan = parse_point(
        "kagome",
        "1.275";
        allow_experimental_gamma=true,
    )
    @test point_id(tfim_scan) == "tfim-n9-g0p5-d2-lso6-gamma0p258"
    @test point_id(kagome_scan) == "kagome-n13-d3-lso5-gamma1p275"
    @test_throws ArgumentError parse_point(
        "tfim",
        "0";
        allow_experimental_gamma=true,
    )

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
          "99d7938bd8758251042146de518b0a8d4234cebfee8a5ee2856bda569101acf3"
    @test basis_metadata(kagome).fingerprint ==
          "98b97cf96837c4b21ee24f6326eaf1d7bc1626a2cc880c5ccc131e327e56ff54"

    named = adapt_certify_result(
        (
            flag=0,
            termination=:SLOW_PROGRESS,
            primal=:INFEASIBILITY_CERTIFICATE,
            dual=:NO_SOLUTION,
            objective=1.25,
            audit=(requested=false,),
        ),
    )
    @test named.adapter == "namedtuple-v1"
    @test isequal(named.flag, 0)
    @test named.termination == "SLOW_PROGRESS"
    @test named.primal == "INFEASIBILITY_CERTIFICATE"
    @test named.dual == "NO_SOLUTION"
    @test named.objective == 1.25
    @test named.objective_availability == "available"
    @test named.solver_audit == (requested=false,)

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
