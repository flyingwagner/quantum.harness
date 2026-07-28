import JuMP

const RayMOI = JuMP.MOI

function write_synthetic_ray_model(path::AbstractString; sense=RayMOI.MAX_SENSE)
    model = RayMOI.Utilities.Model{Float64}()
    variables = RayMOI.add_variables(model, 4)
    for (ordinal, variable) in enumerate(variables)
        RayMOI.set(model, RayMOI.VariableName(), variable, "C$ordinal")
    end

    # The constants are deliberately nonzero. A homogeneous recession ray must
    # satisfy only the direction x1-x2=0, not the affine equation at the origin.
    equality = RayMOI.ScalarAffineFunction(
        [
            RayMOI.ScalarAffineTerm(1.0, variables[1]),
            RayMOI.ScalarAffineTerm(-1.0, variables[2]),
        ],
        13.0,
    )
    RayMOI.add_constraint(model, equality, RayMOI.EqualTo(13.0))

    # Triangle order is (1,1), (1,2), (2,2). The constant matrix is irrelevant
    # to the recession cone; the ray matrix is [x1 x3; x3 x2].
    psd = RayMOI.VectorAffineFunction(
        [
            RayMOI.VectorAffineTerm(
                1,
                RayMOI.ScalarAffineTerm(1.0, variables[1]),
            ),
            RayMOI.VectorAffineTerm(
                2,
                RayMOI.ScalarAffineTerm(1.0, variables[3]),
            ),
            RayMOI.VectorAffineTerm(
                3,
                RayMOI.ScalarAffineTerm(1.0, variables[2]),
            ),
        ],
        [4.0, -2.0, 9.0],
    )
    RayMOI.add_constraint(
        model,
        psd,
        RayMOI.PositiveSemidefiniteConeTriangle(2),
    )

    objective = RayMOI.ScalarAffineFunction(
        [RayMOI.ScalarAffineTerm(1.0, variables[4])],
        17.0,
    )
    RayMOI.set(model, RayMOI.ObjectiveSense(), sense)
    RayMOI.set(
        model,
        RayMOI.ObjectiveFunction{typeof(objective)}(),
        objective,
    )
    RayMOI.write_to_file(model, path)
    return variables
end

function write_synthetic_ray(path::AbstractString, values)
    open(path, "w") do io
        println(io, "ordinal\tmoi_index\tname\tvalue")
        for (ordinal, value) in enumerate(values)
            println(io, ordinal, '\t', ordinal, "\t\t", value)
        end
    end
    return path
end

function write_duplicate_equality_model(path::AbstractString)
    model = RayMOI.Utilities.Model{Float64}()
    variables = RayMOI.add_variables(model, 2)
    function vector_term(output, coefficient, variable)
        return RayMOI.VectorAffineTerm(
            output,
            RayMOI.ScalarAffineTerm(coefficient, variable),
        )
    end
    duplicate_vector = RayMOI.VectorAffineFunction(
        [
            vector_term(1, 1.0, variables[1]),
            vector_term(1, -1.0, variables[2]),
            vector_term(2, 1.0, variables[1]),
            vector_term(2, -1.0, variables[2]),
            vector_term(3, 1.0, variables[1]),
            vector_term(3, 1.0, variables[2]),
        ],
        [3.0, 3.0, 4.0],
    )
    RayMOI.add_constraint(model, duplicate_vector, RayMOI.Zeros(3))
    RayMOI.set(model, RayMOI.ObjectiveSense(), RayMOI.MAX_SENSE)
    RayMOI.set(
        model,
        RayMOI.ObjectiveFunction{RayMOI.VariableIndex}(),
        variables[1],
    )
    RayMOI.write_to_file(model, path)
end

@testset "solver-free homogeneous conic-ray verifier" begin
    mktempdir() do directory
        max_model = joinpath(directory, "max.mof.json")
        min_model = joinpath(directory, "min.mof.json")
        write_synthetic_ray_model(max_model)
        write_synthetic_ray_model(min_model; sense=RayMOI.MIN_SENSE)

        function replay(values; model=max_model)
            ray_path = joinpath(directory, "ray.tsv")
            write_synthetic_ray(ray_path, values)
            return GapRayVerifier.verify(model, ray_path)
        end

        accepted = replay([1.0, 1.0, 0.0, 1.0])
        @test accepted.verified
        @test accepted.verdict == "accepted_floating_point_ray"
        @test accepted.rigor == "floating_point_replay"
        @test accepted.equality_pass
        @test accepted.cone_pass
        @test accepted.objective_pass
        @test isempty(accepted.rejection_reasons)

        equality_bad = replay([1.0, 0.0, 0.0, 1.0])
        @test !equality_bad.verified
        @test !equality_bad.equality_pass
        @test equality_bad.cone_pass
        @test equality_bad.objective_pass
        @test equality_bad.rejection_reasons == ("equality_residual",)

        psd_bad = replay([1.0, 1.0, 2.0, 1.0])
        @test !psd_bad.verified
        @test psd_bad.equality_pass
        @test !psd_bad.cone_pass
        @test psd_bad.objective_pass
        @test psd_bad.rejection_reasons == ("cone_violation",)

        wrong_max_sign = replay([1.0, 1.0, 0.0, -1.0])
        @test !wrong_max_sign.verified
        @test wrong_max_sign.equality_pass
        @test wrong_max_sign.cone_pass
        @test !wrong_max_sign.objective_pass
        @test wrong_max_sign.rejection_reasons == ("objective_direction",)

        accepted_min = replay([1.0, 1.0, 0.0, -1.0]; model=min_model)
        rejected_min = replay([1.0, 1.0, 0.0, 1.0]; model=min_model)
        @test accepted_min.verified
        @test !rejected_min.verified
        @test rejected_min.rejection_reasons == ("objective_direction",)

        # Scaling a valid ray changes neither homogeneous feasibility nor the
        # normalized verdict. Nonzero affine constants above are subtracted.
        scaled = replay([1.0e16, 1.0e16, 0.0, 1.0e12])
        @test scaled.verified
        @test scaled.max_equality_residual_relative == 0.0
        @test scaled.improving_objective_relative == 1.0e-4

        # Kagome-style pathology: the PSD and objective checks look excellent
        # after normalization, but a 1e-10 relative equality defect must reject.
        pathological = replay([1.0e16, 1.0e16 - 1.0e6, 0.0, 1.0e12])
        @test !pathological.verified
        @test !pathological.equality_pass
        @test pathological.cone_pass
        @test pathological.objective_pass
        @test pathological.max_equality_residual_relative > 9.9e-11
        @test pathological.max_equality_residual_relative < 1.01e-10
        @test pathological.rejection_reasons == ("equality_residual",)

        conditioning = equality_conditioning(
            extract_exact_problem(max_model),
            [1.0e16, 1.0e16 - 1.0e6, 0.0, 1.0e12],
        )
        @test conditioning[1].residual_high_precision == 1.0e6
        @test conditioning[1].backward_error > 4.9e-11
        @test conditioning[1].backward_error < 5.1e-11
        @test conditioning[1].row_scaled_residual == 1.0e-10
        @test conditioning[1].summation_difference == 0.0

        block_scales = psd_block_scales(
            extract_exact_problem(max_model),
            [1.0e16, 1.0e16 - 1.0e6, 0.0, 1.0e12],
        )
        @test only(block_scales).dimension == 2
        @test only(block_scales).maximum == 1.0e16
        @test only(block_scales).minimum_nonzero == 1.0e16 - 1.0e6

        @test_throws ArgumentError GapRayVerifier.verify(
            max_model,
            write_synthetic_ray(joinpath(directory, "bad-tol.tsv"), [1, 1, 0, 1]);
            relative_tolerance=-1,
        )

        exact_problem = extract_exact_problem(max_model; coefficient_tolerance=1e-14)
        @test exact_problem.equality_offsets == BigRational[0]
        exact_candidate, _, _ = normalize_rational_ray(
            [1.0, 1.0, 0.0, 1.0];
            rational_tolerance=1e-12,
        )
        exact_corrected, exact_unresolved, _ = correct_with_private_pivots(
            exact_problem,
            exact_candidate,
        )
        @test isempty(exact_unresolved)
        @test all(iszero, exact_residuals(exact_problem, exact_corrected))
        @test objective_improvement(exact_problem, exact_corrected) == 1
        @test rigorous_psd_proof(
            only(exact_problem.psd_blocks),
            exact_corrected,
        ).proved

        indefinite_candidate = BigRational[1, 1, 2, 1]
        @test !rigorous_psd_proof(
            only(exact_problem.psd_blocks),
            indefinite_candidate,
        ).proved

        zero_row_block = PSDDirectionBlock(
            2,
            [
                ExactAffineRow([1 => BigRational(1)]),
                ExactAffineRow([2 => BigRational(1)]),
                ExactAffineRow([3 => BigRational(1)]),
            ],
        )
        zero_row_proof = rigorous_psd_proof(
            zero_row_block,
            BigRational[1, 0, 0],
        )
        @test zero_row_proof.proved
        @test zero_row_proof.exact_zero_rows == (2,)

        exact_path = joinpath(directory, "exact-ray.tsv")
        write_exact_ray(exact_path, exact_corrected)
        @test readlines(exact_path)[1] == "ordinal\tnumerator\tdenominator"

        duplicate_path = joinpath(directory, "duplicate.mof.json")
        deduplicated_path = joinpath(directory, "deduplicated.mof.json")
        write_duplicate_equality_model(duplicate_path)
        duplicate_model = RayMOI.FileFormats.Model(filename=duplicate_path)
        RayMOI.read_from_file(duplicate_model, duplicate_path)
        deduplication = deduplicate_affine_equalities!(duplicate_model)
        @test deduplication.scalar_equalities_removed == 0
        @test deduplication.vector_zero_coordinates_removed == 1
        @test deduplication.total_removed == 1
        @test length(RayMOI.get(duplicate_model, RayMOI.ListOfVariableIndices())) == 2
        RayMOI.write_to_file(duplicate_model, deduplicated_path)
        deduplicated = extract_exact_problem(deduplicated_path)
        @test length(deduplicated.equalities) == 2
        @test Set(
            (Tuple(row.terms), offset) for
            (row, offset) in zip(
                deduplicated.equalities,
                deduplicated.equality_offsets,
            )
        ) == Set([
            ((1 => BigRational(1), 2 => BigRational(-1)), BigRational(3)),
            ((1 => BigRational(1), 2 => BigRational(1)), BigRational(4)),
        ])
    end
end
