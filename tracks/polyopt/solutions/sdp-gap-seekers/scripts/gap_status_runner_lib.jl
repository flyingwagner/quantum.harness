module GapStatusRunner

using Dates
using SHA

export RunPoint,
       adapt_certify_result,
       allowed_points,
       basis_metadata,
       build_dry_run_record,
       build_hamiltonian_data,
       execute_point,
       hamiltonian_metadata,
       json_encode,
       main,
       parse_point,
       parse_run_spec_point,
       point_id

const SCHEMA_VERSION = "gap-status-result-v1"
const EXPECTED_HARNESS_BASE = "0d2d21bbbb690a49f5e3c3d8142e6f3cdb0c76f7"
const EXPECTED_HARNESS_TREE = "840df11e4507de9a77e78acf22d4334b39759de9"
const EXPECTED_PATCH_SHA256 =
    "5ef9585c71b84b7a07b36610e2bc8aab060a40a8b5062633b070c92dc74fc947"
const EXPECTED_SPECTRALGAP_COMMIT =
    "a1171c906ff2cc2901e58c2426397a2f68c32bb7"
const EXPECTED_SPECTRALGAP_TREE =
    "52d2b037c2275866d482a7dc531198d412c566e1"
const EXPECTED_SPECTRALGAP_FILES = Dict(
    "SpectralGap.jl" =>
        "940cd72b9c4bea39b6daaefd2b9797c54df450cb0afbfb4d9322b2df1b3838bb",
    "basicfunction.jl" =>
        "2095cf7401355f37e9d17915b3ab29d44712d8e40f750eb8449f8c294229b03a",
    "sdp.jl" =>
        "b1fa2280cca51fca38154daf5c767f7538ab68c2297e673eef474da3505f0ccc",
    "strengthening.jl" =>
        "de56b12b17049f81f689d4caef193b9dfd3bf50061fc78b1bf1547a748f7c57b",
)
const RELEVANT_HARNESS_PATHS = (
    "julia-env/Project.toml",
    "tracks/polyopt/solutions/sdp-gap-seekers/spectralgap_a1171c9.patch",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/gap_status_runner_lib.jl",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/gap_status_runner.jl",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/gap_status_array.sbatch",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/prepare_gap_status_run.sh",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/finalize_gap_status_run.sh",
    "skills/using-slurm/profiles/scnet-xh5.toml",
)

const TFIM_GAMMAS = (1 // 4, 13 // 50)
const KAGOME_GAMMAS = (1 // 1, 6 // 5, 63 // 50, 32 // 25)
const KAGOME_TRIPLES = (
    (1, 2, 3),
    (1, 4, 5),
    (2, 6, 7),
    (3, 8, 9),
    (4, 10, 11),
    (5, 12, 13),
)
const KAGOME_INNER_TRIPLES = (KAGOME_TRIPLES[1], KAGOME_TRIPLES[2])

big_rational(value::Rational) =
    BigInt(numerator(value)) // BigInt(denominator(value))

struct RunPoint
    model::Symbol
    N::Int
    g::Union{Nothing, Rational{BigInt}}
    d::Int
    lso::Int
    gamma::Rational{BigInt}
end

function decimal_rational(text::AbstractString)
    occursin(r"^[+-]?[0-9]+(?:\.[0-9]+)?$", text) ||
        throw(ArgumentError("expected an ordinary decimal, got '$text'"))
    sign = startswith(text, "-") ? -1 : 1
    body = startswith(text, ('+', '-')) ? text[2:end] : text
    pieces = split(body, "."; limit=2)
    whole = parse(BigInt, pieces[1])
    if length(pieces) == 1
        return BigInt(sign) * whole // BigInt(1)
    end
    fraction = pieces[2]
    denominator = BigInt(10)^length(fraction)
    numerator = whole * denominator + parse(BigInt, fraction)
    return BigInt(sign) * numerator // denominator
end

function canonical_decimal(value::Rational)
    denominator(value) == 1 && return string(numerator(value))
    for digits in 1:18
        scale = BigInt(10)^digits
        if scale % denominator(value) == 0
            scaled = numerator(value) * (scale ÷ denominator(value))
            sign = scaled < 0 ? "-" : ""
            body = lpad(string(abs(scaled)), digits + 1, '0')
            return string(
                sign,
                body[1:(end - digits)],
                ".",
                body[(end - digits + 1):end],
            )
        end
    end
    return string(numerator(value), "/", denominator(value))
end

function allowed_points()
    points = RunPoint[]
    append!(
        points,
        [
            RunPoint(:tfim, 9, BigInt(1) // BigInt(2), 2, 6, big_rational(gamma))
            for gamma in TFIM_GAMMAS
        ],
    )
    append!(
        points,
        [
            RunPoint(:kagome, 13, nothing, 3, 5, big_rational(gamma))
            for gamma in KAGOME_GAMMAS
        ],
    )
    return points
end

function point_id(point::RunPoint)
    gamma = replace(canonical_decimal(point.gamma), "." => "p")
    return point.model == :tfim ? "tfim-n9-g0p5-d2-lso6-gamma$gamma" :
           "kagome-n13-d3-lso5-gamma$gamma"
end

function run_spec_cell_id(index::Int)
    ids = (
        "01-tfim-gamma-0p25",
        "02-tfim-gamma-0p26",
        "03-kagome-gamma-1",
        "04-kagome-gamma-1p2",
        "05-kagome-gamma-1p26",
        "06-kagome-gamma-1p28",
    )
    1 <= index <= length(ids) || throw(BoundsError(ids, index))
    return ids[index]
end

function locked_run_spec(run_id::AbstractString)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$", run_id) ||
        throw(ArgumentError("invalid run_id in run spec"))
    points = allowed_points()
    return (
        schema_version="gap-status-run-spec-v1",
        run_id=String(run_id),
        cells=[
            (
                cell_id=run_spec_cell_id(index),
                params=(
                    model=string(point.model),
                    gamma=canonical_decimal(point.gamma),
                ),
            ) for (index, point) in enumerate(points)
        ],
    )
end

function parse_run_spec_point(path::AbstractString, cell_index::Int)
    1 <= cell_index <= length(allowed_points()) ||
        throw(ArgumentError("cell index must be between 1 and 6"))
    isfile(path) || throw(ArgumentError("run spec does not exist: $path"))
    text = read(path, String)
    match_result =
        match(r"\"run_id\"\s*:\s*\"([A-Za-z0-9][A-Za-z0-9._-]{0,79})\"", text)
    isnothing(match_result) &&
        throw(ArgumentError("run spec has no valid run_id"))
    run_id = match_result.captures[1]
    normalized_actual = replace(text, r"[ \t\r\n]+" => "")
    normalized_expected = json_encode(locked_run_spec(run_id))
    normalized_actual == normalized_expected ||
        throw(ArgumentError("run spec differs from the locked six-cell contract"))
    return allowed_points()[cell_index], run_id
end

function parse_point(
    model_text::AbstractString,
    gamma_text::AbstractString;
    N_text::Union{Nothing, AbstractString}=nothing,
    g_text::Union{Nothing, AbstractString}=nothing,
    d_text::Union{Nothing, AbstractString}=nothing,
    lso_text::Union{Nothing, AbstractString}=nothing,
)
    model = Symbol(lowercase(model_text))
    model in (:tfim, :kagome) ||
        throw(ArgumentError("model must be exactly 'tfim' or 'kagome'"))
    gamma = decimal_rational(gamma_text)
    if model == :tfim
        N = isnothing(N_text) ? 9 : parse(Int, N_text)
        g = isnothing(g_text) ? BigInt(1) // BigInt(2) : decimal_rational(g_text)
        d = isnothing(d_text) ? 2 : parse(Int, d_text)
        lso = isnothing(lso_text) ? 6 : parse(Int, lso_text)
        N == 9 || throw(ArgumentError("TFIM runner is locked to N=9"))
        g == 1 // 2 || throw(ArgumentError("TFIM runner is locked to g=0.5"))
        d == 2 || throw(ArgumentError("TFIM runner is locked to d=2"))
        lso == 6 || throw(ArgumentError("TFIM runner is locked to lso=6"))
        gamma in TFIM_GAMMAS ||
            throw(ArgumentError("TFIM gamma must be one of 0.25 or 0.26"))
        return RunPoint(:tfim, N, g, d, lso, gamma)
    end

    isnothing(g_text) ||
        throw(ArgumentError("--g is not a parameter of the Kagome point"))
    N = isnothing(N_text) ? 13 : parse(Int, N_text)
    d = isnothing(d_text) ? 3 : parse(Int, d_text)
    lso = isnothing(lso_text) ? 5 : parse(Int, lso_text)
    N == 13 || throw(ArgumentError("Kagome runner is locked to N=13"))
    d == 3 || throw(ArgumentError("Kagome runner is locked to d=3"))
    lso == 5 || throw(ArgumentError("Kagome runner is locked to lso=5"))
    gamma in KAGOME_GAMMAS ||
        throw(ArgumentError("Kagome gamma must be one of 1, 1.2, 1.26, or 1.28"))
    return RunPoint(:kagome, N, nothing, d, lso, gamma)
end

function build_hamiltonian_data(point::RunPoint)
    if point.model == :tfim
        supports = [[3i, 3(i + 1)] for i in 1:(point.N - 1)]
        append!(supports, [[3i - 2] for i in 1:point.N])
        coefficients = Rational{BigInt}[
            fill(BigInt(-1) // BigInt(1), point.N - 1)
            fill(something(point.g), point.N)
        ]
        return supports, coefficients
    end

    supports = Vector{Vector{Int}}()
    for (a, b, c) in KAGOME_TRIPLES
        for (i, j) in ((a, b), (a, c), (b, c))
            push!(supports, [3i - 2, 3j - 2])
            push!(supports, [3i - 1, 3j - 1])
            push!(supports, [3i, 3j])
        end
    end
    coefficients =
        fill(BigInt(1) // BigInt(4), length(supports))
    return supports, coefficients
end

sha256_hex(data::AbstractString) = bytes2hex(sha256(codeunits(data)))
file_sha256(path::AbstractString) = open(path, "r") do io
    bytes2hex(sha256(io))
end

function hamiltonian_metadata(point::RunPoint)
    supports, coefficients = build_hamiltonian_data(point)
    rows = String["hamiltonian-support-coefficient-v1", "model=$(point.model)"]
    for (index, (support, coefficient)) in
        enumerate(zip(supports, coefficients))
        push!(
            rows,
            string(
                index,
                "|",
                numerator(coefficient),
                "/",
                denominator(coefficient),
                "|",
                join(support, ","),
            ),
        )
    end
    convention =
        point.model == :tfim ?
        "H=-sum_{i=1}^8 Z_i Z_{i+1}+0.5 sum_{i=1}^9 X_i; open boundary" :
        "H=sum_{six triangles} sum_{three pairs} 0.25(XX+YY+ZZ)"
    return (
        fingerprint=sha256_hex(join(rows, "\n") * "\n"),
        fingerprint_scheme="hamiltonian-support-coefficient-v1",
        term_count=length(supports),
        convention=convention,
    )
end

function basis_metadata(point::RunPoint)
    positive_blocks =
        point.model == :tfim ? [211, 50] : [271, 104]
    gap_blocks = point.model == :tfim ? [11, 14] : [18, 17]
    geometry =
        point.model == :tfim ?
        "open-chain-n9; sign-symmetric" :
        "kagome-six-triangle-n13; inner-triangles=1,2; edges=[]; sign-symmetric"
    rows = [
        "legacy-basis-metadata-v1",
        "model=$(point.model)",
        "N=$(point.N)",
        "d=$(point.d)",
        "lso=$(point.lso)",
        "geometry=$geometry",
        "positive_blocks=$(join(positive_blocks, ','))",
        "gap_blocks=$(join(gap_blocks, ','))",
        "spectralgap_source_commit=$EXPECTED_SPECTRALGAP_COMMIT",
        "spectralgap_patch_sha256=$EXPECTED_PATCH_SHA256",
    ]
    return (
        fingerprint=sha256_hex(join(rows, "\n") * "\n"),
        fingerprint_scheme="legacy-basis-metadata-v1",
        fingerprint_scope="declared basis inputs and source-gated expected block sizes",
        positive_block_sizes=positive_blocks,
        gap_block_sizes=gap_blocks,
        runtime_basis_export="unavailable",
    )
end

function point_parameters(point::RunPoint)
    if point.model == :tfim
        return (
            N=point.N,
            boundary="open",
            g=canonical_decimal(something(point.g)),
            d=point.d,
            lso=point.lso,
            gamma=canonical_decimal(point.gamma),
            symmetry="sign-symmetric",
        )
    end
    return (
        N=point.N,
        triangles=[collect(triple) for triple in KAGOME_TRIPLES],
        edges=Any[],
        inner_triangles=[collect(triple) for triple in KAGOME_INNER_TRIPLES],
        inner_edges=Any[],
        spin_normalization="S=sigma/2",
        d=point.d,
        lso=point.lso,
        gamma=canonical_decimal(point.gamma),
        symmetry="sign-symmetric",
    )
end

function unavailable(reason::AbstractString)
    return (availability="unavailable", reason=String(reason))
end

function value_or_nothing(value)
    if value isa AbstractFloat
        return isfinite(value) ? value : nothing
    end
    return value
end

function adapt_certify_result(raw)
    if raw isa NamedTuple
        required = (:flag, :termination, :primal, :dual, :objective)
        missing = [name for name in required if !hasproperty(raw, name)]
        isempty(missing) ||
            throw(
                ArgumentError(
                    "NamedTuple certify result is missing fields: " *
                    join(string.(missing), ", "),
                ),
            )
        return (
            adapter="namedtuple-v1",
            flag=getproperty(raw, :flag),
            termination=string(getproperty(raw, :termination)),
            primal=string(getproperty(raw, :primal)),
            dual=string(getproperty(raw, :dual)),
            objective=value_or_nothing(getproperty(raw, :objective)),
            objective_availability=(
                hasproperty(raw, :objective_available) ?
                (
                    getproperty(raw, :objective_available) ?
                    "available" : "unavailable"
                ) :
                (
                    isnothing(value_or_nothing(getproperty(raw, :objective))) ?
                    "unavailable" : "available"
                )
            ),
        )
    end
    if raw isa Integer
        return (
            adapter="legacy-int-explicit-unavailable-status",
            flag=Int(raw),
            termination="UNAVAILABLE_LEGACY_INT",
            primal="UNAVAILABLE_LEGACY_INT",
            dual="UNAVAILABLE_LEGACY_INT",
            objective=nothing,
            objective_availability="unavailable-legacy-api",
        )
    end
    throw(
        ArgumentError(
            "unsupported certify return type $(typeof(raw)); expected NamedTuple or legacy Integer",
        ),
    )
end

function exception_record(err, bt)
    rendered = sprint(showerror, err, bt)
    return (
        type=string(typeof(err)),
        message=sprint(showerror, err),
        backtrace_sha256=sha256_hex(rendered),
        stacktrace=rendered,
    )
end

function runtime_metadata()
    run_spec_path = get(ENV, "HARNESS_RUN_SPEC", nothing)
    absolute_run_spec =
        isnothing(run_spec_path) ? nothing : abspath(run_spec_path)
    julia_executable =
        realpath(joinpath(Sys.BINDIR, Base.julia_exename()))
    return (
        host=gethostname(),
        julia_version=string(VERSION),
        julia_executable=julia_executable,
        julia_executable_sha256=file_sha256(julia_executable),
        startup_file_disabled=Base.JLOptions().startupfile == 2,
        history_file_disabled=Base.JLOptions().historyfile == 0,
        slurm_job_id=get(ENV, "SLURM_JOB_ID", nothing),
        slurm_array_job_id=get(ENV, "SLURM_ARRAY_JOB_ID", nothing),
        slurm_array_task_id=get(ENV, "SLURM_ARRAY_TASK_ID", nothing),
        cpus_per_task=get(ENV, "SLURM_CPUS_PER_TASK", nothing),
        julia_num_threads=string(Threads.nthreads()),
        run_spec_path=absolute_run_spec,
        run_spec_sha256=(
            !isnothing(absolute_run_spec) && isfile(absolute_run_spec) ?
            file_sha256(absolute_run_spec) : nothing
        ),
    )
end

function result_record(
    point::RunPoint,
    adapted,
    elapsed::Real,
    source;
    exception=nothing,
)
    return (
        schema_version=SCHEMA_VERSION,
        record_type="solver_result",
        point_id=point_id(point),
        model=string(point.model),
        parameters=point_parameters(point),
        flag=adapted.flag,
        termination=adapted.termination,
        primal=adapted.primal,
        dual=adapted.dual,
        objective=adapted.objective,
        objective_availability=adapted.objective_availability,
        walltime=Float64(elapsed),
        walltime_unit="seconds",
        exception=exception,
        result_adapter=adapted.adapter,
        residual=unavailable(
            "the current SpectralGap certify API does not export independently recomputable residuals",
        ),
        witness=unavailable(
            "the current SpectralGap certify API does not export a Farkas/infeasibility witness",
        ),
        source=source,
        basis=basis_metadata(point),
        hamiltonian=hamiltonian_metadata(point),
        runtime=runtime_metadata(),
        recorded_at=Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"),
    )
end

function execute_point(
    point::RunPoint,
    certifier::Function;
    source=(mode="fixture",),
)
    start_ns = time_ns()
    try
        raw = certifier(point)
        adapted = adapt_certify_result(raw)
        elapsed = (time_ns() - start_ns) / 1.0e9
        return result_record(point, adapted, elapsed, source)
    catch err
        bt = catch_backtrace()
        elapsed = (time_ns() - start_ns) / 1.0e9
        adapted = (
            adapter="exception",
            flag=nothing,
            termination="UNAVAILABLE_EXCEPTION",
            primal="UNAVAILABLE_EXCEPTION",
            dual="UNAVAILABLE_EXCEPTION",
            objective=nothing,
            objective_availability="unavailable-exception",
        )
        return result_record(
            point,
            adapted,
            elapsed,
            source;
            exception=exception_record(err, bt),
        )
    end
end

function git_read(repo::AbstractString, args::Vector{String})
    try
        return readchomp(Cmd(["git", "-C", repo, args...]))
    catch
        return nothing
    end
end

function environment_metadata(environment_dir::AbstractString)
    project = joinpath(environment_dir, "Project.toml")
    manifest = joinpath(environment_dir, "Manifest.toml")
    return (
        lock_state=isfile(manifest) ? "locked-by-manifest" : "unlocked-manifest-missing",
        project_path=abspath(project),
        project_sha256=isfile(project) ? file_sha256(project) : nothing,
        manifest_path=abspath(manifest),
        manifest_sha256=isfile(manifest) ? file_sha256(manifest) : nothing,
    )
end

function source_metadata(
    repo_root::AbstractString,
    environment_dir::AbstractString;
    verify::Bool,
)
    @eval import SpectralGap
    module_file = pathof(SpectralGap)
    source_dir = dirname(module_file)
    package_root = dirname(source_dir)
    files = Dict{String, Union{Nothing, String}}()
    for filename in sort(collect(keys(EXPECTED_SPECTRALGAP_FILES)))
        path = joinpath(source_dir, filename)
        files[filename] = isfile(path) ? file_sha256(path) : nothing
    end
    patch_path = joinpath(
        repo_root,
        "tracks",
        "polyopt",
        "solutions",
        "sdp-gap-seekers",
        "spectralgap_a1171c9.patch",
    )
    patch_sha = isfile(patch_path) ? file_sha256(patch_path) : nothing
    package_commit = git_read(package_root, ["rev-parse", "HEAD"])
    package_tree = git_read(package_root, ["rev-parse", "HEAD^{tree}"])
    harness_commit = git_read(repo_root, ["rev-parse", "HEAD"])
    harness_tree = git_read(repo_root, ["rev-parse", "HEAD^{tree}"])
    locked_harness_commit =
        get(ENV, "GAP_STATUS_EXPECTED_HARNESS_COMMIT", nothing)
    locked_harness_tree =
        get(ENV, "GAP_STATUS_EXPECTED_HARNESS_TREE", nothing)
    exact_locked_harness =
        !isnothing(locked_harness_commit) &&
        !isnothing(locked_harness_tree) &&
        harness_commit == locked_harness_commit &&
        harness_tree == locked_harness_tree
    julia_executable =
        realpath(joinpath(Sys.BINDIR, Base.julia_exename()))
    julia_version = "julia version $(VERSION)"
    julia_sha256 = file_sha256(julia_executable)
    locked_julia_executable =
        get(ENV, "GAP_STATUS_EXPECTED_JULIA_BIN", nothing)
    locked_julia_version =
        get(ENV, "GAP_STATUS_EXPECTED_JULIA_VERSION", nothing)
    locked_julia_sha256 =
        get(ENV, "GAP_STATUS_EXPECTED_JULIA_SHA256", nothing)
    exact_locked_julia =
        julia_executable == locked_julia_executable &&
        julia_version == locked_julia_version &&
        julia_sha256 == locked_julia_sha256 &&
        Base.JLOptions().startupfile == 2 &&
        Base.JLOptions().historyfile == 0
    harness_base_tree =
        git_read(repo_root, ["rev-parse", "$EXPECTED_HARNESS_BASE^{tree}"])
    base_is_ancestor =
        git_read(
            repo_root,
            [
                "merge-base",
                "--is-ancestor",
                EXPECTED_HARNESS_BASE,
                something(harness_commit, "HEAD"),
            ],
        ) !== nothing
    relevant_files = Dict{String, Union{Nothing, String}}()
    relevant_tracked = Dict{String, Bool}()
    for relative_path in RELEVANT_HARNESS_PATHS
        path = joinpath(repo_root, relative_path)
        relevant_files[relative_path] =
            isfile(path) ? file_sha256(path) : nothing
        relevant_tracked[relative_path] =
            git_read(
                repo_root,
                ["ls-files", "--error-unmatch", "--", relative_path],
            ) !== nothing
    end
    relevant_status = something(
        git_read(
            repo_root,
            [
                "status",
                "--porcelain",
                "--untracked-files=all",
                "--",
                RELEVANT_HARNESS_PATHS...,
            ],
        ),
        "git-status-error",
    )
    relevant_files_clean =
        isempty(relevant_status) &&
        all(values(relevant_tracked)) &&
        all(!isnothing, values(relevant_files))
    env = environment_metadata(environment_dir)
    file_match = all(
        get(files, filename, nothing) == expected
        for (filename, expected) in EXPECTED_SPECTRALGAP_FILES
    )
    exact_match =
        patch_sha == EXPECTED_PATCH_SHA256 &&
        package_commit == EXPECTED_SPECTRALGAP_COMMIT &&
        package_tree == EXPECTED_SPECTRALGAP_TREE &&
        file_match &&
        base_is_ancestor &&
        harness_base_tree == EXPECTED_HARNESS_TREE &&
        exact_locked_harness &&
        exact_locked_julia &&
        relevant_files_clean &&
        env.lock_state == "locked-by-manifest"
    source_fingerprint = sha256_hex(
        join(
            [
                "gap-status-source-v1",
                "harness_commit=$(something(harness_commit, "unavailable"))",
                "harness_tree=$(something(harness_tree, "unavailable"))",
                "locked_harness_commit=$(something(locked_harness_commit, "unavailable"))",
                "locked_harness_tree=$(something(locked_harness_tree, "unavailable"))",
                "julia_executable=$julia_executable",
                "julia_version=$julia_version",
                "julia_sha256=$julia_sha256",
                "harness_base_tree=$(something(harness_base_tree, "unavailable"))",
                [
                    "$relative_path=$(something(relevant_files[relative_path], "unavailable"))"
                    for relative_path in sort(collect(keys(relevant_files)))
                ]...,
                "spectralgap_commit=$(something(package_commit, "unavailable"))",
                "spectralgap_tree=$(something(package_tree, "unavailable"))",
                "patch=$(something(patch_sha, "unavailable"))",
                [
                    "$filename=$(something(files[filename], "unavailable"))"
                    for filename in sort(collect(keys(files)))
                ]...,
                "project=$(something(env.project_sha256, "unavailable"))",
                "manifest=$(something(env.manifest_sha256, "unavailable"))",
            ],
            "\n",
        ) * "\n",
    )
    metadata = (
        fingerprint=source_fingerprint,
        fingerprint_scheme="gap-status-source-v1",
        exact_source_match=exact_match,
        harness=(
            commit=harness_commit,
            tree=harness_tree,
            required_base_commit=EXPECTED_HARNESS_BASE,
            required_base_tree=EXPECTED_HARNESS_TREE,
            required_base_is_ancestor=base_is_ancestor,
            actual_required_base_tree=harness_base_tree,
            locked_commit=locked_harness_commit,
            locked_tree=locked_harness_tree,
            exact_locked_harness=exact_locked_harness,
            relevant_files=relevant_files,
            relevant_files_tracked=relevant_tracked,
            relevant_files_clean=relevant_files_clean,
            relevant_status=relevant_status,
        ),
        julia=(
            executable=julia_executable,
            version=julia_version,
            sha256=julia_sha256,
            locked_executable=locked_julia_executable,
            locked_version=locked_julia_version,
            locked_sha256=locked_julia_sha256,
            startup_file_disabled=Base.JLOptions().startupfile == 2,
            history_file_disabled=Base.JLOptions().historyfile == 0,
            exact_locked_julia=exact_locked_julia,
        ),
        spectralgap=(
            package_root=abspath(package_root),
            base_commit=package_commit,
            base_tree=package_tree,
            expected_base_commit=EXPECTED_SPECTRALGAP_COMMIT,
            expected_base_tree=EXPECTED_SPECTRALGAP_TREE,
            files=files,
            expected_files=EXPECTED_SPECTRALGAP_FILES,
        ),
        patch=(
            path=abspath(patch_path),
            sha256=patch_sha,
            expected_sha256=EXPECTED_PATCH_SHA256,
        ),
        environment=env,
    )
    if verify && !exact_match
        throw(
            ArgumentError(
                "source/environment gate failed; inspect preflight JSON (no solver was called)",
            ),
        )
    end
    return metadata
end

function source_metadata_or_error(
    repo_root::AbstractString,
    environment_dir::AbstractString,
)
    try
        return source_metadata(repo_root, environment_dir; verify=false), nothing
    catch err
        bt = catch_backtrace()
        fallback = (
            exact_source_match=false,
            source_load_error=exception_record(err, bt),
            expected=(
                harness_base_commit=EXPECTED_HARNESS_BASE,
                harness_base_tree=EXPECTED_HARNESS_TREE,
                spectralgap_base_commit=EXPECTED_SPECTRALGAP_COMMIT,
                spectralgap_base_tree=EXPECTED_SPECTRALGAP_TREE,
                patch_sha256=EXPECTED_PATCH_SHA256,
                spectralgap_files=EXPECTED_SPECTRALGAP_FILES,
            ),
            environment=environment_metadata(environment_dir),
        )
        return fallback, err
    end
end

function build_dry_run_record(
    point::RunPoint,
    repo_root::AbstractString,
    environment_dir::AbstractString,
)
    patch_path = joinpath(
        repo_root,
        "tracks",
        "polyopt",
        "solutions",
        "sdp-gap-seekers",
        "spectralgap_a1171c9.patch",
    )
    return (
        schema_version=SCHEMA_VERSION,
        record_type="dry_run_plan",
        point_id=point_id(point),
        model=string(point.model),
        parameters=point_parameters(point),
        solver_called=false,
        expected_source=(
            harness_base_commit=EXPECTED_HARNESS_BASE,
            harness_base_tree=EXPECTED_HARNESS_TREE,
            spectralgap_base_commit=EXPECTED_SPECTRALGAP_COMMIT,
            spectralgap_base_tree=EXPECTED_SPECTRALGAP_TREE,
            patch_sha256=EXPECTED_PATCH_SHA256,
            local_patch_sha256=isfile(patch_path) ? file_sha256(patch_path) : nothing,
            spectralgap_files=EXPECTED_SPECTRALGAP_FILES,
        ),
        environment=environment_metadata(environment_dir),
        basis=basis_metadata(point),
        hamiltonian=hamiltonian_metadata(point),
    )
end

function actual_certifier(point::RunPoint)
    supports, coefficients = build_hamiltonian_data(point)
    H = SpectralGap.ncpoly(supports, Float64.(coefficients))
    if point.model == :tfim
        return SpectralGap.certify_Ising_gap(
            point.N,
            H,
            Float64(point.gamma),
            point.d;
            lso=point.lso,
            QUIET=true,
        )
    end
    triples = [collect(triple) for triple in KAGOME_TRIPLES]
    inner_triples = [collect(triple) for triple in KAGOME_INNER_TRIPLES]
    return SpectralGap.certify_Heisenberg_kagome_gap(
        point.N,
        H,
        triples,
        Any[],
        inner_triples,
        Any[],
        Float64(point.gamma),
        point.d;
        lso=point.lso,
        QUIET=true,
    )
end

function json_escape(text::AbstractString)
    io = IOBuffer()
    for character in text
        if character == '"'
            print(io, "\\\"")
        elseif character == '\\'
            print(io, "\\\\")
        elseif character == '\b'
            print(io, "\\b")
        elseif character == '\f'
            print(io, "\\f")
        elseif character == '\n'
            print(io, "\\n")
        elseif character == '\r'
            print(io, "\\r")
        elseif character == '\t'
            print(io, "\\t")
        elseif Int(character) < 0x20
            print(io, "\\u", lpad(string(Int(character); base=16), 4, '0'))
        else
            print(io, character)
        end
    end
    return String(take!(io))
end

function json_encode(value)
    if value === nothing || value === missing
        return "null"
    elseif value isa Bool
        return value ? "true" : "false"
    elseif value isa Integer
        return string(value)
    elseif value isa AbstractFloat
        return isfinite(value) ? repr(Float64(value)) : "null"
    elseif value isa AbstractString || value isa Symbol
        return "\"$(json_escape(string(value)))\""
    elseif value isa NamedTuple
        entries = [
            string(json_encode(string(name)), ":", json_encode(getproperty(value, name)))
            for name in propertynames(value)
        ]
        return "{" * join(entries, ",") * "}"
    elseif value isa AbstractDict
        entries = [
            string(json_encode(string(key)), ":", json_encode(value[key]))
            for key in sort(collect(keys(value)); by=string)
        ]
        return "{" * join(entries, ",") * "}"
    elseif value isa Tuple || value isa AbstractVector
        return "[" * join(json_encode.(collect(value)), ",") * "]"
    end
    return json_encode(string(value))
end

function write_record(record, output::Union{Nothing, AbstractString})
    line = json_encode(record) * "\n"
    if isnothing(output)
        print(stdout, line)
        flush(stdout)
        return
    end
    directory = dirname(abspath(output))
    mkpath(directory)
    temporary, io = mktemp(directory)
    try
        write(io, line)
        flush(io)
        close(io)
        mv(temporary, output; force=true)
    catch
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
        rethrow()
    end
end

function usage(io::IO=stdout)
    println(
        io,
        """
Usage:
  julia --project=ENV gap_status_runner.jl --run-spec PATH --cell-index 1 [--dry-run|--preflight] [--output PATH]
  julia --project=julia-env gap_status_runner.jl --model tfim --gamma 0.25 [--dry-run|--preflight] [--output PATH]
  julia --project=julia-env gap_status_runner.jl --model kagome --gamma 1.28 [--dry-run|--preflight] [--output PATH]
  julia gap_status_runner.jl --list-points

Locked cells:
  TFIM:   N=9, g=0.5, open, d=2, lso=6, gamma in {0.25,0.26}
  Kagome: N=13, six triangles, inner first two, d=3, lso=5,
          gamma in {1,1.2,1.26,1.28}

--dry-run validates and fingerprints inputs without importing SpectralGap or
calling Mosek. --preflight imports SpectralGap and verifies source/patch/
Manifest fingerprints without constructing or optimizing a model.
""",
    )
end

function parse_cli(args)
    values = Dict{String, String}()
    switches = Set{String}()
    value_options = Set([
        "--model",
        "--gamma",
        "--N",
        "--g",
        "--d",
        "--lso",
        "--output",
        "--repo-root",
        "--environment-dir",
        "--run-spec",
        "--cell-index",
    ])
    switch_options = Set(["--dry-run", "--preflight", "--list-points", "--help"])
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument in switch_options
            push!(switches, argument)
            index += 1
        elseif argument in value_options
            index < length(args) ||
                throw(ArgumentError("$argument requires a value"))
            haskey(values, argument) &&
                throw(ArgumentError("$argument may only be specified once"))
            values[argument] = args[index + 1]
            index += 2
        else
            throw(ArgumentError("unknown argument '$argument'"))
        end
    end
    return values, switches
end

function default_repo_root()
    return normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
end

function main(args=ARGS)
    values, switches = parse_cli(args)
    if "--help" in switches
        usage()
        return 0
    end
    modes = count(
        mode -> mode in switches,
        ("--dry-run", "--preflight", "--list-points"),
    )
    modes <= 1 ||
        throw(ArgumentError("choose at most one of --dry-run, --preflight, --list-points"))

    repo_root = abspath(get(values, "--repo-root", default_repo_root()))
    environment_dir =
        abspath(get(values, "--environment-dir", joinpath(repo_root, "julia-env")))
    output = get(values, "--output", nothing)

    if "--list-points" in switches
        point_options = intersect(
            Set(keys(values)),
            Set([
                "--model",
                "--gamma",
                "--N",
                "--g",
                "--d",
                "--lso",
                "--run-spec",
                "--cell-index",
            ]),
        )
        isempty(point_options) ||
            throw(ArgumentError("--list-points does not accept point parameters"))
        records = [
            (
                schema_version=SCHEMA_VERSION,
                record_type="allowed_point",
                point_id=point_id(point),
                model=string(point.model),
                parameters=point_parameters(point),
            ) for point in allowed_points()
        ]
        if isnothing(output)
            for record in records
                write_record(record, nothing)
            end
        else
            mkpath(dirname(abspath(output)))
            open(output, "w") do io
                for record in records
                    println(io, json_encode(record))
                end
            end
        end
        return 0
    end

    haskey(values, "--model") ||
        haskey(values, "--run-spec") ||
        throw(ArgumentError("provide --run-spec/--cell-index or --model/--gamma"))
    run_spec_mode = haskey(values, "--run-spec") || haskey(values, "--cell-index")
    if run_spec_mode
        haskey(values, "--run-spec") && haskey(values, "--cell-index") ||
            throw(ArgumentError("--run-spec and --cell-index must be used together"))
        point_args = intersect(
            Set(keys(values)),
            Set(["--model", "--gamma", "--N", "--g", "--d", "--lso"]),
        )
        isempty(point_args) ||
            throw(ArgumentError("run-spec mode does not accept point parameters"))
        cell_index = parse(Int, values["--cell-index"])
        ENV["HARNESS_RUN_SPEC"] = abspath(values["--run-spec"])
        ENV["GAP_STATUS_CELL_INDEX"] = string(cell_index)
        point, _ = parse_run_spec_point(values["--run-spec"], cell_index)
    else
        haskey(values, "--model") ||
            throw(ArgumentError("--model is required"))
        haskey(values, "--gamma") ||
            throw(ArgumentError("--gamma is required"))
        point = parse_point(
            values["--model"],
            values["--gamma"];
            N_text=get(values, "--N", nothing),
            g_text=get(values, "--g", nothing),
            d_text=get(values, "--d", nothing),
            lso_text=get(values, "--lso", nothing),
        )
    end

    if "--dry-run" in switches
        write_record(build_dry_run_record(point, repo_root, environment_dir), output)
        return 0
    end

    if "--preflight" in switches
        source, _ = source_metadata_or_error(repo_root, environment_dir)
        record = (
            schema_version=SCHEMA_VERSION,
            record_type="source_preflight",
            point_id=point_id(point),
            solver_called=false,
            source=source,
        )
        write_record(record, output)
        return source.exact_source_match ? 0 : 2
    end

    source, source_error =
        source_metadata_or_error(repo_root, environment_dir)
    if !source.exact_source_match
        reason =
            isnothing(source_error) ?
            ArgumentError(
                "source/environment gate failed; no optimizer was constructed",
            ) : source_error
        record = execute_point(
            point,
            _ -> throw(reason);
            source=source,
        )
        write_record(record, output)
        return 1
    end
    println(
        "START point=$(point_id(point)) at ",
        Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
    )
    flush(stdout)
    record = execute_point(point, actual_certifier; source=source)
    write_record(record, output)
    println(
        "DONE point=$(point_id(point)) termination=$(record.termination) result=$(isnothing(output) ? "stdout" : output)",
    )
    flush(stdout)
    return isnothing(record.exception) ? 0 : 1
end

end # module
