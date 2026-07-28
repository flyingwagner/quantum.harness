module SquareStatusEnvelope

using SHA
using ..SharedCoreWire:
    CanonicalObject,
    canonical_framed_bytes

export build_square_status_envelope,
    write_square_status_envelope

object(fields::Pair...) = CanonicalObject(fields...)

function require_sha256(value::AbstractString, label::AbstractString)
    occursin(r"^[0-9a-f]{64}$", value) ||
        error("$label is not a lowercase SHA-256 digest")
    return String(value)
end

function require_git_id(value::AbstractString, label::AbstractString)
    occursin(r"^[0-9a-f]{40}$", value) ||
        error("$label is not a full lowercase Git object ID")
    return String(value)
end

"""
Construct a canonical envelope for a solver-free Square conic status gate.

The caller supplies independently validated core and MOF audit records. This
function fails closed unless both records explicitly report complete
validation and no optimizer invocation. The envelope deliberately records
`status=unsolved`; it cannot itself support a gap bound.
"""
function build_square_status_envelope(
    identity;
    source_commit::AbstractString,
    source_tree::AbstractString,
    project_sha256::AbstractString,
    manifest_sha256::AbstractString,
    core_validation,
    mof_sha256::AbstractString,
    mof_audit,
)
    core_validation.canonical_roundtrip === true ||
        error("core canonical roundtrip did not pass")
    core_validation.envelope_valid === true ||
        error("core envelope validation did not pass")
    core_validation.derived_ids_valid === true ||
        error("core derived-ID validation did not pass")
    core_validation.coverage_valid === true ||
        error("core coverage validation did not pass")
    core_validation.source_rebuild_match === true ||
        error("core artifact does not match the current source rebuild")
    core_validation.optimizer_invoked === false ||
        error("core validation unexpectedly invoked an optimizer")
    mof_audit.exact_coefficient_match === true ||
        error("MOF exact coefficient audit did not pass")
    mof_audit.optimizer_invoked === false ||
        error("MOF audit unexpectedly invoked an optimizer")
    mof_audit.objective_sense == :feasibility ||
        error("Square MOF is not a feasibility problem")
    core_validation.scalar_rows == identity.scalar_rows ||
        error("core scalar-row count differs from source identity")
    mof_audit.variables == identity.scalar_rows ||
        error("MOF variable count differs from source identity")
    mof_audit.affine_equalities == identity.affine_equalities ||
        error("MOF equality count differs from source identity")
    mof_audit.psd_dimensions == identity.psd_real_dimensions ||
        error("MOF PSD dimensions differ from source identity")

    envelope = object(
        "schema_id" => "square-gap-status-envelope",
        "schema_version" => BigInt(1),
        "status" => "unsolved",
        "claim_boundary" =>
            "source-gated conic model; no Square bulk-gap bound",
        "source" => object(
            "commit" => require_git_id(source_commit, "source commit"),
            "tree" => require_git_id(source_tree, "source tree"),
        ),
        "environment" => object(
            "project_sha256" =>
                require_sha256(project_sha256, "Project.toml"),
            "manifest_sha256" =>
                require_sha256(manifest_sha256, "Manifest.toml"),
        ),
        "physical_setup" => object(
            "model" => "square-j1-j2",
            "hamiltonian" =>
                "H=1/4 sum_J1(XX+YY+ZZ)+(g/4) sum_J2(XX+YY+ZZ)",
            "l" => BigInt(identity.L),
            "d" => BigInt(identity.d),
            "g" => identity.g,
            "gamma" => identity.gamma,
            "state_class" => identity.state_class,
            "symmetry" => "unrestricted; no quotient",
            "volume_semantics" =>
                "local consistency window for infinite-volume KMS ground states",
        ),
        "source_identity" => object(
            "problem_sha256" =>
                require_sha256(identity.problem_sha256, "problem"),
            "positive_basis_sha256" => require_sha256(
                identity.positive_basis_sha256,
                "positive basis",
            ),
            "gap_basis_sha256" =>
                require_sha256(identity.gap_basis_sha256, "gap basis"),
        ),
        "exact_core" => object(
            "math_sha256" => require_sha256(
                core_validation.math_sha256,
                "core math",
            ),
            "envelope_sha256" => require_sha256(
                core_validation.envelope_sha256,
                "core envelope",
            ),
            "math_byte_count" => BigInt(core_validation.math_byte_count),
            "blocks" => BigInt(core_validation.blocks),
            "scalar_rows" => BigInt(core_validation.scalar_rows),
            "wiring_pairs" => BigInt(core_validation.wiring_pairs),
            "component_records" =>
                BigInt(core_validation.component_records),
            "nonzero_coefficients" =>
                BigInt(core_validation.nonzero_coefficients),
        ),
        "rendered_conic" => object(
            "mof_sha256" => require_sha256(mof_sha256, "MOF"),
            "variables" => BigInt(mof_audit.variables),
            "affine_equalities" => BigInt(mof_audit.affine_equalities),
            "affine_coefficients" =>
                BigInt(mof_audit.affine_coefficients),
            "psd_real_dimensions" =>
                BigInt.(mof_audit.psd_dimensions),
            "psd_coefficients" => BigInt(mof_audit.psd_coefficients),
            "objective_sense" => "feasibility",
            "exact_coefficient_match" => true,
        ),
        "audit" => object(
            "core_validated" => true,
            "core_source_rebuilt" => true,
            "mof_exactly_replayed" => true,
            "optimizer_invoked" => false,
        ),
    )
    bytes = canonical_framed_bytes("AISQSTATUS1", envelope)
    return (
        envelope=envelope,
        bytes=bytes,
        sha256=bytes2hex(sha256(bytes)),
    )
end

function write_square_status_envelope(path::AbstractString, result)
    ispath(path) && error("refusing to overwrite status envelope: $path")
    open(path, "x") do io
        write(io, result.bytes)
    end
    return path
end

end
