module SquareCoreInventory

using SHA
using ..SquareJ1J2Prototype: PauliWord
using ..GenericGapModel:
    BasisManifest,
    GapProblem,
    StateMonomial,
    validate_basis_manifest
import ..CoreMGK
using ..CoreMGK:
    CoreMGKPlan,
    GaussianRational,
    ScalarMoment,
    core_mgk_pair,
    core_mgk_plan
import ..SharedCoreWire
using ..SharedCoreWire:
    CanonicalObject,
    canonical_framed_bytes,
    canonical_scalar_row,
    canonical_state_monomial,
    canonical_support_word,
    decode_framed_bytes,
    domain_sha256,
    encode_value,
    entry_id,
    row_id,
    term_id

const BigRational = Rational{BigInt}

export build_square_core_inventory,
    validate_square_core_inventory,
    write_square_core_inventory

object(fields::Pair...) = CanonicalObject(fields...)
exact_rational(value::Integer) = BigInt(value) // BigInt(1)
exact_rational(value::Rational) =
    BigInt(numerator(value)) // BigInt(denominator(value))
exact_rational(value::AbstractFloat) =
    throw(ArgumentError("shared core inventory requires exact rational inputs"))

gaussian(value::Complex) = object(
    "real" => exact_rational(real(value)),
    "imag" => exact_rational(imag(value)),
)
gaussian(value::Real) = gaussian(complex(value, zero(value)))

function byte_vector_less(left, right)
    for index in 1:min(length(left), length(right))
        left[index] == right[index] || return left[index] < right[index]
    end
    return length(left) < length(right)
end

encoded_less(left, right) =
    byte_vector_less(encode_value(left), encode_value(right))

function canonical_hamiltonian(source::CoreMGKPlan)
    accumulated = Dict{PauliWord,GaussianRational}()
    for term in source.hamiltonian_terms
        coefficient = GaussianRational(
            exact_rational(real(term.coefficient)),
            exact_rational(imag(term.coefficient)),
        )
        value = get(accumulated, term.word, zero(GaussianRational)) + coefficient
        iszero(value) ? delete!(accumulated, term.word) : (accumulated[term.word] = value)
    end
    any(!iszero(imag(value)) for value in values(accumulated)) &&
        error("Square Hamiltonian is not exactly Hermitian")
    records = [
        object(
            "term_id" => term_id(word),
            "support_word" => canonical_support_word(word),
            "coefficient" => gaussian(coefficient),
        )
        for (word, coefficient) in accumulated
    ]
    sort!(records; lt=(left, right) -> encoded_less(
        only(filter(pair -> first(pair) == "support_word", left.fields)).second,
        only(filter(pair -> first(pair) == "support_word", right.fields)).second,
    ))
    return records
end

function selector_spec(manifest::BasisManifest)
    selector_parameters = object(
        "finite_formal_basis_complete" => manifest.is_complete,
        "max_degree" => BigInt(manifest.max_degree),
        "role" => string(manifest.role),
        "site_ids" => BigInt.(manifest.site_ids),
    )
    selection_rule = object(
        "rule_id" => "structured-one-symbol-lift-v1",
        "rule_version" => BigInt(1),
        "parameters" => object(
            "family" => string(manifest.family),
            "family_version" => BigInt(manifest.family_version),
            "max_degree" => BigInt(manifest.max_degree),
            "role" => string(manifest.role),
            "site_ids" => BigInt.(manifest.site_ids),
        ),
    )
    preimage = object(
        "basis_selector_id" => "structured-one-symbol-lift",
        "basis_selector_version" => BigInt(1),
        "selector_parameters" => selector_parameters,
        "selection_rule" => selection_rule,
    )
    digest = domain_sha256("AISEL1", preimage)
    return object(
        preimage.fields...,
        "basis_selector_sha256" => digest,
    )
end

function field(object_value::CanonicalObject, key::AbstractString)
    matches = [last(pair) for pair in object_value.fields if first(pair) == key]
    return only(matches)
end

function without_field(object_value::CanonicalObject, key::AbstractString)
    return CanonicalObject(Pair{String,Any}[
        pair for pair in object_value.fields if first(pair) != key
    ])
end

function unrestricted_action(manifests::Vector{BasisManifest})
    ids = unique(vcat(
        ([entry_id(entry) for entry in manifest.entries] for manifest in manifests)...,
    ))
    sort!(ids)
    word_rule_entries = Any[
        object(
            "input_entry_id" => id,
            "output_entry_id" => id,
            "phase" => gaussian(one(BigRational)),
            "forced_zero" => false,
            "zero_reason" => nothing,
        )
        for id in ids
    ]
    preimage = object(
        "action_id" => "unrestricted-identity",
        "action_version" => BigInt(1),
        "mode" => "unrestricted",
        "generators" => Any[],
        "word_rule" => object(
            "rule_schema_id" => "explicit-word-rule-v1",
            "rule_schema_version" => BigInt(1),
            "input_domain_entry_ids" => ids,
            "entries" => word_rule_entries,
        ),
        "physical_group_claim" => nothing,
        "proof_status" => "not_claimed",
    )
    digest = domain_sha256("AIACT1", preimage)
    action = object(preimage.fields..., "action_sha256" => digest)
    reference = object(
        "action_id" => field(action, "action_id"),
        "action_version" => field(action, "action_version"),
        "action_sha256" => digest,
    )
    return action, reference
end

function action_set_sha256(actions)
    refs = [
        object(
            "action_id" => field(action, "action_id"),
            "action_version" => field(action, "action_version"),
            "action_sha256" => field(action, "action_sha256"),
        )
        for action in actions
    ]
    sort!(refs; lt=encoded_less)
    return domain_sha256("AIACTSET1", refs)
end

function block_selector_spec(
    positive_ids::Vector{String},
    gap_ids::Vector{String},
    action_reference::CanonicalObject,
)
    rules = Any[
        object(
            "role" => "positive",
            "block_label" => "all",
            "instance_key" => "main",
            "state_action" => action_reference,
            "selector_kind" => "explicit_ordered_entries",
            "selector_payload" => positive_ids,
        ),
        object(
            "role" => "gap",
            "block_label" => "all",
            "instance_key" => "main",
            "state_action" => action_reference,
            "selector_kind" => "explicit_ordered_entries",
            "selector_payload" => gap_ids,
        ),
    ]
    preimage = object(
        "block_selector_id" => "square-flat-positive-gap",
        "block_selector_version" => BigInt(1),
        "block_selector_parameters" => object(
            "decomposition_mode" => "flat",
        ),
        "block_selection_rule" => object(
            "rule_id" => "explicit-flat-positive-gap-v1",
            "rule_version" => BigInt(1),
            "parameters" => object(),
        ),
        "labelled_block_rules" => rules,
    )
    return object(
        preimage.fields...,
        "block_selector_sha256" => domain_sha256("AIBLOCKSEL1", preimage),
    )
end

coefficient_convention() = object(
    "mgk_definition_id" => "symmetrized-k-covariance-g-v1",
    "matrix_orientation_id" => "row-j-column-k-v1",
    "inner_product_id" => "trace-aq-v1",
    "packing_id" => "plus-imag-upper-triangle-v1",
    "gamma_domain_id" => "exact-nonnegative-rational-v1",
)

normalization_spec() = object(
    "normalization_id" => "state-functional-unit-v1",
    "normalization_version" => BigInt(1),
    "parameters" => object("identity_value" => one(BigRational)),
)

function exact_model_config(problem::GapProblem)
    tagged = Dict{Symbol,Set{BigRational}}()
    for template in problem.model.templates
        coefficient = 4 * exact_rational(template.coefficient)
        push!(get!(tagged, template.tag, Set{BigRational}()), coefficient)
    end
    j1 = only(get(tagged, :J1, Set{BigRational}()))
    j2 = only(get(tagged, :J2, Set{BigRational}()))
    j1 == one(BigRational) ||
        throw(ArgumentError("Square inventory requires the J1=1 convention"))
    return object(
        "config_schema_id" => "square-j1-j2-local-window-v1",
        "config_schema_version" => BigInt(1),
        "parameters" => object(
            "j1" => j1,
            "j2" => j2,
            "j2_over_j1" => j2 / j1,
            "outer_patch_kind" => problem.patch.name,
            "outer_patch_level" => BigInt(problem.patch.level),
            "inner_site_ids" => BigInt.(sort(problem.patch.inner_ids)),
            "interaction_range_linf" =>
                BigInt(problem.model.interaction_range_linf),
            "spin_operator_scale" => BigInt(1) // BigInt(2),
        ),
    )
end

function model_site_table(problem::GapProblem)
    return Any[
        object(
            "site_namespace" => "model",
            "site_id" => BigInt(site_id),
            "site_label" => "($(site.x),$(site.y))",
            "coordinate_frame_id" => "square-lattice-xy",
            "coordinates" => BigRational[
                exact_rational(site.x),
                exact_rational(site.y),
            ],
            "metadata" => object(),
        )
        for (site_id, site) in enumerate(problem.patch.sites)
    ]
end

function relaxation_id(
    exact_config,
    normalization,
    sites,
    hamiltonian_terms,
    problem::GapProblem,
    positive_selector,
    gap_selector,
    action_set_digest,
    block_selector,
    convention,
)
    return domain_sha256(
        "AIREL1",
        object(
            "schema_id" => "shared-core-mgk-inventory",
            "schema_version" => BigInt(1),
            "math_scope" => "core_mgk",
            "model_id" => "square-j1-j2",
            "model_version" => BigInt(1),
            "exact_model_config" => exact_config,
            "normalization" => normalization,
            "model_site_table" => sites,
            "exact_hamiltonian_terms" => hamiltonian_terms,
            "degree" => BigInt(problem.d),
            "positive_basis_selector_sha256" =>
                field(positive_selector, "basis_selector_sha256"),
            "gap_basis_selector_sha256" =>
                field(gap_selector, "basis_selector_sha256"),
            "action_set_sha256" => action_set_digest,
            "block_selector_sha256" =>
                field(block_selector, "block_selector_sha256"),
            "coefficient_convention" => convention,
        ),
    )
end

function basis_order_sha256(block_id::String, entry_ids::Vector{String})
    return domain_sha256(
        "AIBO1",
        object(
            "block_id" => block_id,
            "ordered_entry_ids" => entry_ids,
        ),
    )
end

function block_record(
    relaxation::String,
    role::Symbol,
    manifest::BasisManifest,
    selector::CanonicalObject,
    block_selector::CanonicalObject,
    action_reference::CanonicalObject,
    convention::CanonicalObject,
)
    role_text = string(role)
    block_id = "blk:" * domain_sha256(
        "AIBLK1",
        object(
            "relaxation_id" => relaxation,
            "role" => role_text,
            "block_label" => "all",
            "instance_key" => "main",
            "basis_selector_id" => field(selector, "basis_selector_id"),
            "basis_selector_version" => field(
                selector,
                "basis_selector_version",
            ),
        ),
    )
    entry_ids = entry_id.(manifest.entries)
    entries = Any[]
    for (index, entry) in enumerate(manifest.entries)
        payload = canonical_state_monomial(entry)
        push!(entries, object(
            "block_id" => block_id,
            "entry_id" => entry_ids[index],
            "matrix_index" => BigInt(index),
            "state_symbols" => field(payload, "state_symbols"),
            "operator_word" => field(payload, "operator_word"),
            "source_ordinal" => nothing,
        ))
    end
    spec = object(
        "block_id" => block_id,
        "role" => role_text,
        "block_label" => "all",
        "instance_key" => "main",
        "decomposition_mode" => "flat",
        "source_label" => nothing,
        "source_block_ordinal" => nothing,
        "state_action" => action_reference,
        "symmetry_applied" => false,
        "irrep_label" => nothing,
        "irrep_dimension" => nothing,
        "irrep_proof_status" => "not_claimed",
        "field" => "complex_hermitian",
        "dimension" => BigInt(length(entries)),
        "render_policy" => isempty(entries) ? "skip_empty" : "materialize",
        "coefficient_convention" => convention,
        "block_selector_id" => field(block_selector, "block_selector_id"),
        "block_selector_version" =>
            field(block_selector, "block_selector_version"),
        "block_selector_sha256" =>
            field(block_selector, "block_selector_sha256"),
        "basis_selector_id" => field(selector, "basis_selector_id"),
        "basis_selector_version" =>
            field(selector, "basis_selector_version"),
        "selector_parameters" => field(selector, "selector_parameters"),
        "selection_rule" => field(selector, "selection_rule"),
        "basis_selector_sha256" =>
            field(selector, "basis_selector_sha256"),
        "basis_manifest_schema" => "structured-basis-manifest-v1",
        "basis_manifest_sha256" => manifest.sha256,
        "ordered_entry_ids" => entry_ids,
        "basis_order_sha256" => basis_order_sha256(block_id, entry_ids),
    )
    return (
        id=block_id,
        role=role,
        manifest=manifest,
        entry_ids=entry_ids,
        record=object("spec" => spec, "basis_entries" => entries),
    )
end

const COMPONENT_ORDER = Dict(
    :M => 0,
    :K => 1,
    :G_moment => 2,
    :G_product => 3,
)

function report_pair_progress(
    callback,
    progress_every::Int,
    phase::Symbol,
    completed::Int,
    total::Int,
)
    if progress_every > 0 &&
       (completed % progress_every == 0 || completed == total)
        callback((
            phase=phase,
            completed=completed,
            total=total,
        ))
    end
    return nothing
end

function collect_rows(
    blocks,
    source::CoreMGKPlan,
    progress_every::Int,
    progress_callback,
    total_pairs::Int,
)
    origins = Dict{ScalarMoment,Set{Symbol}}()
    component_records = 0
    completed_pairs = 0
    for block in blocks
        dimension = length(block.manifest.entries)
        for j in 1:dimension, k in j:dimension
            wiring = core_mgk_pair(source, block.role, j, k)
            component_records += length(wiring.component_records)
            for component in wiring.component_records
                for coefficient in component.coefficients
                    push!(
                        get!(origins, coefficient.row, Set{Symbol}()),
                        component.component,
                    )
                end
            end
            completed_pairs += 1
            report_pair_progress(
                progress_callback,
                progress_every,
                :collect_rows,
                completed_pairs,
                total_pairs,
            )
        end
    end
    rows = collect(keys(origins))
    sort!(rows; lt=(left, right) -> begin
        left_bytes = encode_value(canonical_scalar_row(left))
        right_bytes = encode_value(canonical_scalar_row(right))
        left_bytes == right_bytes ?
            row_id(left) < row_id(right) :
            byte_vector_less(left_bytes, right_bytes)
    end)
    return rows, origins, component_records
end

function scalar_row_records(rows, origins, row_ids)
    records = Any[]
    for row in rows
        payload = canonical_scalar_row(row)
        push!(records, object(
            "row_id" => row_ids[row],
            "state_symbol_multiset" =>
                field(payload, "state_symbol_multiset"),
            "scope" => "core_mgk",
            "origins" => [
                string(origin)
                for origin in sort!(
                    collect(origins[row]);
                    by=origin -> COMPONENT_ORDER[origin],
                )
            ],
            "source_index" => nothing,
            "source_row_count" => nothing,
        ))
    end
    return records
end

function wiring_records(
    blocks,
    source::CoreMGKPlan,
    row_positions,
    row_ids,
    progress_every::Int,
    progress_callback,
    total_pairs::Int,
)
    result = Any[]
    nonzero_coefficients = 0
    completed_pairs = 0
    for block in blocks
        dimension = length(block.manifest.entries)
        for j in 1:dimension, k in j:dimension
            wiring = core_mgk_pair(source, block.role, j, k)
            components = Any[]
            for component in wiring.component_records
                coefficients = sort!(
                    copy(component.coefficients);
                    by=record -> row_positions[record.row],
                )
                nonzero_coefficients += length(coefficients)
                push!(
                    components,
                    object(
                        "component" => string(component.component),
                        "status" => string(component.status),
                        "zero_reason" => isnothing(component.zero_reason) ?
                            nothing : string(component.zero_reason),
                        "coefficients" => Any[
                            object(
                                "row_id" => row_ids[record.row],
                                "coefficient" => gaussian(record.coefficient),
                            )
                            for record in coefficients
                        ],
                    ),
                )
            end
            push!(
                result,
                object(
                    "block_id" => block.id,
                    "j_index" => BigInt(j),
                    "k_index" => BigInt(k),
                    "j_entry_id" => block.entry_ids[j],
                    "k_entry_id" => block.entry_ids[k],
                    "component_records" => components,
                ),
            )
            completed_pairs += 1
            report_pair_progress(
                progress_callback,
                progress_every,
                :wiring_records,
                completed_pairs,
                total_pairs,
            )
        end
    end
    return result, nonzero_coefficients
end

function section_sha256(name::String, value)
    io = IOBuffer()
    write(io, "AISECTION1\n")
    write(io, encode_value(name))
    write(io, encode_value(value))
    write(io, "\n")
    return bytes2hex(sha256(take!(io)))
end

"""
Build the complete gamma-independent native Square `core_mgk` inventory and
its hash envelope. No conic model or optimizer is constructed.
"""
function build_square_core_inventory(
    problem::GapProblem;
    progress_every::Integer=0,
    progress_callback=nothing,
)
    progress_every >= 0 ||
        throw(ArgumentError("progress_every must be nonnegative"))
    progress_every <= typemax(Int) ||
        throw(ArgumentError("progress_every is too large"))
    progress_interval = Int(progress_every)
    progress_interval == 0 || progress_callback isa Function ||
        throw(ArgumentError(
            "a progress callback is required when progress_every is nonzero",
        ))
    source = core_mgk_plan(problem)
    validate_basis_manifest(source.positive_basis, problem, :positive) ||
        error("positive structured manifest failed recomputation")
    validate_basis_manifest(source.gap_basis, problem, :gap) ||
        error("gap structured manifest failed recomputation")

    hamiltonian_terms = canonical_hamiltonian(source)
    exact_config = exact_model_config(problem)
    normalization = normalization_spec()
    sites = model_site_table(problem)
    positive_selector = selector_spec(source.positive_basis)
    gap_selector = selector_spec(source.gap_basis)
    action, action_reference = unrestricted_action(
        [source.positive_basis, source.gap_basis],
    )
    actions = Any[action]
    action_set_digest = action_set_sha256(actions)
    positive_ids = entry_id.(source.positive_basis.entries)
    gap_ids = entry_id.(source.gap_basis.entries)
    block_selector = block_selector_spec(
        positive_ids,
        gap_ids,
        action_reference,
    )
    convention = coefficient_convention()
    relaxation = relaxation_id(
        exact_config,
        normalization,
        sites,
        hamiltonian_terms,
        problem,
        positive_selector,
        gap_selector,
        action_set_digest,
        block_selector,
        convention,
    )
    blocks = [
        block_record(
            relaxation,
            :positive,
            source.positive_basis,
            positive_selector,
            block_selector,
            action_reference,
            convention,
        ),
        block_record(
            relaxation,
            :gap,
            source.gap_basis,
            gap_selector,
            block_selector,
            action_reference,
            convention,
        ),
    ]
    ordered_block_ids = [block.id for block in blocks]
    block_order_digest = domain_sha256("AIBLOCKORDER1", ordered_block_ids)

    total_pairs = sum(
        div(
            length(block.manifest.entries) *
            (length(block.manifest.entries) + 1),
            2,
        )
        for block in blocks
    )
    rows, origins, component_record_count = collect_rows(
        blocks,
        source,
        progress_interval,
        progress_callback,
        total_pairs,
    )
    row_positions = Dict(row => index for (index, row) in enumerate(rows))
    row_ids = Dict(row => row_id(row) for row in rows)
    row_records = scalar_row_records(rows, origins, row_ids)
    wiring, nonzero_coefficient_count =
        wiring_records(
            blocks,
            source,
            row_positions,
            row_ids,
            progress_interval,
            progress_callback,
            total_pairs,
        )
    expected_pairs = total_pairs
    expected_components = div(
        length(source.positive_basis.entries) *
        (length(source.positive_basis.entries) + 1),
        2,
    ) + 3div(
        length(source.gap_basis.entries) *
        (length(source.gap_basis.entries) + 1),
        2,
    )
    length(wiring) == expected_pairs ||
        error("Square core wiring pair coverage is incomplete")
    component_record_count == expected_components ||
        error("Square core component coverage is incomplete")

    setup = object(
        "model_id" => "square-j1-j2",
        "model_version" => BigInt(1),
        "exact_model_config" => exact_config,
        "normalization" => normalization,
        "model_site_table" => sites,
        "degree" => BigInt(problem.d),
        "positive_basis_selector" => positive_selector,
        "gap_basis_selector" => gap_selector,
        "action_set_sha256" => action_set_digest,
        "block_selector" => block_selector,
        "coefficient_convention" => convention,
        "source_mapping_mode" => "native_null",
        "source_mapping_spec_sha256" => nothing,
    )
    hamiltonian = object("terms" => hamiltonian_terms)
    block_records = Any[block.record for block in blocks]
    completeness = object(
        "block_basis_completeness" => Any[
            object(
                "block_id" => block.id,
                "finite_formal_basis_complete" => block.manifest.is_complete,
                "hierarchy_exhaustive" => false,
            )
            for block in blocks
        ],
        "wiring_pair_coverage_complete" => true,
        "component_coverage_complete" => true,
    )
    source_anchors = Any[
        object(
            "anchor_kind" => "structured_basis_manifest",
            "block_id" => block.id,
            "anchor_schema_id" => "structured-basis-manifest-v1",
            "anchor_schema_version" => BigInt(1),
            "sha256" => block.manifest.sha256,
        )
        for block in blocks
    ]
    sort!(source_anchors; lt=encoded_less)

    inventory = object(
        "schema_id" => "shared-core-mgk-inventory",
        "schema_version" => BigInt(1),
        "math_scope" => "core_mgk",
        "relaxation_id" => relaxation,
        "setup" => setup,
        "actions" => actions,
        "hamiltonian" => hamiltonian,
        "ordered_block_ids" => ordered_block_ids,
        "block_order_sha256" => block_order_digest,
        "blocks" => block_records,
        "rows" => row_records,
        "wiring" => wiring,
        "completeness" => completeness,
        "source_anchors" => source_anchors,
    )
    math_bytes = canonical_framed_bytes("AICORE1", inventory)
    math_digest = bytes2hex(sha256(math_bytes))
    section_values = [
        "setup" => object(
            "schema_id" => "shared-core-mgk-inventory",
            "schema_version" => BigInt(1),
            "math_scope" => "core_mgk",
            "relaxation_id" => relaxation,
            "setup" => setup,
        ),
        "actions" => actions,
        "hamiltonian" => hamiltonian,
        "blocks" => object(
            "ordered_block_ids" => ordered_block_ids,
            "block_order_sha256" => block_order_digest,
            "blocks" => block_records,
        ),
        "rows" => row_records,
        "wiring" => wiring,
        "completeness" => completeness,
        "source_anchors" => source_anchors,
    ]
    section_hashes = Any[
        object(
            "section_name" => name,
            "section_sha256" => section_sha256(name, value),
        )
        for (name, value) in section_values
    ]
    envelope = object(
        "envelope_schema_id" => "shared-core-mgk-envelope",
        "envelope_schema_version" => BigInt(1),
        "math_byte_count" => BigInt(length(math_bytes)),
        "math_sha256" => math_digest,
        "section_hashes" => section_hashes,
    )
    envelope_bytes = canonical_framed_bytes("AICOREENV1", envelope)
    return (
        source=source,
        inventory=inventory,
        envelope=envelope,
        math_bytes=math_bytes,
        envelope_bytes=envelope_bytes,
        relaxation_id=relaxation,
        math_sha256=math_digest,
        envelope_sha256=bytes2hex(sha256(envelope_bytes)),
        scalar_rows=length(rows),
        wiring_pairs=length(wiring),
        component_records=component_record_count,
        nonzero_coefficients=nonzero_coefficient_count,
    )
end

function write_new_file(path::AbstractString, bytes::Vector{UInt8})
    ispath(path) && error("refusing to overwrite existing artifact: $path")
    mkpath(dirname(abspath(path)))
    temporary, io = mktemp(dirname(abspath(path)))
    try
        write(io, bytes)
        close(io)
        mv(temporary, path)
    catch
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary)
        rethrow()
    end
    return path
end

function write_square_core_inventory(
    problem::GapProblem,
    math_path::AbstractString,
    envelope_path::AbstractString,
)
    result = build_square_core_inventory(problem)
    write_new_file(math_path, result.math_bytes)
    write_new_file(envelope_path, result.envelope_bytes)
    return result
end

function require_equal(actual, expected, message::AbstractString)
    actual == expected || error(
        "$message: expected $(repr(expected)), got $(repr(actual))",
    )
    return actual
end

function nonzero_gaussian(value::CanonicalObject)
    return !iszero(field(value, "real")) || !iszero(field(value, "imag"))
end

"""
Independently parse and byte-for-byte re-encode a Square core math artifact,
then recompute its envelope, IDs, section hashes, references, and exact
pair/component coverage. This validator never constructs an optimizer.
"""
function validate_square_core_inventory(
    math_path::AbstractString,
    envelope_path::AbstractString,
)
    math_bytes = read(math_path)
    envelope_bytes = read(envelope_path)
    inventory = decode_framed_bytes(math_bytes, "AICORE1")
    inventory isa CanonicalObject || error("core inventory is not an object")
    envelope = decode_framed_bytes(envelope_bytes, "AICOREENV1")
    envelope isa CanonicalObject || error("core envelope is not an object")
    math_digest = bytes2hex(sha256(math_bytes))
    require_equal(
        field(envelope, "envelope_schema_id"),
        "shared-core-mgk-envelope",
        "envelope schema",
    )
    require_equal(field(envelope, "envelope_schema_version"), BigInt(1), "envelope version")
    require_equal(field(envelope, "math_byte_count"), BigInt(length(math_bytes)), "math byte count")
    require_equal(field(envelope, "math_sha256"), math_digest, "math SHA-256")
    require_equal(field(inventory, "schema_id"), "shared-core-mgk-inventory", "inventory schema")
    require_equal(field(inventory, "schema_version"), BigInt(1), "inventory version")
    require_equal(field(inventory, "math_scope"), "core_mgk", "inventory scope")

    setup = field(inventory, "setup")
    actions = field(inventory, "actions")
    hamiltonian = field(inventory, "hamiltonian")
    ordered_block_ids = field(inventory, "ordered_block_ids")
    block_order_digest = field(inventory, "block_order_sha256")
    blocks = field(inventory, "blocks")
    rows = field(inventory, "rows")
    wiring = field(inventory, "wiring")
    completeness = field(inventory, "completeness")
    source_anchors = field(inventory, "source_anchors")
    relaxation = field(inventory, "relaxation_id")

    action_refs = Any[]
    for action in actions
        expected = domain_sha256("AIACT1", without_field(action, "action_sha256"))
        require_equal(field(action, "action_sha256"), expected, "action SHA-256")
        push!(
            action_refs,
            object(
                "action_id" => field(action, "action_id"),
                "action_version" => field(action, "action_version"),
                "action_sha256" => expected,
            ),
        )
    end
    sort!(action_refs; lt=encoded_less)
    require_equal(
        field(setup, "action_set_sha256"),
        domain_sha256("AIACTSET1", action_refs),
        "action-set SHA-256",
    )
    for selector_key in ("positive_basis_selector", "gap_basis_selector")
        selector = field(setup, selector_key)
        require_equal(
            field(selector, "basis_selector_sha256"),
            domain_sha256(
                "AISEL1",
                without_field(selector, "basis_selector_sha256"),
            ),
            "$selector_key SHA-256",
        )
    end
    block_selector = field(setup, "block_selector")
    require_equal(
        field(block_selector, "block_selector_sha256"),
        domain_sha256(
            "AIBLOCKSEL1",
            without_field(block_selector, "block_selector_sha256"),
        ),
        "block-selector SHA-256",
    )
    relaxation_preimage = object(
        "schema_id" => field(inventory, "schema_id"),
        "schema_version" => field(inventory, "schema_version"),
        "math_scope" => field(inventory, "math_scope"),
        "model_id" => field(setup, "model_id"),
        "model_version" => field(setup, "model_version"),
        "exact_model_config" => field(setup, "exact_model_config"),
        "normalization" => field(setup, "normalization"),
        "model_site_table" => field(setup, "model_site_table"),
        "exact_hamiltonian_terms" => field(hamiltonian, "terms"),
        "degree" => field(setup, "degree"),
        "positive_basis_selector_sha256" => field(
            field(setup, "positive_basis_selector"),
            "basis_selector_sha256",
        ),
        "gap_basis_selector_sha256" => field(
            field(setup, "gap_basis_selector"),
            "basis_selector_sha256",
        ),
        "action_set_sha256" => field(setup, "action_set_sha256"),
        "block_selector_sha256" => field(
            block_selector,
            "block_selector_sha256",
        ),
        "coefficient_convention" => field(setup, "coefficient_convention"),
    )
    require_equal(
        relaxation,
        domain_sha256("AIREL1", relaxation_preimage),
        "relaxation ID",
    )

    for term in field(hamiltonian, "terms")
        require_equal(
            field(term, "term_id"),
            "h:" * domain_sha256("AIHT1", field(term, "support_word")),
            "Hamiltonian term ID",
        )
        nonzero_gaussian(field(term, "coefficient")) ||
            error("Hamiltonian contains an exact-zero term")
        iszero(field(field(term, "coefficient"), "imag")) ||
            error("Hamiltonian coefficient is not real")
    end
    require_equal(
        block_order_digest,
        domain_sha256("AIBLOCKORDER1", ordered_block_ids),
        "block-order SHA-256",
    )
    require_equal(length(blocks), length(ordered_block_ids), "block count")

    block_by_id = Dict{String,CanonicalObject}()
    entry_ids_by_block = Dict{String,Vector{String}}()
    expected_pair_count = 0
    expected_component_count = 0
    for (position, block) in enumerate(blocks)
        spec = field(block, "spec")
        block_id = field(spec, "block_id")
        require_equal(block_id, ordered_block_ids[position], "canonical block order")
        haskey(block_by_id, block_id) && error("duplicate block ID")
        block_by_id[block_id] = block
        expected_block_id = "blk:" * domain_sha256(
            "AIBLK1",
            object(
                "relaxation_id" => relaxation,
                "role" => field(spec, "role"),
                "block_label" => field(spec, "block_label"),
                "instance_key" => field(spec, "instance_key"),
                "basis_selector_id" => field(spec, "basis_selector_id"),
                "basis_selector_version" =>
                    field(spec, "basis_selector_version"),
            ),
        )
        require_equal(block_id, expected_block_id, "block ID")
        entries = field(block, "basis_entries")
        dimension = Int(field(spec, "dimension"))
        require_equal(length(entries), dimension, "block dimension")
        ids = String[]
        for (index, entry) in enumerate(entries)
            require_equal(field(entry, "block_id"), block_id, "basis parent block")
            require_equal(field(entry, "matrix_index"), BigInt(index), "basis matrix index")
            entry_payload = object(
                "state_symbols" => field(entry, "state_symbols"),
                "operator_word" => field(entry, "operator_word"),
            )
            expected_entry_id = "be:" * domain_sha256("AIBE1", entry_payload)
            require_equal(field(entry, "entry_id"), expected_entry_id, "basis entry ID")
            push!(ids, expected_entry_id)
        end
        require_equal(field(spec, "ordered_entry_ids"), ids, "block entry order")
        require_equal(
            field(spec, "basis_order_sha256"),
            domain_sha256(
                "AIBO1",
                object(
                    "block_id" => block_id,
                    "ordered_entry_ids" => ids,
                ),
            ),
            "basis-order SHA-256",
        )
        entry_ids_by_block[block_id] = ids
        pair_count = div(dimension * (dimension + 1), 2)
        expected_pair_count += pair_count
        expected_component_count +=
            field(spec, "role") == "positive" ? pair_count : 3pair_count
    end

    row_position = Dict{String,Int}()
    for (position, row) in enumerate(rows)
        payload = object(
            "state_symbol_multiset" => field(row, "state_symbol_multiset"),
        )
        expected_row_id = "row:" * domain_sha256("AIROW1", payload)
        require_equal(field(row, "row_id"), expected_row_id, "scalar row ID")
        haskey(row_position, expected_row_id) && error("duplicate scalar row")
        row_position[expected_row_id] = position
        isempty(field(row, "origins")) && error("scalar row has no origin")
    end

    observed_pairs = Set{Tuple{String,Int,Int}}()
    referenced_rows = Set{String}()
    component_count = 0
    coefficient_count = 0
    for record in wiring
        block_id = field(record, "block_id")
        block = get(block_by_id, block_id, nothing)
        isnothing(block) && error("wiring references an unknown block")
        spec = field(block, "spec")
        j = Int(field(record, "j_index"))
        k = Int(field(record, "k_index"))
        ids = entry_ids_by_block[block_id]
        1 <= j <= k <= length(ids) || error("wiring matrix pair is out of range")
        require_equal(field(record, "j_entry_id"), ids[j], "wiring j entry")
        require_equal(field(record, "k_entry_id"), ids[k], "wiring k entry")
        pair_key = (block_id, j, k)
        pair_key in observed_pairs && error("duplicate wiring pair")
        push!(observed_pairs, pair_key)
        components = field(record, "component_records")
        expected_names = field(spec, "role") == "positive" ?
            ["M"] : ["K", "G_moment", "G_product"]
        require_equal(
            [field(component, "component") for component in components],
            expected_names,
            "component coverage/order",
        )
        component_count += length(components)
        for component in components
            coefficients = field(component, "coefficients")
            status = field(component, "status")
            if status == "computed_nonzero"
                isempty(coefficients) &&
                    error("nonzero component has no coefficients")
                isnothing(field(component, "zero_reason")) ||
                    error("nonzero component has a zero reason")
            elseif status == "computed_exact_zero"
                isempty(coefficients) ||
                    error("exact-zero component has coefficients")
                isnothing(field(component, "zero_reason")) &&
                    error("exact-zero component lacks a reason")
            else
                error("unknown component status")
            end
            previous_position = 0
            for coefficient in coefficients
                row_id_value = field(coefficient, "row_id")
                position = get(row_position, row_id_value, 0)
                position > previous_position ||
                    error("component row coefficients are not canonical")
                previous_position = position
                nonzero_gaussian(field(coefficient, "coefficient")) ||
                    error("stored coefficient is exact zero")
                j == k &&
                    !iszero(field(field(coefficient, "coefficient"), "imag")) &&
                    error("Hermitian diagonal coefficient is not real")
                push!(referenced_rows, row_id_value)
                coefficient_count += 1
            end
        end
    end
    require_equal(length(wiring), expected_pair_count, "wiring pair coverage")
    require_equal(length(observed_pairs), expected_pair_count, "unique pair coverage")
    require_equal(component_count, expected_component_count, "component coverage")
    require_equal(referenced_rows, Set(keys(row_position)), "referenced scalar rows")
    require_equal(
        field(completeness, "wiring_pair_coverage_complete"),
        true,
        "stored pair completeness",
    )
    require_equal(
        field(completeness, "component_coverage_complete"),
        true,
        "stored component completeness",
    )

    sections = [
        "setup" => object(
            "schema_id" => field(inventory, "schema_id"),
            "schema_version" => field(inventory, "schema_version"),
            "math_scope" => field(inventory, "math_scope"),
            "relaxation_id" => relaxation,
            "setup" => setup,
        ),
        "actions" => actions,
        "hamiltonian" => hamiltonian,
        "blocks" => object(
            "ordered_block_ids" => ordered_block_ids,
            "block_order_sha256" => block_order_digest,
            "blocks" => blocks,
        ),
        "rows" => rows,
        "wiring" => wiring,
        "completeness" => completeness,
        "source_anchors" => source_anchors,
    ]
    stored_section_hashes = field(envelope, "section_hashes")
    require_equal(length(stored_section_hashes), length(sections), "section count")
    for (index, (name, value)) in enumerate(sections)
        stored = stored_section_hashes[index]
        require_equal(field(stored, "section_name"), name, "section order")
        require_equal(
            field(stored, "section_sha256"),
            section_sha256(name, value),
            "$name section SHA-256",
        )
    end
    return (
        canonical_roundtrip=true,
        envelope_valid=true,
        derived_ids_valid=true,
        coverage_valid=true,
        math_byte_count=length(math_bytes),
        math_sha256=math_digest,
        envelope_sha256=bytes2hex(sha256(envelope_bytes)),
        blocks=length(blocks),
        scalar_rows=length(rows),
        wiring_pairs=length(wiring),
        component_records=component_count,
        nonzero_coefficients=coefficient_count,
        optimizer_invoked=false,
    )
end

end
