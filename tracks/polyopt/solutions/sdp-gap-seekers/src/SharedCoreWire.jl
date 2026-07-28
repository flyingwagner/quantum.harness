module SharedCoreWire

using SHA
using Unicode
using ..SquareJ1J2Prototype: PauliWord
using ..GenericGapModel: StateMonomial
using ..CoreMGK: ScalarMoment

export CanonicalObject,
    canonical_math_bytes,
    canonical_framed_bytes,
    canonical_sha256,
    canonical_state_monomial,
    canonical_scalar_row,
    canonical_support_word,
    decode_framed_bytes,
    decode_math_bytes,
    domain_sha256,
    encode_value,
    entry_id,
    row_id,
    term_id

const BigRational = Rational{BigInt}

"""An object with an explicit field list; encoding sorts keys canonically."""
struct CanonicalObject
    fields::Vector{Pair{String,Any}}

    function CanonicalObject(fields::Vector{Pair{String,Any}})
        keys = first.(fields)
        length(unique(keys)) == length(keys) ||
            throw(ArgumentError("canonical object contains a duplicate key"))
        new(fields)
    end
end

CanonicalObject(fields::Pair...) =
    CanonicalObject(Pair{String,Any}[String(first(pair)) => last(pair) for pair in fields])

function valid_string(value::AbstractString; key::Bool=false)
    text = String(value)
    isvalid(text) || throw(ArgumentError("canonical string is not valid UTF-8"))
    occursin('\0', text) &&
        throw(ArgumentError("canonical string contains forbidden U+0000"))
    Unicode.normalize(text; compose=true, stable=true) == text ||
        throw(ArgumentError("canonical string is not NFC"))
    if key
        all(isascii, text) ||
            throw(ArgumentError("canonical object key is not ASCII"))
        occursin(r"^[a-z][a-z0-9_.-]*$", text) ||
            throw(ArgumentError("invalid canonical object key '$text'"))
    end
    return text
end

function rational_parts(value::Integer)
    return BigInt(value), BigInt(1)
end

function rational_parts(value::Rational)
    return BigInt(numerator(value)), BigInt(denominator(value))
end

function rational_parts(value)
    throw(ArgumentError("canonical exact rational cannot encode $(typeof(value))"))
end

signed_decimal(value::Integer) = string(BigInt(value))

function encode_rational_body(io::IO, value)
    numerator_value, denominator_value = rational_parts(value)
    write(
        io,
        signed_decimal(numerator_value),
        "/",
        string(denominator_value),
    )
    return io
end

function sorted_fields(object::CanonicalObject)
    fields = [
        valid_string(first(pair); key=true) => last(pair)
        for pair in object.fields
    ]
    sort!(fields; by=first)
    return fields
end

function sorted_fields(object::AbstractDict)
    return sorted_fields(CanonicalObject(
        Pair{String,Any}[String(key) => value for (key, value) in object],
    ))
end

function byte_vector_less(left_bytes, right_bytes)
    for index in 1:min(length(left_bytes), length(right_bytes))
        left_bytes[index] == right_bytes[index] || return left_bytes[index] < right_bytes[index]
    end
    return length(left_bytes) < length(right_bytes)
end

byte_less(left::AbstractString, right::AbstractString) =
    byte_vector_less(codeunits(left), codeunits(right))

function write_value(io::IO, value)
    if value === nothing
        write(io, "N;")
    elseif value isa Bool
        write(io, value ? "B1;" : "B0;")
    elseif value isa Integer
        write(io, "I", signed_decimal(value), ";")
    elseif value isa Rational
        write(io, "Q")
        encode_rational_body(io, value)
        write(io, ";")
    elseif value isa Complex
        write(io, "C")
        encode_rational_body(io, real(value))
        write(io, ",")
        encode_rational_body(io, imag(value))
        write(io, ";")
    elseif value isa AbstractFloat
        throw(ArgumentError("floats are forbidden in canonical math bytes"))
    elseif value isa AbstractString
        text = valid_string(value)
        bytes = Vector{UInt8}(codeunits(text))
        write(io, "S", string(length(bytes)), ":")
        write(io, bytes)
    elseif value isa AbstractVector
        write(io, "L", string(length(value)), ":")
        for item in value
            write_value(io, item)
        end
    elseif value isa Tuple
        throw(ArgumentError("tuples are not a canonical wire type; use a list"))
    elseif value isa NamedTuple
        write_value(
            io,
            CanonicalObject(
                Pair{String,Any}[
                    string(key) => getproperty(value, key)
                    for key in keys(value)
                ],
            ),
        )
    elseif value isa CanonicalObject || value isa AbstractDict
        fields = sorted_fields(value)
        write(io, "O", string(length(fields)), ":")
        for (key, field_value) in fields
            write_value(io, key)
            write_value(io, field_value)
        end
    else
        throw(ArgumentError("unsupported canonical value type $(typeof(value))"))
    end
    return io
end

function encode_value(value)
    io = IOBuffer()
    write_value(io, value)
    return take!(io)
end

function canonical_framed_bytes(header::AbstractString, value)
    header_text = valid_string(header)
    occursin(r"^[A-Z][A-Z0-9]*$", header_text) ||
        throw(ArgumentError("canonical frame header is not an uppercase identifier"))
    io = IOBuffer()
    write(io, header_text, "\n")
    write_value(io, value)
    write(io, "\n")
    return take!(io)
end

canonical_math_bytes(value) = canonical_framed_bytes("AICORE1", value)

canonical_sha256(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))
canonical_sha256(value) = canonical_sha256(canonical_math_bytes(value))

mutable struct Parser
    bytes::Vector{UInt8}
    index::Int
end

at_end(parser::Parser) = parser.index > length(parser.bytes)

function take_byte!(parser::Parser)
    at_end(parser) && error("unexpected end of canonical bytes")
    byte = parser.bytes[parser.index]
    parser.index += 1
    return byte
end

function expect_byte!(parser::Parser, expected::UInt8)
    actual = take_byte!(parser)
    actual == expected ||
        error("expected byte $(Char(expected)), got $(Char(actual))")
    return nothing
end

function take_until!(parser::Parser, delimiter::UInt8)
    start = parser.index
    while !at_end(parser) && parser.bytes[parser.index] != delimiter
        parser.index += 1
    end
    at_end(parser) && error("unterminated canonical scalar")
    bytes = parser.bytes[start:(parser.index - 1)]
    parser.index += 1
    return String(bytes)
end

function parse_uint(text::AbstractString)
    occursin(r"^(0|[1-9][0-9]*)$", text) ||
        error("noncanonical unsigned decimal '$text'")
    return parse(BigInt, text)
end

function parse_signed(text::AbstractString)
    occursin(r"^(0|-?[1-9][0-9]*)$", text) ||
        error("noncanonical signed decimal '$text'")
    return parse(BigInt, text)
end

function parse_rational_text(text::AbstractString)
    pieces = split(text, "/"; keepempty=true)
    length(pieces) == 2 || error("malformed rational body '$text'")
    numerator_value = parse_signed(pieces[1])
    denominator_value = parse_uint(pieces[2])
    denominator_value > 0 || error("rational denominator is not positive")
    gcd(abs(numerator_value), denominator_value) == 1 ||
        error("rational body is not reduced")
    iszero(numerator_value) && denominator_value != 1 &&
        error("zero rational must be 0/1")
    return numerator_value // denominator_value
end

function parse_count!(parser::Parser)
    value = parse_uint(take_until!(parser, UInt8(':')))
    value <= typemax(Int) || error("canonical container count exceeds Int")
    return Int(value)
end

function parse_string!(parser::Parser)
    length_value = parse_count!(parser)
    last_index = parser.index + length_value - 1
    last_index <= length(parser.bytes) || error("truncated canonical string")
    bytes = parser.bytes[parser.index:last_index]
    parser.index = last_index + 1
    isvalid(String, bytes) || error("canonical string contains invalid UTF-8")
    return valid_string(String(bytes))
end

function parse_value!(parser::Parser)
    tag = take_byte!(parser)
    if tag == UInt8('N')
        expect_byte!(parser, UInt8(';'))
        return nothing
    elseif tag == UInt8('B')
        digit = take_byte!(parser)
        digit in (UInt8('0'), UInt8('1')) || error("invalid canonical bool")
        expect_byte!(parser, UInt8(';'))
        return digit == UInt8('1')
    elseif tag == UInt8('I')
        return parse_signed(take_until!(parser, UInt8(';')))
    elseif tag == UInt8('Q')
        return parse_rational_text(take_until!(parser, UInt8(';')))
    elseif tag == UInt8('C')
        body = take_until!(parser, UInt8(';'))
        pieces = split(body, ","; keepempty=true)
        length(pieces) == 2 || error("malformed complex rational")
        return complex(
            parse_rational_text(pieces[1]),
            parse_rational_text(pieces[2]),
        )
    elseif tag == UInt8('S')
        return parse_string!(parser)
    elseif tag == UInt8('L')
        count = parse_count!(parser)
        return Any[parse_value!(parser) for _ in 1:count]
    elseif tag == UInt8('O')
        count = parse_count!(parser)
        fields = Pair{String,Any}[]
        previous = nothing
        for _ in 1:count
            take_byte!(parser) == UInt8('S') ||
                error("canonical object key is not a string")
            key = parse_string!(parser)
            valid_string(key; key=true)
            if !isnothing(previous)
                byte_less(previous, key) ||
                    error("canonical object keys are not strictly increasing")
            end
            push!(fields, key => parse_value!(parser))
            previous = key
        end
        return CanonicalObject(fields)
    end
    error("unknown canonical type tag $(Char(tag))")
end

function decode_framed_bytes(
    input::AbstractVector{UInt8},
    header_text::AbstractString,
)
    bytes = Vector{UInt8}(input)
    header = Vector{UInt8}(codeunits(valid_string(header_text) * "\n"))
    length(bytes) >= length(header) + 1 || error("truncated canonical math file")
    bytes[1:length(header)] == header || error("canonical math header mismatch")
    last(bytes) == UInt8('\n') || error("canonical math file lacks final newline")
    parser = Parser(bytes, length(header) + 1)
    value = parse_value!(parser)
    parser.index == length(bytes) ||
        error("canonical math file has trailing or missing bytes")
    canonical_framed_bytes(header_text, value) == bytes ||
        error("decoded canonical value does not re-encode byte-identically")
    return value
end

decode_math_bytes(input::AbstractVector{UInt8}) =
    decode_framed_bytes(input, "AICORE1")

function preimage_sha256(prefix::AbstractString, value)
    io = IOBuffer()
    write(io, prefix, "\n")
    write(io, encode_value(value))
    write(io, "\n")
    return bytes2hex(sha256(take!(io)))
end

domain_sha256(prefix::AbstractString, value) =
    preimage_sha256(prefix, value)

function canonical_support_word(
    word::PauliWord;
    site_namespace::AbstractString="model",
)
    namespace = valid_string(site_namespace)
    axes = ("X", "Y", "Z")
    factor_records = [
        (
            site=site,
            axis=Int(axis),
            value=CanonicalObject(
                "axis" => axes[Int(axis)],
                "site_id" => BigInt(site),
                "site_namespace" => namespace,
            ),
        ) for (site, axis) in word.ops
    ]
    sort!(factor_records; by=record -> (record.site, record.axis))
    factors = Any[record.value for record in factor_records]
    return CanonicalObject("ordered_factors" => factors)
end

function canonical_state_monomial(
    monomial::StateMonomial;
    site_namespace::AbstractString="model",
)
    symbols = [
        canonical_support_word(word; site_namespace=site_namespace)
        for word in monomial.state_symbols
    ]
    sort!(symbols; lt=(left, right) ->
        byte_vector_less(encode_value(left), encode_value(right)))
    return CanonicalObject(
        "operator_word" => canonical_support_word(
            monomial.operator_word;
            site_namespace=site_namespace,
        ),
        "state_symbols" => symbols,
    )
end

function canonical_scalar_row(
    row::ScalarMoment;
    site_namespace::AbstractString="model",
)
    symbols = [
        canonical_support_word(word; site_namespace=site_namespace)
        for word in row.state_symbol_multiset
    ]
    sort!(symbols; lt=(left, right) ->
        byte_vector_less(encode_value(left), encode_value(right)))
    return CanonicalObject("state_symbol_multiset" => symbols)
end

term_id(word::PauliWord; site_namespace::AbstractString="model") =
    "h:" * preimage_sha256(
        "AIHT1",
        canonical_support_word(word; site_namespace=site_namespace),
    )

entry_id(monomial::StateMonomial; site_namespace::AbstractString="model") =
    "be:" * preimage_sha256(
        "AIBE1",
        canonical_state_monomial(monomial; site_namespace=site_namespace),
    )

row_id(row::ScalarMoment; site_namespace::AbstractString="model") =
    "row:" * preimage_sha256(
        "AIROW1",
        canonical_scalar_row(row; site_namespace=site_namespace),
    )

end
