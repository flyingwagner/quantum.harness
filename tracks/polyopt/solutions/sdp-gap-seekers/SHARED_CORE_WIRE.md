# Canonical shared-core wire

`src/SharedCoreWire.jl` implements the solver-independent `AICORE1` typed byte
grammar used to identify exact mathematical objects before any floating-point
or solver rendering.

The implemented scalar types are null, Boolean, integer, reduced rational,
Gaussian rational, UTF-8 string, list, and canonical object. Strings must be
NFC and may not contain U+0000. Object keys are restricted to the contract's
ASCII identifier grammar and are encoded in byte order. Floats and tuples fail
closed. A decoded file is accepted only when re-encoding returns the identical
bytes.

Adapters cover:

- canonical Pauli support words;
- state-monomial basis entries;
- scalar moment rows;
- domain-separated `h:`, `be:`, and `row:` content IDs.

`SquareCoreInventory.jl` builds the native unrestricted Square block,
coefficient, wiring, completeness, source-anchor, section-hash, and envelope
records over these primitives. Its independent validator replays all derived
IDs and full coverage.

The unit tests reproduce the reviewed 70-byte grammar golden
`b75e0d4c...958975` and 181-byte `H=Z` golden
`4c8f8281...70ad4`, reject noncanonical object/rational encodings, and prove
that scalar state-symbol multisets are order independent.

The native Square artifact is sufficient for a byte-for-byte full-tensor diff
against another conforming native emitter. A legacy-aligned comparison still
requires the separately reviewed SpectralGap source mapping and complete
source-event trace; the native emitter does not fabricate that evidence.
