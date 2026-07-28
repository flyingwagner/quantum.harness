# Kagome exact affine-projection boundary

The six completed γ=1.272 xH5 formulations all fail the independent affine
equality gate. `affine_peeling_analysis` separates the part of that exact
rational system that can be projected without a general sparse solve from the
irreducibly coupled remainder.

The analysis first removes identical homogeneous rows. It then repeatedly
selects a variable that occurs in exactly one active row and removes that row.
The resulting order is triangular: reverse-order rational corrections solve
every peeled row exactly without changing the final coupled core.

## Reproduce

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/analyze_affine_projection.jl \
  .bohr-handoff/artifacts/certificate-audit-b1a1cad-20260728T102942Z/\
kagome-1.272/audit.mof.json.gz \
  > .bohr-handoff/kagome-affine-projection-structure.tsv
```

The script exits `2` when a coupled core remains. That is a diagnostic
outcome, not an optimizer or certificate failure. No optimizer is invoked.

## Measured structure

| quantity | value |
|---|---:|
| variables | 54,944 |
| original equality rows | 15,671 |
| unique homogeneous rows | 10,784 |
| exact duplicate rows removed | 4,887 |
| peelable unique/original rows | 5,806 / 5,806 |
| coupled unique/original rows | 4,978 / 9,865 |
| columns in coupled core | 12,283 |
| maximum pattern matching | 4,978 |
| full row structural rank | true |
| columns outside coupled core | 42,661 |

The matching-extended output is
`.bohr-handoff/kagome-affine-projection-structure-matching.tsv`, SHA-256
`9d10bc6d7d49be8b3eb3e43fd74ddd54037818c076717074911b6faab5f00491`.
The earlier peel-only output is preserved with SHA-256
`0a6a03e2d5672a966ac55730bf4a5773e593dc91a37189ee385bf50f61bd9ec7`.
Synthetic tests prove complete reverse correction on a triangular fixture and
distinguish a full coupled matching from a structurally rank-deficient
fixture.

## Exact coefficient-rank proof

The complete coupled coefficient matrix can be exported and ranked over prime
fields:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/export_affine_coupled_core.jl \
  .bohr-handoff/artifacts/certificate-audit-b1a1cad-20260728T102942Z/\
kagome-1.272/audit.mof.json.gz \
  .bohr-handoff/kagome-affine-coupled-core.tsv

python tracks/polyopt/solutions/sdp-gap-seekers/scripts/\
rank_affine_matrix_mod.py \
  .bohr-handoff/kagome-affine-coupled-core.tsv \
  1000000007 1000000009 \
  --pivots-output \
  .bohr-handoff/kagome-affine-coupled-core-pivots.tsv
```

The Python script uses `python-flint==0.9.0` / FLINT 3.6.0. Installation and
wheel provenance are recorded in
`.bohr-handoff/dependency-install-log.md`.

The 4,978×12,283 export has 704,704 nonzero exact rational coefficients.
Its SHA-256 is
`18f10deccfcba97005ada3321599fc530c501255c990794cc70a2b0813f73c40`.
FLINT returned rank 4,978 modulo both 1,000,000,007 and 1,000,000,009 in
185.34s. The rank log SHA-256 is
`2fb8f81211e1a02a4c6f7becc36ed897d63bb1e7418c8c4f11f3f3066de28aa8`.

A separate exact modular RREF over 1,000,000,007 selected 4,978 distinct
coefficient-valid pivot columns in 121.91s. The canonical pivot table has
SHA-256
`8c7da3d6f68e5ad54e0caea9e460a21153fafd737ab5745b2845c6042580e1fb`;
its log has SHA-256
`450397fece31c7c952346d683d7aaf446d0d8d56c121e0ca04f312874209a543`.
The selected rational square minor is nonsingular because its reduction
modulo the prime is nonsingular. Pivot output is write-once: an existing
destination is rejected.

One full modular rank is already a proof of full row rank over the rationals:
clearing denominators yields a row-size minor that is nonzero modulo the
prime, hence nonzero over the integers and rationals. The second prime is an
independent implementation check.

The first attempt ranked only the arbitrary 4,978-column matching minor. It
had rank 4,916 over both primes and therefore was a bad pivot selection, not a
rank obstruction. That failed attempt is preserved in
`.bohr-handoff/kagome-affine-coupled-matching-minor-rank.tsv`, SHA-256
`c2cb2895e111d9968182a5e4c8fcf754662a9a005b347ac43eb9be29dfaaf06d`.

## Closed projection result and remaining source boundary

A later exact construction avoided the general 4,978-row solve: one
rationalization cell already satisfied the coupled core exactly, reverse
peeling solved the remaining rows, and an exact six-kernel perturbation
repaired the only indefinite PSD block. Reconstructing the intended MOF
coefficients then required 46 exact private-pivot corrections. The final ray
has zero residual on all 15,671 rows, the positive exact objective
`119089//14841408527`, and rigorous PSD proofs for all nine blocks.

The construction, artifact hashes, and replay command are in
`KAGOME_EXACT_EXPORTED_MODEL_CERTIFICATE.md`. This closes the exported-model
projection/PSD problem. γ=1.272 remains withheld as a physical bound until a
separate source-assembly audit reproduces the Kagome rows, objective, all nine
blocks, geometry, Hamiltonian convention, and sign-symmetric state class.
