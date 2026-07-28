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

## Remaining obstruction and smallest next experiment

Private-pivot correction alone cannot close this Kagome point. All 4,887
duplicate copies belong to representatives in the coupled core, and 4,978
unique equations remain jointly supported on 12,283 columns. Exact
coefficient rank now shows that a rational affine correction exists for every
rationalized residual.

Existence is not yet a certificate: an arbitrary solution may destroy the
very small strengthening-block PSD margins or the improving objective. The
modular pivot-selection step is complete. The smallest next experiment is:

1. form the selected exact rational 4,978×4,978 minor and solve its correction
   system for a normalized supplied ray;
2. reverse the exact leaf-peeling corrections;
3. verify every one of the 15,671 original rational rows exactly and measure
   objective change;
4. test the corrected nine PSD blocks with directed interval factorization.

The correction must retain a rigorously nonnegative cone margin. Full affine
row rank alone does not change γ=1.272 from numerical unknown.
