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

## Obstruction and smallest next experiment

Private-pivot correction alone cannot close this Kagome point. All 4,887
duplicate copies belong to representatives in the coupled core, and 4,978
unique equations remain jointly supported on 12,283 columns. An exact
certificate therefore needs a coefficient-rank-revealing solve of that sparse
rational subsystem, followed by exact residual verification and a rigorous
PSD proof after correction. Maximum bipartite matching covers all 4,978 rows,
so the pattern has full structural row rank; exact coefficient cancellation
is the remaining rank question.

The smallest next experiment is:

1. start from the complete structural matching and check the coefficient rank
   of the 4,978×12,283 matrix modulo at least two large primes;
2. solve the correction equations by modular reconstruction or fraction-free
   sparse elimination;
3. verify every one of the 15,671 original rational rows exactly;
4. test the corrected nine PSD blocks with directed interval factorization.

Failure of modular full row rank would expose exact dependencies that must be
removed first. Success still would not prove PSD membership; the correction
must retain a rigorously nonnegative cone margin.
