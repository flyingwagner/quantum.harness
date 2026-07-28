# Square J1-J2 solver-free conic rendering

`src/SquareGapConic.jl` closes the primal-model fields that are not part of the
gamma-independent `M/G/K` tensor:

- exact `L(1)=1` normalization with right-hand side one;
- exact stationarity rows for selected `q` of degree at most `2d-2`, emitted
  as independent real and imaginary equalities after zero removal and exact
  deduplication;
- exact evaluation of `A_gamma = K - gamma G` only at the render boundary;
- the real PSD embedding
  `[Re(A) -Im(A); Im(A) Re(A)]` for complex Hermitian blocks;
- a feasibility objective with an explicitly empty objective vector;
- deterministic scalar-row order and content IDs.

The source plan remains exact over Gaussian rationals. `render_mof` is the only
Float64 boundary. It rejects coefficients that cannot be recovered as the
same rational at tolerance `1e-14`, writes through the JuMP/MOI file API, and
does not instantiate or invoke an optimizer.

## Reproduce the bounded source gate

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/emit_square_conic_mof.jl \
  .bohr-handoff/square-l1-d2-g0p5-gamma0p1.mof.json.gz

julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/solve_exported_mof.jl \
  --preflight \
  .bohr-handoff/square-l1-d2-g0p5-gamma0p1.mof.json.gz
```

For the locked `L=1, d=2, g=1/2, gamma=1/10` structured source, the measured
render took 6m38s locally and produced:

| quantity | value |
|---|---:|
| scalar variables | 74,602 |
| normalization equalities | 1 |
| stationarity selector entries | 7 |
| nonzero stationarity equalities | 3 |
| exact duplicate equalities removed | 0 |
| complex PSD dimensions | 703, 7 |
| real PSD dimensions | 1406, 14 |
| objective | feasibility |
| optimizer invoked | false |

The 25 MiB ignored MOF artifact has SHA-256
`8cc83c7ed497a940823b9b57433c3b6b50f2aaa20ae9330746fa0ba0ede60685`.
An independent MOI preflight recovered 74,602 variables, four affine
equalities, two PSD constraints with dimensions 1406 and 14, and
`optimization_invoked=false`.

The stronger exact replay command is:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/audit_square_conic_mof.jl \
  .bohr-handoff/square-l1-d2-g0p5-gamma0p1.mof.json.gz
```

It rebuilds the exact source plan, independently parses the MOF, and compares
every variable name/order, right-hand side, affine coefficient, real PSD
coordinate/coefficient, cone dimension, and objective sense. The measured
audit matched all 49 affine and 495,560 PSD coefficients and reported
`exact_coefficient_match=true`, `optimizer_invoked=false`.

## Interpretation boundary

The patch is a local-consistency window for unrestricted infinite-volume KMS
ground states with a flat structured basis and no symmetry quotient; it is not
an open-boundary finite system. The artifact is a source-gated primal
feasibility model, not a solve and not a Square bulk-gap bound.
`SQUARE_STATUS_ENVELOPE.md` documents the implemented write-once evaluation
envelope that rebuilds the core, replays the MOF, and binds
source/environment while reporting `status=unsolved`. A future solve still
requires a separate immutable result envelope and independent replay of any
returned infeasibility evidence.
