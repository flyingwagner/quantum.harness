# Normalized homogeneous recession model

## Mathematical equivalence

For an affine conic program

```text
maximize cᵀx + c₀
subject to A x + a ∈ K,
```

an improving recession direction `d` satisfies

```text
A d ∈ rec(K),   cᵀd > 0.
```

Because the recession constraints are homogeneous, such a direction exists if
and only if one exists with `cᵀd=1`. The minimization normalization is
`cᵀd=-1`. The resulting problem is finite conic feasibility rather than an
unbounded optimization trajectory.

`normalize_recession_problem!` implements this equivalence:

1. scalar equality constants and right-hand sides become zero;
2. vector equality and affine PSD constants become zero;
3. direct PSD cones are retained unchanged;
4. maximization adds `cᵀd=1`, minimization adds `-cᵀd=1`;
5. the new objective sense is feasibility.

Unsupported constraint types and zero objectives fail closed.

## Command

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/normalize_mof_recession.jl \
  ORIGINAL.mof.json.gz NORMALIZED-RECESSION.mof.json.gz
```

This command parses and writes MOF only; it invokes no optimizer.

An isolated positive row-equilibration variant can then be generated without
variable scaling:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/equilibrate_mof_rows.jl \
  NORMALIZED-RECESSION.mof.json.gz ROW-SCALED.mof.json.gz
```

Each scalar equality and both sides receive one positive power-of-two factor.
This preserves the normalized feasible set exactly in real arithmetic.

## Kagome γ=1.272 prototype

The immutable Kagome model becomes:

- 54,944 variables;
- 15,672 scalar equalities: 15,671 homogeneous source rows plus one objective
  normalization;
- nine direct PSD blocks with dimensions
  `271,104,18,17,1,9,36,84,126`;
- feasibility objective.

The normalized-recession MOF SHA-256 is
`ce6ec074f7bee5e0a01b2b29b240758e945317c90709357469606304d28acceb`;
its uncompressed canonical JSON SHA-256 is
`450cf07c215d42a54d4f49822eaa33d56249379d6720cd232a8c45f0d18fe0be`.

Any returned point must still be replayed against the immutable original MOF.
The original verifier checks homogeneous equalities, PSD membership, and
objective direction without trusting the normalization model or solver
status. A floating replay is still not a strict certificate.

Job `22988185` tested this formulation. It produced a finite unit-improvement
point but the immutable-original replay rejected normalized equality residual
`9.991205557433147e-11`. See `KAGOME_RECESSION_RESULT.md`.
