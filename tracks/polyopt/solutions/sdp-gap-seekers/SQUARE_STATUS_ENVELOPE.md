# Square J1-J2 source-gated status envelope

`src/SquareStatusEnvelope.jl` writes a canonical solver-free status record for
the locked Square source point:

- `H = ¼Σ_J1(XX+YY+ZZ) + (g/4)Σ_J2(XX+YY+ZZ)`;
- `L=1`, `d=2`, `g=1/2`, and `gamma=1/10`;
- unrestricted infinite-volume KMS ground states;
- flat structured positive/gap bases, with no symmetry quotient.

The emitter fails closed unless the tracked worktree is clean. It binds the
full Git commit and tree, Project/Manifest hashes, physical semantics, exact
source-plan and basis hashes, canonical core hashes and counts, and rendered
MOF hash and inventory. It also:

1. independently parses and validates the complete 161 MB canonical core;
2. rebuilds that core from the current source and requires both artifact
   hashes to match;
3. rebuilds the exact conic plan and compares every rendered affine and PSD
   coefficient, variable name/order, cone dimension, and objective sense.

No optimizer is constructed or invoked.

## Emit or replay

From the repository root, using the isolated Challenge 88 Git database:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/emit_square_status_envelope.jl \
  .bohr-handoff/square-l1-d2-g0p5.core.aicore \
  .bohr-handoff/square-l1-d2-g0p5.core.aicoreenv \
  .bohr-handoff/square-l1-d2-g0p5-gamma0p1.mof.json.gz \
  . .bohr-handoff/local-git \
  julia-env/Project.toml julia-env/Manifest.toml \
  .bohr-handoff/square-l1-d2-g0p5-gamma0p1.aisqstatus
```

For an ordinary checkout, replace `.bohr-handoff/local-git` with `.git`.
Output is write-once: a completed same-directory temporary file is published
by an atomic hard link, and an existing or concurrently created destination is
rejected without modification.

## Measured source-bound result

The completed solver-free audit used:

- source commit `e5c5249047cd0df7c6387483ee2e745eabed2514`;
- source tree `3a13ccb1089d0e5dfc0763c5cc463d4bf23e9bdc`;
- core math SHA
  `1d466cea7256a71d03d50f84d58f2d088d69360bfeb5b18161108d10c16e8549`;
- core envelope SHA
  `8944d5e41afd77a187678a2a3981c898ec9cce4772592fd28e3eb623ccec66f5`;
- MOF SHA
  `8cc83c7ed497a940823b9b57433c3b6b50f2aaa20ae9330746fa0ba0ede60685`.

The complete core source rebuild matched, every MOF coefficient replayed, and
no optimizer was invoked. The 1,885-byte output is:

```text
.bohr-handoff/square-l1-d2-g0p5-gamma0p1.aisqstatus
SHA-256 6f8d97f69a3e1f3aabc9b41fed5826573f81e25bf62a108922ecf944c57a3989
```

Wall time was 1154.82s. An independent decode/canonical-reencode was
byte-identical and recovered `status=unsolved`,
`optimizer_invoked=false`, and the same source commit and claim boundary.
The first full audit at source `46b7bd8` passed every mathematical gate but
failed only because Julia does not support `open(path,"x")`; it created no
output. The atomic hard-link publisher fixed that isolated failure.

## Interpretation boundary

The hard-coded status is `unsolved`, and the claim boundary is
`source-gated conic model; no Square bulk-gap bound`. The envelope establishes
that the committed source, exact core, and Float64 MOF are the same declared
mathematical relaxation. It does not establish feasibility or infeasibility.

A later solve must use a separate immutable run/result envelope and preserve
the same source identity. An upper-bound statement `Δ_bulk ≤ gamma` requires
independently replayed strict infeasibility evidence; a raw solver status is
insufficient. Any symmetry-restricted run requires a different, explicitly
hashed state-class declaration.
