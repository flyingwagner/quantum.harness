# Square J1-J2 canonical core inventory

`src/SquareCoreInventory.jl` serializes the exact, gamma-independent Square
`M/G/K` tensor into the reviewed `shared-core-mgk-inventory` v1 wire format.
It includes:

- exact model, normalization, sites, and accumulated Hamiltonian terms;
- an executable unrestricted identity action;
- hashed positive/gap basis and flat two-block selectors;
- canonical block, basis-entry, and scalar-row content IDs;
- every upper-triangle pair and every required `M`, `K`, `G_moment`, and
  `G_product` component, including explicit exact-zero records;
- row origins, finite-basis versus hierarchy completeness, and exact coverage;
- eight canonical section hashes and a separate hash envelope.

Numeric gamma is not part of the math bytes. It belongs to a later conic
evaluation that derives `A_gamma = K - gamma(G_moment + G_product)`.
This was tested by independent full reconstructions at gamma `1/10` and `1/5`:
both math files and both envelopes are byte-for-byte identical, not merely
equal by digest.

## Emit and validate

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/emit_square_core_inventory.jl \
  .bohr-handoff/square-l1-d2-g0p5.core.aicore \
  .bohr-handoff/square-l1-d2-g0p5.core.aicoreenv

julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/validate_square_core_inventory.jl \
  .bohr-handoff/square-l1-d2-g0p5.core.aicore \
  .bohr-handoff/square-l1-d2-g0p5.core.aicoreenv
```

An optional third emitter argument selects a different exact gamma solely to
test this separation, for example `1/5`.

For `L=1,d=2,g=1/2`, emission took about 7½ minutes. The validator independently
parsed and byte-identically re-encoded the 161,886,794-byte artifact, then
recomputed its envelope, section hashes, all derived IDs, references,
Hermitian diagonal checks, and pair/component coverage in about two minutes.

| quantity | value |
|---|---:|
| relaxation ID | `f714640f678a09f62c55e123bbecd25f309018f335e1b00112d6c41abe8df58a` |
| math SHA-256 | `1d466cea7256a71d03d50f84d58f2d088d69360bfeb5b18161108d10c16e8549` |
| envelope SHA-256 | `8944d5e41afd77a187678a2a3981c898ec9cce4772592fd28e3eb623ccec66f5` |
| blocks | 2 |
| scalar rows | 74,602 |
| pair records | 247,484 |
| component records | 247,540 |
| nonzero coefficients | 247,824 |
| optimizer invoked | false |

The first serialization attempt is preserved under
`.bohr-handoff/failed-square-core-inventory-order-v1/`. Its independent
validation correctly failed because row payloads used display order while
their content IDs used canonical byte order. The corrected emitter uses the
same canonical payload object for both the stored record and ID preimage.

## Remaining boundary

This closes the native shared `core_mgk` inventory and envelope, but not legacy
compatibility Gate C: no complete pinned SpectralGap source-event trace or
frozen legacy mapping exists yet. The distinct solver-free evaluation envelope
is now implemented in `SQUARE_STATUS_ENVELOPE.md`; it binds gamma,
stationarity/cone rendering, MOF, source, and environment with hard
`status=unsolved` semantics. Runtime solver evidence still requires a separate
immutable result envelope and independent replay.
