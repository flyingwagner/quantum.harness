# Kagome γ=1.272 conditioning audit

The supplied `N=13, d=3, lso=5`, sign-symmetric Kagome ray remains rejected.
Nothing in this audit changes the `1e-12` independent replay tolerance or
supports calling γ=1.272 infeasible.

## Reproduce the solver-free diagnosis

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/diagnose_ray_conditioning.jl \
  /path/to/audit.mof.json.gz \
  /path/to/audit.variables.tsv
```

The diagnostic evaluates the exact binary Float64 model coefficients and ray
coordinates at 256-bit precision. It therefore distinguishes a defective
exported ray from error introduced only by Float64 summation.

For the supplied γ=1.272 pair it finds:

| quantity | measured value |
|---|---:|
| variables / affine equalities | 54,944 / 15,671 |
| ray max norm | `8.589574013609442e16` |
| largest 256-bit equality residual | `5.682359651454748e6` |
| residual / ray max norm | `6.615414969882717e-11` |
| largest Float64-vs-256-bit summation difference | `1.195931051382795e2` |
| exact duplicate affine equalities | 4,887 |

The 256-bit residual is essentially the same as the exported residual, while
the summation difference is about 47,000 times smaller. The failure is in the
represented solver ray, not in the verifier's summation order.

Rows 1, 6,494, 16, and 8,277 carry the largest defects. Their dominant terms
come from the 17×17 second gap block, including variable 42,488 at
`7.343587887752901e16`. Rows 23 and 1,613 are instead dominated by the 104×104
second positive block. The leading γ-dependent coefficient in row 1 is the
binary Float64 representation of `−1.272`; this is consistent with the locked
γ convention and is not evidence for a sign or normalization mismatch.

The four certificate blocks span roughly `8.1e11` to `8.6e16`. The five
`posepsd9!` strengthening blocks span only about `1.6e-4` to `1.8e-3`, leaving
about 20 decimal orders between active conic blocks in the same model. This,
together with the 4,887 duplicate equality rows, is a concrete conditioning
mechanism compatible with Mosek's `SLOW_PROGRESS` / unknown status.

## Safe equivalent transformation

`deduplicate_mof_equalities.jl` removes only bit-for-bit identical affine
equalities, including the affine offset. It uses no tolerance and preserves
the feasible set, objective, variables, and cones exactly.

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/deduplicate_mof_equalities.jl \
  /path/to/audit.mof.json.gz \
  /path/to/audit.deduplicated.mof.json.gz
```

On the supplied model this reduces the equality inventory from 15,671 to
10,784 while retaining all 54,944 variables and nine PSD cones. Replaying the
old ray against the reduced model still rejects it with normalized equality
residual `6.615275739340028e-11`. Deduplication is a solver-conditioning A/B
treatment, not post-hoc acceptance of the old ray.

For the xH5 A/B run, `solve_exported_mof.jl` copies either exported model into
the same pinned MosekTools optimizer and writes raw statuses plus a canonical
variable TSV. Its preflight parses and inventories a model without creating an
optimizer or invoking optimization:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/solve_exported_mof.jl \
  --preflight /path/to/audit.deduplicated.mof.json.gz
```

Uniform row normalization also cannot validate this ray: the largest
row-scaled residual is `5.200797932297734e-11`. Per-block positive rescaling is
mathematically equivalent under a corresponding inverse coefficient scaling,
but the appropriate factors require a fresh solver A/B and should not be
chosen from a failed ray alone.

The completed xH5 run is reported in `KAGOME_AB_RESULT.md`. Deduplication
reduced solve time, peak memory, and ray scale, but both new rays were rejected
for normalized equality residual at `1e-12`; the deduplicated residual was
worse. One exactly round-tripped block-equilibration A/B is therefore the next
isolated experiment.
