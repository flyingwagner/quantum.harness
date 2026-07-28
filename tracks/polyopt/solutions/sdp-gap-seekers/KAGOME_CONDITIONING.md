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

Uniform row normalization alone cannot validate this ray: the largest
row-scaled residual is `5.200797932297734e-11`.

The completed xH5 run is reported in `KAGOME_AB_RESULT.md`. Deduplication
reduced solve time, peak memory, and ray scale, but both new rays were rejected
for normalized equality residual at `1e-12`; the deduplicated residual was
worse. One exactly round-tripped block-equilibration A/B is therefore the next
isolated experiment.

## Power-of-two block/row equilibration

`equilibrate_mof_columns.jl` now constructs that isolated experiment from the
completed original A/B ray:

1. every direct PSD matrix receives one uniform positive power-of-two variable
   scale, so `X ⪰ 0` is equivalent to `Z ⪰ 0` under `X=sZ`;
2. variables outside direct PSD blocks may receive independent positive
   power-of-two scales;
3. every scalar equality and both sides receive one positive power-of-two row
   factor;
4. each affine PSD block receives only one common positive factor;
5. the objective receives one positive factor, preserving improvement
   direction.

Arbitrary coordinate-wise scaling inside a direct PSD matrix is explicitly
rejected. A structured coordinate map is allowed only when its triangular
scales factor as `sᵢⱼ=dᵢdⱼ`, because then it is the invertible diagonal
congruence `X=DZD` and preserves PSD membership in both directions.
The first block-only prototype was mathematically equivalent but produced
coefficients as large as `1.24e35`; it was not submitted and is preserved as a
failed conditioning attempt. Composing row/objective normalization reduces the
largest transformed column coefficient to `1.272`.

For the original job `22986235` ray, the final map has nine direct PSD blocks,
variable exponents `[-13,58]`, equality-row exponents `[-59,10]`, and objective
exponent `-41`. It maps the old ray scale
`4.319157242331815e17` to `1.9419280324476962`; the inverse map reproduces the
original ray TSV byte-for-byte.

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/equilibrate_mof_columns.jl \
  ORIGINAL.mof.json.gz EQUILIBRATED.mof.json.gz SCALE.tsv \
  ORIGINAL.variables.tsv

julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/backtransform_mof_ray.jl \
  SCALE.tsv SCALED.variables.tsv ORIGINAL-COORDINATES.variables.tsv
```

The equilibrated MOF SHA-256 is
`830de74414ea3f1f16037376a9dfb3ec38fa01e31e5e80420f34c981a6411a83`.
Its preflight preserves 54,944 variables, 15,671 equalities, and PSD
dimensions `271,104,18,17,1,9,36,84,126`.

The mapped old ray is deliberately not accepted in scaled coordinates: its
small strengthening-block PSD defect becomes order one after balancing,
instead of being hidden by the global `4.3e17` scale. This is precisely why a
fresh solve can be informative. Any returned scaled ray must be inverse-mapped
and replayed against the immutable original MOF at `1e-12`; the transformed
model's own residual is not acceptance evidence.

## Completed uniform-block xH5 result

Job `22988046`, source `57befde`, ran the declared uniform-block/row map with
36 CPU and 128 GiB on `xhacnormalb`. It completed in 4m53s with batch MaxRSS
15,957,768 KiB; Mosek spent 66.56s and again returned `SLOW_PROGRESS`,
primal/dual `UNKNOWN_RESULT_STATUS`.

The back-transformed ray passes PSD and improving-objective checks but fails
the original-model equality check:

- normalized equality residual: `9.914901660490167e-12`;
- normalized PSD violation: between `1.27e-17` and `1.65e-17` across the xH5
  and local LAPACK replays, both safely below tolerance;
- normalized improving objective: `6.452369852106409e-8`;
- original-coordinate scale: `5.808433831135018e34`.

The equality residual is about 3.3 times better than the original A/B ray but
still about 9.9 times the accepted tolerance. The complete result and
independent local replay are documented in
`KAGOME_EQUILIBRATION_RESULT.md`. No infeasibility or gap bound follows.

## Next isolated map: diagonal congruence

`ray_congruence_equilibration` uses the diagonal of the supplied ray to choose
one positive power-of-two row/column factor per PSD matrix index. Stored
triangle coordinates receive `sᵢⱼ=dᵢdⱼ`; the code verifies that factorization
before omitting the factors from the direct PSD constraint.

This matters most for the 104×104 block, whose reference diagonal spans
approximately `1.1e-3` to `4.0e16`. The generated map:

- preserves all 54,944 variables, 15,671 equalities, and nine PSD dimensions;
- uses coordinate exponents `[-12,58]`, equality-row exponents `[-58,10]`,
  and objective exponent `-41`;
- maps reference-ray scale `4.319157242331815e17 → 1.9645621823673163`;
- limits the maximum transformed column coefficient to `1.272`;
- round-trips the reference ray byte-for-byte.

The mapped old ray still rejects, as expected for an inexact old candidate.
A fresh solve is decision-relevant because this map equilibrates within the
large PSD blocks, which the completed uniform-block experiment could not do.

The fresh diagonal-congruence job `22988126` nevertheless stalled and its
back-transformed original-model equality residual was
`7.049624122338534e-11`. PSD and objective checks passed, but the unchanged
`1e-12` verifier rejected it. Full details are in
`KAGOME_CONGRUENCE_RESULT.md`.

Both reference-ray-derived maps have now been tested. The next route is not
another coordinate heuristic: explicitly construct the homogeneous recession
system and normalize its improving objective to one. That turns existence of
an improving ray into finite conic feasibility while preserving the exact
mathematical certificate condition.

That route is implemented by `normalize_mof_recession.jl` and documented in
`NORMALIZED_RECESSION.md`. The solver-free Kagome prototype has 15,672
equalities (the 15,671 homogeneous source rows plus one objective
normalization), retains all nine PSD blocks, and invokes no optimizer.

The corresponding job `22988185` returned a finite unit-improvement point,
but original-model equality residual `9.991205557433147e-11` still failed.
See `KAGOME_RECESSION_RESULT.md`. Its scale is no longer pathological; the
next isolated question is row equilibration of this finite feasibility model.
