# Kagome γ=1.272 row-equilibrated recession result

## Claim boundary

Positive row equilibration slightly improved the finite unit-ray equality
residual, but the returned point still failed the unchanged `1e-12`
independent replay.

Kagome γ=1.272 remains numerical instability. No tested formulation supports
infeasibility or a certified bulk-gap upper bound.

## Reproducibility

- Slurm job: `22988265`, state `COMPLETED`, exit `0:0`.
- Source commit/tree:
  `7edf6ef0ec7cca1db9f1ddcf2c8f00ac54e76628` /
  `7b1d364413a042c05663f88180b6e1624c2a3d3d`.
- Immutable original MOF SHA-256:
  `2e59c395f3af1cdaf762cd98af3507e3b4ee7b75e315af338d8bef5a8bbb3a53`.
- Row-scaled normalized-recession MOF SHA-256:
  `dce12cc883bfc0292af749d1b6a7e93f21f9c36db8de37d3f93849aaf6737753`.
- Environment: Julia 1.11.5, MosekTools 0.15.10.
- Resources: one `xhacnormalb` node, 36 CPU, 128 GiB, one-hour limit.
- Accounting: 5m01s elapsed, batch MaxRSS 19,933,888 KiB.
- Evidence:
  `.bohr-handoff/xh5-results/kagome-recession-row-7edf6ef-20260728T174000Z/`.

The source, environment, immutable input, normalized intermediate, row-scaled
model, and inventory gates all passed. Row exponents were only `[-2,0]`; no
variables or constraints were removed or rescaled.

## Solver and replay

Mosek terminated after 82.99s with `SLOW_PROGRESS` and primal/dual
`FEASIBLE_POINT`. Original-model replay found:

| quantity | value | test |
|---|---:|---|
| normalized equality residual | `9.133500587283089e-11` | fail |
| normalized PSD violation | `1.1956643932013813e-17` | pass |
| objective direction | `1.0` | pass |
| normalized improving objective | `8.0241036275869e-6` | pass |
| ray scale | `124624.51214637807` | diagnostic |

The local replay reproduced equality/objective and rejection; its
below-tolerance PSD eigenvalue differed but also passed.

## Decision

Row scaling improved the unscaled normalized-recession residual from
`9.9912e-11` to `9.1335e-11`, only about nine percent and still about 91× the
acceptance tolerance. High-precision diagnosis remains dominated by rows 8,839
and 1 and the 17×17/271×271 blocks.

The original, exact-deduplicated, uniform-block, diagonal-congruence,
unit-improvement, and row-equilibrated formulations have now all been tested.
Further tolerance-preserving scaling variants are no longer decision-relevant.
The remaining strict-certificate route is exact affine projection followed by
a rigorous PSD argument; the present candidate has no demonstrated PSD margin
robust enough for that correction. The implementation focus therefore moves
to the Square source-gated status/audit boundary.
